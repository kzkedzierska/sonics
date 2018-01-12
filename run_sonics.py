import re
import os
import argparse
import logging
import time
import shutil
from itertools import repeat
from multiprocessing import Pool
import numpy as np
import sonics

DESC = """SONiCS - Stutter mONte Carlo Simulation
Monte Carlo simulation of PCR based sequencing of Short Tandem
Repeats with one or two alleles as a starting condition.
"""

def run_sonics(feed_in, constants, ranges, sonics_run_options):
    """
    save_intermediate => create tmp directory with files for each starting
    conditions. Each file would have one PCR pool per line.
    """
    out_path = sonics_run_options['out_path']
    file_name = sonics_run_options['file_name']

    try:
        os.stat(out_path)
    except FileNotFoundError:
        os.mkdir(out_path)

    if sonics_run_options['save_intermediate']:
        intermediate = os.path.join(out_path, "tmp")
        try:
            os.stat(intermediate)
        except FileNotFoundError:
            os.mkdir(intermediate)
    else:
        intermediate = None

    #output file path
    out_file_path = os.path.join(out_path, file_name)

    #options for processing one genotype
    options = (sonics_run_options['repetitions'],
               intermediate,
               sonics_run_options['verbose'],
               out_file_path)

    if sonics_run_options['vcf_mode']:
        #delete output file if exists
        try:
            os.remove(out_file_path)
        except OSError:
            pass

        parser_out = parse_vcf(feed_in, sonics_run_options['strict'])

        if not parser_out:
            raise Exception("No valid genotype after parsing VCF file.")

        if sonics_run_options['processes'] > 0:
            pool = Pool(sonics_run_options['processes'])
            pool.starmap(process_one_genotype, zip(parser_out,
                                                   repeat(constants),
                                                   repeat(ranges),
                                                   repeat(options)))
        else:

            print(parser_out)
            list(map(process_one_genotype,
                     parser_out,
                     repeat(constants),
                     repeat(ranges),
                     repeat(options)))

    else:
        genotype = feed_in
        input_tuple = (
            sonics_run_options['name'],
            sonics_run_options['block'],
            genotype
        )

        process_one_genotype(
            input_tuple,
            constants,
            ranges,
            options
        )

def parse_vcf(file_path, strict=1):
    """
    Transforms genotype from lobSTR VCF file (i.e. all1:#;all2:#;all3:#,
    where allele is presented as number of nucleotides difference from
    reference). To do so, the function first extracts number of
    repetitions in the reference, length of the reference and genotype
    to be transformed and then each value from allele field is dived
    by the length of the motif and the reference is added.

    Example of transformation

    in: REF=10, MOTIF=TCTA, GT=-4|1;0|48
    out: GT=9|1;0:48.

    If support for partial alleles encountered strict argument
    determines the procedure. 0 - support from partial alleles is
    transfered to one of the two closest alleles; 1 - partial allele
    is excluded; 2 - particular genotype is excluded from simulation.
    """
    genotypes_list = []
    regexp_ref = re.compile(".*REF=([0-9.]*);.*")
    regexp_motif = re.compile(".*MOTIF=([TCGA]*);.*")
    with open(file_path, "r") as file:
        for line in file:
            sample_ind = 9 #sample index
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    samples = line.strip("\n").split("\t")[sample_ind:]
                continue
            line_list = line.split("\t")
            block_name = line_list[0]
            info_field = line_list[7]
            try:
                ref = int(regexp_ref.search(info_field).group(1).split(".")[0])
            except AttributeError:
                raise Exception(("There's something wrong with the vcf file, "
                                 "expecting REF in the INFO field "
                                 "in the 8th column."))
            try:
                motif_length = len(regexp_motif.search(info_field).group(1))
            except AttributeError:
                raise Exception(("There's something wrong with the vcf file, "
                                 "expecting MOTIF in the INFO field "
                                 "in the 8th column."))

            all_reads_loc = line_list[8].split(":").index("ALLREADS")

            for sample in samples:
                gt_string = line_list[sample_ind].split(":")[all_reads_loc]
                gt_list = gt_string.split(";")
                alleles_list = [
                    int(f.split("|")[0]) / motif_length + ref for f in gt_list
                ]
                alleles = np.array(alleles_list)
                support = np.array([f.split("|")[1] for f in gt_list],
                                   dtype=int)
                diff = alleles % 1
                if sum(diff != 0) > 0:
                    if strict == 2:
                        logging.warning(("Partial repetition of the motif in "
                                         "block: %s sample: %s. "
                                         "Excluding."), block_name, sample)
                        continue
                    elif strict == 0:
                        logging.warning(("Partial repetition of the motif in "
                                         "block: %s, sample: %s. Adding "
                                         "support from partial allele to one "
                                         "of the closest allele on random."),
                                        block_name, sample)

                        for partial in alleles[diff != 0]:
                            choice = np.random.choice([0, 1])
                            support[int(partial) + choice] += support[partial]
                            support[partial] = 0
                    else:
                        logging.warning(("Partial repetition of the motif in "
                                         "block: %s sample: %s. Excluding "
                                         "allele."), block_name, sample)
                genot_list = [
                    "%i|%i" %(a, s) for a, s in zip(alleles, support) if s > 0
                ]
                genot = ";".join(genot_list)
                genotypes_list.append((sample, block_name, genot))
                sample_ind += 1
    return genotypes_list

def process_one_genotype(input_tuple, constants, ranges, options):
    """runs one Monte Carlo simulation"""
    repetitions, intermediate, verbose, out_file_path = options
    name, block, genotype = input_tuple
    alleles, constants['max_allele'] = sonics.get_alleles(genotype)
    constants['genotype_total'] = sum(alleles)

    if constants['genotype_total'] == 0:
        logging.warning(("Less than 2 alleles provided in input,"
                         "skipping this genotype: %s"
                         "for sample: %s."), block, name)
        return
    #determine if can add noise
    frac = (np.amin(alleles[alleles.nonzero()])
            / sum(alleles)
            * np.count_nonzero(alleles))
    noise_coef = frac / (frac + 1)

    constants['alleles'] = alleles
    constants['noise_coef'] = noise_coef

    logging.info("Initiating simulation for %s %s", name, block)

    all_simulation_params = {
        'intermediate': intermediate,
        'verbose': verbose,
        'block': block,
        'name': name
    }

    start = time.time()
    result = sonics.monte_carlo(
        repetitions,
        constants,
        ranges,
        all_simulation_params
    )
    elapsed = time.time() - start

    logging.info(result)
    with open(out_file_path, "a+") as out_file:
        out_file.write("{}\t{}\t{}\n".format(name, block, result))

    logging.debug("Monte Carlo simulation took: %f second(s)", elapsed)

#TODO: Rewrite help.

def main():
    parser = argparse.ArgumentParser(
        prog="SONiCS",
        description=DESC
    )
    parser.add_argument(
        "PCR_CYCLES",
        type=int,
        help="Number of PCR cycles before introducing capture step."
    )
    parser.add_argument(
        "INPUT",
        type=str,
        help=("Either:"
              "allele composition - number of reads per allele,"
              "example: all1|reads;all2|reads or path to a VCF file"
              "if run with --vcf_mode option")
    )

    ### OPTIONS
    parser.add_argument(
        "-o", "--out_path",
        type=str,
        default=".",
        help=("Directory where output files will be stored. Warning: if file"
              "with the given name exists, SONiCS will overwrite it. Default: ."),
        metavar="OUT_PATH"
    )
    parser.add_argument(
        "-n", "--file_name",
        type=str,
        default="sonics_out.txt",
        metavar="FILE_NAME",
        help="Output file name. Default: sonics_out.txt"
    )
    parser.add_argument(
        "-p", "--processes",
        type=int,
        default=0,
        metavar="PROCESSES",
        help=("Number of sub-processes used in multiprocessing mode. "
              "It can be understood as additional processes lunched by "
              "the main sonics process and 0 means that SONiCS won't "
              "lunch additional processes. "
              "This option is valid only in VCF mode. Default: 0"))
    parser.add_argument(
        "-t", "--strict",
        type=int,
        default=1,
        metavar="N",
        help=("Procedure when encountered partial repetitions of the motif"
              "while parsing the VCF file. Options: 0 - pick on random one of"
              "the closest alleles, 1 - exclude given STR, 2 -exclude given"
              "genotype. Default: 1")
    )
    parser.add_argument(
        "-c", "--after_capture",
        type=int,
        default=12,
        metavar="AFTER_CAPTURE",
        help=("How many cycles of PCR amplification were performed after"
              "introducing capture step. Default: 12")
    )
    parser.add_argument(
        "-b", "--block",
        type=str,
        default="Block",
        metavar="BLOCK",
        help=("Block name, valid only with string genotype as input."
              "Default: Block")
    )
    parser.add_argument(
        "-a", "--name",
        type=str,
        default="sample",
        metavar="NAME",
        help=("Sample name, valid only with string genotype as input."
              "Default: sample")
    )
    parser.add_argument(
        "-r", "--repetitions",
        default=1000,
        type=int,
        metavar="REPS",
        help="Number of maximum repetitions in the simulations. Default: 1000"
    )
    parser.add_argument(
        "-g", "--pvalue_threshold",
        metavar="PVALUE",
        type=float,
        default=0.01,
        help="P-value threshold for selecting the allele. Default: 0.01"
    )
    parser.add_argument(
        "-s", "--start_copies",
        default=3e5,
        type=int,
        metavar="START_COPIES",
        help="Number of start copies. Default: 3e5"
    )
    parser.add_argument(
        "-l", "--noise_threshold", #noise ~ loud
        default=0.05,
        type=float,
        metavar="NOSIE_THRESHOLD",
        help=("What fraction of the staring pool should noise comprise."
              "Default: 0.05") #the worst sentence ever written in English
    )
    parser.add_argument(
        "-f", "--floor",
        type=int,
        default=5,
        metavar="floor",
        help=("Number that will be subtracted from the lowest number of STRs"
              "in the input genotype to set the threshold for the minimum"
              "number of STRs in a molecule for it to be included in the"
              "simulations. Default: 5")
    )
    parser.add_argument(
        "-e", "--efficiency",
        nargs=2,
        action='append',
        default=(0.001, 0.1),
        metavar=('MIN', 'MAX'),
        type=float,
        help=("PCR efficiency before-after per cycle, i.e. probability of the"
              "amplification. Default: (0.001, 0.1)")
    )
    parser.add_argument(
        "-d", "--down",
        nargs=2,
        action='append',
        metavar=('MIN', 'MAX'),
        type=float,
        default=(0, 0.1),
        help=("Per unit probability of down-slippage - generating a molecule"
              "with less repeats (min and max). Default: (0, 0.1)")
    )
    parser.add_argument(
        "-u", "--up",
        nargs=2,
        action='append',
        metavar=('MIN', 'MAX'),
        type=float,
        default=(0, 0.1),
        help=("Per unit probability of up-slippage - generating a molecule"
              "with more repeats (min and max). Default: (0, 0.1)")
    )
    parser.add_argument(
        "-k", "--capture",
        nargs=2,
        action='append',
        default=(-0.25, 0.25),
        type=float,
        metavar=('MIN', 'MAX'),
        help=("Capture parameter (min and max). Values below zero favor"
              "capturing short alleles. Default: (-0.25, 0.25)")
    )
    parser.add_argument(
        "-m", "--vcf_mode",
        action="store_true",
        default=False,
        help=("VCF file provided. Assuming that different samples, if more"
              "than one is present, are put in the consecutive columns,"
              "starting with 10th. Default: string mode.")
    )
    parser.add_argument(
        "-i", "--save_intermediate",
        default=False,
        action="store_true",
        help=("Save intermediate PCR pools when run in vcf mode to save up"
              "time on simulating same initial conditions."
              "Default: not saving intermediate results and simulating PCR"
              "for each run from scratch")
    )
    parser.add_argument(
        "-x", "--random",
        action="store_true",
        help=("Randomly select alleles for simulations from the input."
              "Default: selects both alleles based on the support information"
              "from the input.")
    )
    parser.add_argument(
        "-y", "--half_random",
        action="store_true",
        help=("Randomly select only the second allele, the first is the most"
              "supported one. Default: selects both alleles based on the"
              "support information from the input.")
    )
    parser.add_argument(
        "-z", "--up_preference",
        action="store_true",
        help=("Up-stutter doesn't have to be less probable than down stutter."
              "Default: probability of down-stutter higher than the"
              "probability of up-stutter.")
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose mode."
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 0.0.4'
    )

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(levelname)s %(message)s"
        )
    else:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)s %(message)s"
        )

    logging.debug('SONiCS run with the following options:')
    logging.debug(args)

    #first cycle of reamplification after capture is the capture step
    capture_cycle = args.PCR_CYCLES + 1
    n_cycles = args.PCR_CYCLES + args.after_capture

    if (args.half_random) & (args.random):
        raise Exception(("It is not possible to execute the modes random and"
                         "half_random at the same time."))

    constants = {
        'random': args.random,
        'half_random': args.half_random,
        'n_cycles': n_cycles,
        'capture_cycle': capture_cycle,
        'up_preference': args.up_preference,
        'floor': args.floor,
        'start_copies': args.start_copies,
        'genotype': args.INPUT,
        'pvalue_threshold': args.pvalue_threshold,
        'noise_threshold': args.noise_threshold
    }

    ranges = (
        args.down,
        args.up,
        args.capture,
        args.efficiency
    )

    verb = True if args.verbose else False

    sonics_run_options = {
        'repetitions': args.repetitions,
        'out_path': args.out_path,
        'strict': args.strict,
        'file_name': args.file_name,
        'vcf_mode': args.vcf_mode,
        'save_intermediate': args.save_intermediate,
        'name': args.name,
        'block': args.block,
        'processes': args.processes,
        'verbose': verb
    }

    run_sonics(
        feed_in=args.INPUT,
        constants=constants,
        ranges=ranges,
        sonics_run_options=sonics_run_options
    )

    if args.save_intermediate:
        shutil.rmtree(os.path.join(args.out_path, "tmp"))

if __name__ == '__main__':
    main()
