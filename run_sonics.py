#!PYTHON_PATH_TO_FILL

"""SONiCS - Stutter mONte Carlo Simulation
"""
desc = """Monte Carlo simulation of PCR based sequencing of Short Tandem Repeats 
with one or two alleles as a starting condition.
"""
import re
import os
import argparse
import logging
import time
import sonics
import shutil
import numpy as np

def run_sonics(feed_in, constants, ranges, strict, repetitions, out_path, file_name="sonics_out", vcf_mode=False, save_intermediate=False, name="sample", block="Block"):
    """
    save_intermediate => create tmp directory with files for each starting conditions. Each file would have one PCR pool per line. 
    """
    if save_intermediate:
        intermediate = os.path.join(out_path, "tmp")
        try:
            os.stat(intermediate)
        except:
            os.mkdir(intermediate)  
    else:
        intermediate = None

    if vcf_mode:
        parser_out = parse_vcf(feed_in, strict)
        if type(parser_out) == str:
            raise Exception("There is something wrong with VCF file: {} record: {}".format(vcf_file, parser_out))
        else:
            for name, block, genotype in parser_out:
                constants['alleles'], constants['max_allele'] = sonics.get_alleles(genotype)
                constants['genotype_total'] = sum(constants['alleles'])
                if constants['genotype_total'] == 0:
                    logging.warning("Less than 2 alleles provided in input, skipping this genotype: {} for the sample: {}.".format(block, name))
                    continue
                logging.info("Initiating simulation for {} {}".format(name, block))
                start = time.time()
                result = sonics.monte_carlo(repetitions, constants, ranges, intermediate, block, name)
                elapsed = time.time() - start
                logging.info(result)
                op = os.path.join(out_path, file_name)
                with open(op, "a+") as of:
                    of.write("{}\t{}\t{}\n".format(name, block, result))
                logging.debug("Monte Carlo simulation took: %f second(s)" %elapsed)
    else:
        genotype = feed_in
        constants['alleles'], constants['max_allele'] = sonics.get_alleles(genotype)
        constants['genotype_total'] = sum(constants['alleles'])
        if constants['genotype_total'] == 0:
            raise Exception("Less than 2 alleles provided in input, skipping this genotype: {} for the sample: {}.".format(genotype, name))
        logging.info("Initiating simulation")
        start = time.time()
        result = sonics.monte_carlo(repetitions, constants, ranges)
        elapsed = time.time() - start
        logging.info(result)
        op = os.path.join(out_path, file_name)
        with open(op, "a+") as of:
            of.write("{}\t{}\t{}\n".format(name, block, result))
        logging.debug("Monte Carlo simulation took: %f second(s)" %elapsed)


def parse_vcf(file_path, strict=1):
    """
    Transforms genotype from lobSTR VCF file (i.e. allele1:#;allele2:#;allele3:#, where allele is presented as number of nucleotides difference from reference). To do so, the function first extracts number of repetitions in the reference, length of the reference and genotype to be transformed and then each value from allele field is dived by the length of the motif and the reference is added.

    Example of transformation
    
    in: REF=10, MOTIF=TCTA, GT=-4|1;0|48
    out: GT=9|1;0:48.

    If support for partial alleles encountered strict argument determines the procedure. 0 - support from partial alleles is transfered to one of the two closest alleles; 1 - partial allele is excluded; 2 - particular genotype is excluded from simulation.
    """        
    genotypes_list=[]
    regexp_ref = re.compile(".*REF=([0-9.]*);.*")
    regexp_motif = re.compile(".*MOTIF=([TCGA]*);.*")
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    samples = line.strip("\n").split("\t")[9:]
                continue
            tmp = line.split("\t")
            block_name = tmp[0]   
            ref = int(regexp_ref.search(tmp[7]).group(1).split(".")[0])
            motif_length = len(regexp_motif.search(tmp[7]).group(1))
            all_reads_loc = tmp[8].split(":").index("ALLREADS")
            n = 9 #sample index
            for sample in samples:   
                gt_string = tmp[n].split(":")[all_reads_loc] 
                gt = gt_string.split(";")
                alleles = np.array([int(f.split("|")[0]) / motif_length + ref for f in gt])
                support = np.array([f.split("|")[1] for f in gt], dtype = int)
                diff = alleles % 1 
                if sum(diff != 0) > 0:
                    if strict == 2:
                        logging.warning("Partial repetition of the motif in %s sample: %s. Excluding block." %(block_name, sample))
                        continue
                    elif strict == 0:
                        logging.warning("Partial repetition of the motif in %s sample: %s. Adding support from partial allele to one of the closest allele on random." %(block_name, sample))
                        for partial in alleles[diff != 0]:
                            choice = np.random.choice([0,1])
                            support[alleles == int(partial) + choice] += support[alleles == partial]
                    else:
                        logging.warning("Partial repetition of the motif in %s sample: %s. Excluding partial allele(s)." %(block_name, sample))
                genot = ";".join(["{}|{}".format(int(a), s) for a,s in zip(alleles[diff == 0], support[diff == 0])])
                genotypes_list.append((sample, block_name, genot))
                n+=1
    return genotypes_list

def main():
    parser = argparse.ArgumentParser(
        prog="SONiCS",
        description=desc
    )
    parser.add_argument(
        "PCR_CYCLES",
        type=int,
        help="Number of PCR cycles before introducing capture step."
    )
    parser.add_argument(
        "INPUT",
        type=str,
        help="""Either: 
        allele composition - number of reads per allele, example: allele1|reads;allele2|reads 
        or path to a VCF file if run with --vcf_mode option"""
    )

    ### OPTIONS
    parser.add_argument(
        "-o", "--out_path",
        type=str,
        default=".",
        help="Directory where output files will be stored. Warning: SONiCS appends results to a file, so if an output file exits it will add results to it, rather than overwriting it.",
        metavar="out_path"
    )
    parser.add_argument(
        "-n", "--file_name",
        type=str,
        nargs=1,
        default="sonics_out.txt",
        help="Output file name. Default: sonics_out.txt",
        metavar="FILE_NAME"
    )
    parser.add_argument(
        "-t", "--strict",
        metavar="N",
        default=0,
        help="""Procedure when encountered partial repetitions of the motif while parsing the VCF file. Default: 1
        * 0 pick on random one of the closest alleles
        * 1 exclude given STR
        * 2 exclude given genotype

        Example: REF=11, MOTIF=GT, GT=0|54;1|27;2|914;4|15

        * 0: 11|81;12|914;13|15 or 11|54;12|941;13|15
        * 1: 11|54;12|914;13|15
        * 2: exclude genotype from simulation pool
        """
    )
    parser.add_argument(
        "-c", "--after_capture",
        type=int,
        default=12,
        metavar="after_capture",
        help="How many cycles of PCR amplification were performed after introducing capture step. Default: 12"
    )
    parser.add_argument(
        "-r", "--repetitions",
        default=1000,
        type=int,
        metavar="REPS",
        help="Number of maximum repetitions in the simulations Default: 1000"
    )
    parser.add_argument(
        "-p", "--pvalue_threshold",
        nargs=1,
        metavar="PVALUE",
        type=float,
        default=0.01,
        help="P-value threshold for selecting the allele Default: 0.01"
    )
    parser.add_argument(
        "-s", "--start_copies",
        default=3e5,
        type=int,
        metavar="START_COPIES",
        help="Number of start copies. Default: 3e5"
    )
    parser.add_argument(
        "-f", "--floor",
        type=int,
        default=5,
        metavar="floor",
        help="""Number that will be subtracted from the lowest number of STRs in the input genotype to set the threshold for the minimum number of STRs in a molecule for it to be included in the simulations.

        Formula:
            min_strs = max(1, min(alleles_in_input) - floor)"""
    )
    parser.add_argument(
        "-e", "--efficiency",
        nargs=2,
        action='append',
        default=(0.001, 0.1),
        metavar=('MIN', 'MAX'),
        type=float,
        help="PCR efficiency before-after per cycle, i.e. probability of amplification [(0.001, 0.1)]"
    )
    parser.add_argument(
        "-d", "--down",
        nargs=2,
        action='append',
        metavar=('MIN', 'MAX'),
        type=float,
        default=(0, 0.1),
        help="Per unit probability of down-slippage - generating a molecule with less repeats (min and max). [(0, 0.1)]"
    )
    parser.add_argument(
        "-u", "--up",
        nargs=2,
        action='append',
        metavar=('MIN', 'MAX'),
        type=float,
        default=(0, 0.1),
        help="Per unit probability of up-slippage - generating a molecule with more repeats (min and max). [(0, 0.1)]"
    )
    parser.add_argument(
        "-k", "--capture",
        nargs=2,
        action='append',
        default=(-0.25, 0.25),
        type=float,
        metavar=('MIN', 'MAX'),
        help="Capture parameter (min and max). Values below zero favor capturing short alleles. [(-0.25, 0.25)]"
    )
    parser.add_argument(
        "-m", "--vcf_mode",
        action='store_true',
        default=False,
        help="VCF file provided. Assuming that different samples, if more than one is present, are put in the consecutive columns, starting with 10th."
    )
    parser.add_argument(
        "-i", "--save_intermediate",
        default=False,
        action="store_true",
        help="Save intermediate PCR pools when run in vcf mode to save up time on simulating same initial conditions. [not saing intermediate results and simulating PCR for each run from scratch]")
    parser.add_argument(
        "-x", "--random",
        help="Randomly select alleles for simulations. [Select both alleles based on the support information from the input genotype]",
        action='store_true'
    )
    parser.add_argument(
        "-y", "--half_random",
        action='store_true',
        help="Randomly select only the second allele, the first is the most supported one. [Select both alleles based on the support information from the input genotype]"
    )
    parser.add_argument(
        "-z", "--up_preference",
        action='store_true',
        help="Up-stutter doesn't have to be less probable than down stutter. [probability of down-stutter higher than the probability of up-stutter]"
    )
    parser.add_argument(
        "-v", "--verbose",
        help="Verbose mode.",
        action='store_true'
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 0.0.1'
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

    capture_cycle = args.PCR_CYCLES + 1 # first cycle of reamplification after capture is the capture step
    n_cycles = args.PCR_CYCLES + args.after_capture

    if (args.half_random) & (args.random):
        raise Exception("It is not possible to execute the modes random and half_random at the same time.")

    constants = {
        'random': args.random,
        'half_random': args.half_random,
        'n_cycles': n_cycles,
        'capture_cycle': capture_cycle,
        'up_preference': args.up_preference,
        'floor': args.floor,
        'start_copies': args.start_copies,
        'genotype': args.INPUT,
        'pvalue_threshold': args.pvalue_threshold
    }

    ranges = (
        args.down,
        args.up,
        args.capture,
        args.efficiency
    )

    run_sonics(
        feed_in = args.INPUT, 
        constants = constants, 
        ranges = ranges, 
        strict = args.strict, 
        repetitions = args.repetitions, 
        out_path = args.out_path, 
        file_name = args.file_name, 
        vcf_mode = args.vcf_mode, 
        save_intermediate = args.save_intermediate
    )

    if args.save_intermediate:
       shutil.rmtree(os.path.join(args.out_path, "tmp")) 

if __name__ == '__main__':
    main()
