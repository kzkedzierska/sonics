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

def run_sonics(feed_in, constants, ranges, strict, repetitions, out_path, file_name="sonics_out", vcf_mode=False, save_intermediate=False, name="sample", block="Block"):
    """
    save_intermediate => create tmp directory with files for each starting conditions. Each file would have one PCR pool per line. 
    """
    if save_intermediate:
        intermediate = "/".join([out_path, "tmp"])
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


def parse_vcf(file_path, strict=False):
    """
    Transforms genotype from lobSTR VCF file (i.e. allele1:#;allele2:#;allele3:#, where allele is presented as number of nucleotides difference from reference). To do so, the function first extracts number of repetitions in the reference, length of the reference and genotype to be transformed and then each value from allele field is dived by the length of the motif and the reference is added.

    Example of transformation
    
    in: REF=10, MOTIF=TCTA, GT=-4|1;0|48
    out: GT=9|1;0:48.
    """        
    genotypes_list=[]
    regexp_ref = re.compile(".*REF=([0-9]*);.*")
    regexp_motif = re.compile(".*MOTIF=([TCGA]*);.*")
    with open(file_path, "r") as file:
        for line in file.readlines():
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    samples = line.strip("\n").split("\t")[9:]
                continue
            tmp = line.split("\t")
            block_name = tmp[0]   
            ref = int(regexp_ref.search(tmp[7]).group(1))
            motif_length = len(regexp_motif.search(tmp[7]).group(1))
            all_reads_loc = tmp[8].split(":").index("ALLREADS")
            n = 9 #sample index
            for sample in samples:   
                gt = tmp[n].split(":")[all_reads_loc] 
                #sanity check
                diff = [int(f.split("|")[0]) % motif_length for f in gt.split(";")]
                if not all(i == 0 for i in diff):
                    if strict:
                        logging.warning("Partial repetition of the motif in %s sample: %s. Excluding block." %(block_name, sample))
                        pass
                    else:
                        logging.warning("Partial repetition of the motif in %s sample: %s. Excluding partial allele(s)." %(block_name, sample))
                        gt_tmp = [[int(int(f.split("|")[0]) / motif_length + ref), f.split("|")[1]] for f in gt.split(";") if int(f.split("|")[0]) % motif_length == 0]
                        genot = ";".join(["{}|{}".format(tup[0], tup[1]) for tup in gt_tmp])
                        genotypes_list.append((sample, block_name, genot))
                else:
                    gt_tmp = [[int(int(f.split("|")[0]) / motif_length + ref), f.split("|")[1]] for f in gt.split(";")]
                    genot = ";".join(["{}|{}".format(tup[0], tup[1]) for tup in gt_tmp])
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
        help="Number of PCR cycles before capture."
    )
    parser.add_argument(
        "INPUT",
        type=str,
        help="""Either: 
        allele composition - number of reads per allele, example: allele1|reads;allele2|reads 
        or path to VCF file if run with --vcf_mode option"""
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
        "-c", "--after_capture",
        type=int,
        default=12,
        metavar="after_capture",
        help="How many cycles of PCR amplification were performed after introducing capture step. [12]"
    )
    parser.add_argument(
        "-r", "--repetitions",
        default=1000,
        type=int,
        metavar="repetitions",
        help="Number of maximum repetitions in the simulations [1000]"
    )
    parser.add_argument(
        "-p", "--pvalue_threshold",
        nargs=1,
        metavar="pvalue_threshold",
        type=float,
        default=0.01,
        help="P-value threshold for selecting the allele [0.01]"
    )
    parser.add_argument(
        "-s", "--start_copies",
        default=3e5,
        type=int,
        metavar="start_copies",
        help="Number of start copies. [3e5]"
    )
    parser.add_argument(
        "-f", "--floor",
        type=int,
        default=5,
        metavar="floor",
        help="Minimum number of repetitions per STR, default: 5 less of minimum number of repetitions in the starting alleles."
    )
    parser.add_argument(
        "-a", "--amplification",
        nargs=2,
        action='append',
        metavar=('MIN', 'MAX'),
        type=float,
        default=(-0.1, 0.1),
        help="Amplification parameters (min and max). Values below zero favor amplifying shorter molecules. [-0.1 0.1]"
    )
    parser.add_argument(
        "-e", "--efficiency",
        nargs=2,
        action='append',
        default=(0.001, 0.1),
        metavar=('MIN', 'MAX'),
        type=float,
        help="PCR efficiency before-after per cycle, i.e. proportion of molecules taken into each cycle [0.001 0.1]"
    )
    parser.add_argument(
        "-d", "--down",
        nargs=2,
        action='append',
        metavar=('MIN', 'MAX'),
        type=float,
        default=(0, 0.1),
        help="Per unit probability of slippage down (min and max). [0 0.1]"
    )
    parser.add_argument(
        "-u", "--up",
        nargs=2,
        action='append',
        metavar=('MIN', 'MAX'),
        type=float,
        default=(0, 0.1),
        help="Per unit probability of slippage up (min and max). [0 0.1]"
    )
    parser.add_argument(
        "-k", "--capture",
        nargs=2,
        action='append',
        default=(-0.25, 0.25),
        type=float,
        metavar=('MIN', 'MAX'),
        help="Capture parameter (min and max). Values below zero favor capturing short alleles. [-0.25]"
    )
    parser.add_argument(
        "-m", "--vcf_mode",
        action='store_true',
        default=False,
        help="VCF file provided. Assuming that different samples in consecutive columns, starting with 10th."
    )
    parser.add_argument(
        "-i", "--save_intermediate",
        default=False,
        action="store_true",
        help="Save intermediate PCR pools when run in vcf mode to save up time on simulating same initial conditions.")
    parser.add_argument(
        "-x", "--random",
        help="Randomly select alleles for simulations. [Select the most frequent allele]",
        action='store_true'
    )
    parser.add_argument(
        "-y", "--half_random",
        action='store_true',
        help="Randomly select second allele. [Select the second most frequent allele]"
    )
    parser.add_argument(
        "-z", "--up_preference",
        action='store_true',
        help="Up-stutter doesn't have to be less than down"
    )
    parser.add_argument(
        "-t", "--strict",
        action='store_true',
        default=False,
        help="If input alleles include partial repetitions of the motif, exclude given STR"
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
        args.amplification,
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
