#!/usr/bin/env python3

"""SONiCS - Stutter mONte Carlo Simulation
"""
desc = """Monte Carlo simulation of PCR based sequencing of Short Tandem Repeats 
with one or two alleles as a starting condition. Uses least squares 
as a method to asses the agreement between with the provided genotype 
and simulated genotypes.
"""

import argparse
import logging
import time
import re
import sonics

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
                continue
            tmp = line.split("\t")
            name = tmp[0]   
            ref = int(regexp_ref.search(tmp[7]).group(1))
            motif_length = len(regexp_motif.search(tmp[7]).group(1))
            gt = tmp[9].split(":")[1] 
            #sanity check
            diff = [int(f.split("|")[0]) % motif_length for f in gt.split(";")]
            if not all(i == 0 for i in diff):
                #raise ValueError("Partial repetition of the motif in %s" %(name))
                if strict:
                    logging.info("WARNING: Partial repetition of the motif in %s. Excluding block." %(name))
                    pass
                else:
                    logging.info("WARNING: Partial repetition of the motif in %s. Excluding partial allele(s)." %(name))
                    gt_tmp = [[int(int(f.split("|")[0]) / motif_length + ref), f.split("|")[1]] for f in gt.split(";") if int(f.split("|")[0]) % motif_length == 0]
                    genot = ";".join(["{}|{}".format(tup[0], tup[1]) for tup in gt_tmp])
                    genotypes_list.append((name, genot))
            else:
                gt_tmp = [[int(int(f.split("|")[0]) / motif_length + ref), f.split("|")[1]] for f in gt.split(";")]
                genot = ";".join(["{}|{}".format(tup[0], tup[1]) for tup in gt_tmp])
                genotypes_list.append((name, genot))
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
        "GENOTYPE",
        type=str,
        help="Either: allele composition - number of reads per allele. Example: allele1|reads;allele2|reads or path to VCF file if run with --vcf option."
    )

    ### OPTIONS
    parser.add_argument(
        "-m", "--vcf_mode",
        action='store_true',
        help="VCF file provided instead of single genotype."
    )
    parser.add_argument(
        "-o", "--out_file",
        type=str,
        help="If provided, sonics output will be written out to the OUTPUT_FILE",
        metavar="OUTPUT_FILE"
    )
    parser.add_argument(
        "-n", "--run_name",
        type=str,
        nargs=1,
        default="Block",
        help="Name of the run for single genotype run, default: Block.",
        metavar="run_name"
    )

    parser.add_argument(
        "-c", "--after_capture",
        type=int,
        default=12,
        help=""""How many cycles of PCR amplification was performed after introducing capture step. [12]"""
    )
    parser.add_argument(
        "-r", "--repetitions",
        default=1000,
        type=int,
        help="Number of repetitions in the simulations [1000]"
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
        help="Number of start copies. [3e5]"
    )
    parser.add_argument(
        "-f", "--floor",
        default=5,
        type=int,
        help="copy floor parameter"
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
        "-l", "--length",
        nargs=2,
        action='append',
        metavar=('MIN', 'MAX'),
        type=float,
        default=0,
        help="Length penalty (min and max) [0 0.1]"
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
        version='%(prog)s 0.1'
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
        'start_copies': args.start_copies,
        'floor': args.floor,
        'length': args.length,
        'genotype': args.GENOTYPE,
        'pvalue_threshold': args.pvalue_threshold
    }

    ranges = (
        args.down,
        args.up,
        args.capture,
        args.amplification,
        args.efficiency
    )

    if args.vcf_mode:
        logging.info("VCF file mode on.")
        parser_out = parse_vcf(args.GENOTYPE, strict=args.strict)
        if type(parser_out) == str:
            raise Exception("There is something wrong with VCF file, record: ".format(parser_out))
        else:
            for sample, genotype in parser_out:
                constants['alleles'], constants['max_allele'] = sonics.get_alleles(genotype)
                constants['genotype_total'] = sum(constants['alleles'])
                logging.info("Initiating simulation for {}".format(sample))
                start = time.time()
                result = sonics.monte_carlo(args.repetitions, constants, ranges)
                elapsed = time.time() - start
                logging.info(result)
                if args.out_file != None:
                    logging.debug(args.out_file)
                    with open(args.out_file, "a+") as of:
                        of.write("{}\t{}\n".format(sample, result)) #test if not slow, if too slow switch to below
                        #of.write(run_name)
                        #of.write("\t")
                        #of.write(result)
                        #of.write("\n") 
                logging.debug("Monte Carlo simulation took: %f second(s)" %elapsed)

    else:
        constants['alleles'], constants['max_allele'] = sonics.get_alleles(args.GENOTYPE)
        constants['genotype_total'] = sum(constants['alleles'])

        logging.info("Initiating simulation")
        start = time.time()
        result = sonics.monte_carlo(args.repetitions, constants, ranges)
        elapsed = time.time() - start
        logging.info(result)
        if args.out_file != None:
            with open(args.out_file, "a+") as of:
                of.write("{}\t{}\n".format(args.run_name, result))
        logging.debug("Monte Carlo simulation took: %f second(s)" %elapsed)


if __name__ == '__main__':
    main()
