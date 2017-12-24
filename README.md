#SONiCS - Stutter mONte Carlo Simulation

##SUMMARY

SONiCS performs Monte Carlo simulation of the PCR of Short Tandem Repeats, calculates the likelihood of generating given support read out (reads per allele) out of the PCR pool and based on the log likelihood distributions determines the most probable genotype.

##REQUIREMENTS

* Python 3.6 >= #only tested on 3.6, need to test on 3.4
* Cython
* Python modules: logging, numpy, pandas, scipy, pymc, argparse, shutil

##INSTALLATION

In order to configure the script run the install.sh script. If more than one python3 installed supply the path for the python version. Script will check for all dependencies. 

```bash
bash ./install.sh [PATH_TO_PYTHON3]
```


##DESCRIPTION

```bash
usage: SONiCS [-h] [-o out_path] [-n FILE_NAME] [-c after_capture]
              [-r repetitions] [-p pvalue_threshold] [-s start_copies]
              [-f floor] [-a MIN MAX] [-e MIN MAX] [-d MIN MAX] [-u MIN MAX]
              [-k MIN MAX] [-m] [-i] [-x] [-y] [-z] [-t] [-v] [--version]
              PCR_CYCLES INPUT

Monte Carlo simulation of PCR based sequencing of Short Tandem Repeats with
one or two alleles as a starting condition.

positional arguments:
  PCR_CYCLES            Number of PCR cycles before capture.
  INPUT                 Either: allele composition - number of reads per
                        allele, example: allele1|reads;allele2|reads or path
                        to VCF file if run with --vcf_mode option

optional arguments:
  -h, --help            show this help message and exit
  -o out_path, --out_path out_path
                        Directory where output files will be stored. Warning:
                        SONiCS appends results to a file, so if an output file
                        exits it will add results to it, rather than
                        overwriting it.
  -n FILE_NAME, --file_name FILE_NAME
                        Output file name. Default: sonics_out
  -c after_capture, --after_capture after_capture
                        How many cycles of PCR amplification were performed
                        after introducing capture step. [12]
  -r repetitions, --repetitions repetitions
                        Number of maximum repetitions in the simulations
                        [1000]
  -p pvalue_threshold, --pvalue_threshold pvalue_threshold
                        P-value threshold for selecting the allele [0.01]
  -s start_copies, --start_copies start_copies
                        Number of start copies. [3e5]
  -f floor, --floor floor
                        Minimum number of repetitions per STR, default: 5 less
                        of minimum number of repetitions in the starting
                        alleles.
  -a MIN MAX, --amplification MIN MAX
                        Amplification parameters (min and max). Values below
                        zero favor amplifying shorter molecules. [-0.1 0.1]
  -e MIN MAX, --efficiency MIN MAX
                        PCR efficiency before-after per cycle, i.e. proportion
                        of molecules taken into each cycle [0.001 0.1]
  -d MIN MAX, --down MIN MAX
                        Per unit probability of slippage down (min and max).
                        [0 0.1]
  -u MIN MAX, --up MIN MAX
                        Per unit probability of slippage up (min and max). [0
                        0.1]
  -k MIN MAX, --capture MIN MAX
                        Capture parameter (min and max). Values below zero
                        favor capturing short alleles. [-0.25]
  -m, --vcf_mode        VCF file provided. Assuming that different samples in
                        consecutive columns, starting with 10th.
  -i, --save_intermediate
                        Save intermediate PCR pools when run in vcf mode to
                        save up time on simulating same initial conditions.
  -x, --random          Randomly select alleles for simulations. [Select the
                        most frequent allele]
  -y, --half_random     Randomly select second allele. [Select the second most
                        frequent allele]
  -z, --up_preference   Up-stutter doesn't have to be less than down
  -t, --strict          If input alleles include partial repetitions of the
                        motif, exclude given STR
  -v, --verbose         Verbose mode.
  --version             show program's version number and exit
```