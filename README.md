# SONiCS - Stutter mONte Carlo Simulation

## SUMMARY

SONiCS performs Monte Carlo simulations of the PCR of Short Tandem Repeats, calculates the likelihood of generating given support read out (reads per allele) out of the PCR pool and determines the most probable genotype based on the log likelihood distributions.

## REQUIREMENTS

* Python version >= 3.4 
* Cython
* Python modules: 
  * numpy, 
  * pandas, 
  * scipy, 
  * pymc

## INSTALLATION

In order to configure the script run the install.sh script. If more than one python3 installed supply the path for the python version. Script will check for all dependencies and inform you if some are missing. After installation it will run a simple test and report on whether everything is working as expected.

```bash
bash install.sh [PATH_TO_PYTHON3]
```

## RUNNING SONiCS

```bash
#example run
python3 sonics "8|9;9|20;10|24;12|1"
```

## DESCRIPTION

SONiCS - Stutter mONte Carlo Simulation Monte Carlo simulation of PCR based sequencing of Short Tandem Repeats with one or two alleles as a starting condition. Each simulation represents independent PCR reaction. Parameters controlling the PCR like amplification efficiency, parameters describing slippage up and down and capture step are randomly selected before each simulation. Initial alleles are chosen randomly or half randomly (in this mode only one allele is random, the other one will always be allele with maximum support). For each simulation sonics assigns log likelihood value describing the probability of generating the input readout from the PCR pool created in a given simulation. When certain number of simulations are run sonics compares the distributions of the log likelihood between the allele with highest median and all others. If adjusted highest p-value is below the threshold and median likelihood of the best genotype is at least 10 times greater than second best sonics calls genotype. As a result, the program reports the r squared and fraction of properly identified reads parameters from the best simulation, median log likelihood, log likelihood ratio between top two genotypes, highest adjusted p_value from Mann-Whitney U test, number of run repeats and called genotype.

### INPUT

SONiCS accepts either the string genotype having alleles, marked as number of repeats per molecule, and number of supporting them reads in the following format: allele1|reads;allele2|reads;allele3|reads or a **VCF file** with one genotype per line and genotypes for samples in the consecutive columns starting with 9th.

Examples:

* string genotype: "5|1;6|1;7|5;8|11;9|20;10|24;11|2;12|1"
* VCF file (lorem_ipsum marks additional information that can be in the file but won't be used)

```bash
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1    sample2    
genotype1 1001    .       GTGTGTGTGTGTGTGTGTGT    GTGTGTGTGTGTGTGTGTGTGT  0       .       END=1020;MOTIF=GT;LOREM_IPSUM=1;REF=10;LOREM_IPSUM=9;LOREM_IPSUM=11 GT:ALLREADS:LOREM_IPSUM    1/1:-2|1;0|54;2|914;4|15:lorem_ipsum   0/0:-6|2;-4|22;-2|150;0|1215;2|11:lorem_ipsum  
```

### MODES

#### Choosing starting alleles
SONiCS can choose input alleles randomly or half randomly - one allele is always the one with the highest support. By default sonics run in half_random mode, that can be changed to random by specifying **--random**.

#### Support for partial alleles
* **-t, --strict N** - Procedure when encountered partial repetitions of the motif while parsing the VCF file. Options: 0 - pick on random one of the closest alleles, 1 - exclude given STR, 2 -exclude given genotype. Default: 1

*Example:*
REF=11, MOTIF=GT, GT=0|54;1|27;2|914;4|15
* 0: 11|81;12|914;13|15 or 11|54;12|941;13|15
* 1: 11|54;12|914;13|15
* 2: exclude genotype from simulation pool

### PARAMETERS FOR SIMULATIONS

* **-r, --repetitions N** - maximum number of repetitions per genotype, default: 1000
* **PCR_CYCLES** - **required** parameter, number of PCR cycles before introducing the capture step
* **-c, --after_capture N** - number of PCR cycles after the introduction of the capture step; default: 12
* **-p, --padjusted N** - threshold for Bonferroni corrected p-value of Mann Whitney test of the likelihood distributions; default: 0.01 

### PARAMETERS SPECIFIC FOR PCR

Parameters are chosen on random from the given open interval before each simulation. 

* **-e, --efficiency N N** - PCR efficiency default: (0.001, 0.1)
* **-d, --down N N** - probability of slippage down - synthesizing molecule with less repeats; default: (0, 0.1)
* **-u, --up N N** - probability of slippage up - synthesizing molecule with more repeats; default: (0, 0.1)
* **-k, --capture N N** - parameter describing efficiency of the capture step. Values below zero favor capturing shorter molecules; default: (-0.25, 0.25)
* **--up_preference** - up-stutter doesn't have to be more probable than down-stutter; default slippage_down > slippage_up
* **-f, --floor N** - number that will be subtracted from the lowest number of STRs in the input genotype to set the threshold for the minimum number of STRs in a molecule for it to be included in the simulations. Default: 5

*Formula:*
min_strs = max(1, min(alleles_in_input) - floor)


### ALL OPTIONS

```
usage: sonics [-h] [-o OUT_PATH] [-n FILE_NAME] [-p PROCESSES] [-t N]
              [-c AFTER_CAPTURE] [-b BLOCK] [-a NAME] [-r REPS] [-g PADJUST]
              [-i LOGLIKE] [-s START_COPIES] [-l NOISE_THRESHOLD] [-f FLOOR]
              [-e MIN MAX] [-d MIN MAX] [-u MIN MAX] [-k MIN MAX] [-m] [-x]
              [-z] [-v] [--version]
              PCR_CYCLES INPUT

SONiCS - Stutter mONte Carlo Simulation Monte Carlo simulation of PCR based
sequencing of Short Tandem Repeats with one or two alleles as a starting
condition. Each simulation represents independent PCR reaction. Parameters
controlling the PCR like amplification efficiency, parameters describing
slippage up and down and capture step are randomly selected before each
simulation. Initial alleles are chosen randomly or half randomly (in this mode
only one allele is random, the other one will always be allele with maximum
support). For each simulation sonics assigns log likelihood value describing
the probability of generating the input readout from the PCR pool created in a
given simulation. When certain number of simulations are run sonics compares
the distributions of the log likelihood between the allele with highest median
and all others. If adjusted highest p-value is below the threshold and median
likelihood of the best genotype is at least e^8 times greater than second best
sonics calls genotype. As a result, the program reports the r squared and
fraction of properly identified reads parameters from the best simulation,
median log likelihood, log likelihood ratio between top two genotypes, highest
adjusted p_value from Mann-Whitney U test, number of run repeats and called
genotype.

positional arguments:
  PCR_CYCLES            Number of PCR cycles before introducing capture step.
  INPUT                 Either:allele composition - number of reads per
                        allele, example: all1|reads;all2|reads or path to a
                        VCF file if run with --vcf_mode option

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_PATH, --out_path OUT_PATH
                        Directory where output files will be stored. Warning:
                        SONiCS will overwrite results from previous runs
                        stored in the same directory if the file name won't be
                        changed. Default: .
  -n FILE_NAME, --file_name FILE_NAME
                        Output file name. Default: sonics_out.txt
  -p PROCESSES, --processes PROCESSES
                        Number of sub-processes used in multiprocessing mode.
                        It can be understood as additional processes lunched
                        by the main sonics process and 0 means that SONiCS
                        won't lunch additional processes. This option is valid
                        only in VCF mode. Default: 0
  -t N, --strict N      What to do when encountered partial repetitions of the
                        motif while parsing the VCF file. Options: 0 - pick on
                        random one of the closest alleles, 1 - exclude given
                        STR, 2 - exclude given genotype. Default: 1
  -c AFTER_CAPTURE, --after_capture AFTER_CAPTURE
                        How many cycles of PCR amplification were performed
                        after introducing capture step. Default: 12
  -b BLOCK, --block BLOCK
                        Block name, valid only with string genotype as input.
                        Default: Block
  -a NAME, --name NAME  Sample name, valid only with string genotype as input.
                        Default: sample
  -r REPS, --repetitions REPS
                        Number of maximum repetitions in the simulations.
                        Default: 1000
  -g PADJUST, --padjust PADJUST
                        Adjusted p-value threshold for selecting the allele.
                        Default: 0.001
  -i LOGLIKE, --loglike LOGLIKE
                        Log likelihood ratio threshold. Default: 2.3
  -s START_COPIES, --start_copies START_COPIES
                        Number of start copies. Default: 3e5
  -l NOISE_THRESHOLD, --noise_threshold NOISE_THRESHOLD
                        How much noise should be added? Noise shouldn't
                        comprise more than the maximum fraction threshold of
                        the initial PCR pool (before first amplification).
                        Default: 0.05
  -f FLOOR, --floor FLOOR
                        Number that will be subtracted from the lowest number
                        of STRs in the input genotype to set the threshold for
                        the minimum number of STRs in a molecule for it to be
                        included in the simulations. If set to -1 minimum
                        number of STRs set to 1. Default: 5
  -e MIN MAX, --efficiency MIN MAX
                        PCR efficiency before-after per cycle, i.e.
                        probability of the amplification. Default: (0.001,
                        0.1)
  -d MIN MAX, --down MIN MAX
                        Per unit probability of down-slippage - generating a
                        molecule with less repeats (min and max). Default: (0,
                        0.1)
  -u MIN MAX, --up MIN MAX
                        Per unit probability of up-slippage - generating a
                        molecule with more repeats (min and max). Default: (0,
                        0.1)
  -k MIN MAX, --capture MIN MAX
                        Capture parameter (min and max). Values below zero
                        favor capturing short alleles. Default: (-0.25, 0.25)
  -m, --vcf_mode        VCF file provided. Assuming that different samples, if
                        more than one is present, are put in the consecutive
                        columns, starting with 10th. Default: string mode.
  -x, --random          Randomly select alleles for simulations from the
                        input. Default: randomly select only one of the
                        alleles, the other is chosen based on the support
                        information from input.
  -z, --up_preference   Up-stutter doesn't have to be less probable than down
                        stutter. Default: probability of down-stutter must be
                        higher than the probability of up-stutter.
  -v, --verbose         Verbose mode.
  --version             show program's version number and exit
```