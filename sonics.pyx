import logging
import os
from itertools import repeat, chain
import cython
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from pymc.distributions import multivariate_hypergeometric_like as mhl
cimport numpy as np
DTYPE = np.int
ctypedef np.int_t DTYPE_t

# inspired by Perl script stutterSim.pl by Logan Kistler

__author__ = "Katarzyna Kedzierska"
__email__ = "kzk5f@virginia.edu"



# cython: linetrace=True

# FUNCTIONS RUN ONLY ONCE/TWICE PER GENOTYPE
def get_alleles(genot_input):
    """get a dictionary with alleles as keys and number_of_reads as values
    from the genotype ('allele1|#;allele2|#').
    """
    max_allele = 500
    alleles = np.zeros(max_allele, dtype=DTYPE)
    for f in genot_input.split(';'):
        if f != "":
            pair = f.split("|")
            alleles[int(pair[0])] = int(pair[1])
    #alleles[alleles < max(alleles) * 0.01] = 0
    if len(alleles.nonzero()[0]) < 2:
        alleles = np.zeros(max_allele, dtype=DTYPE)

    #logging.info(alleles)
    return alleles, max_allele

def generate_params(r, pref):
    """generates random parameters from given ranges"""
    down, up, cap, eff = r
    small_number = 1e-16 # to make the range exclusive, instead of inclusive
    d = np.random.uniform(down[0] + small_number, down[1])
    if pref:
        u = np.random.uniform(up[0] + small_number, up[1])
    else:
        u = np.random.uniform(up[0] + small_number, d)
    c = np.random.uniform(cap[0] + small_number, cap[1])
    p = np.random.uniform(eff[0] + small_number, eff[1])  # pcr-efficiency
    # logging.debug({'down': d, 'up': u, 'capture': c, 'efficiency': p})
    return {'down': d, 'up': u, 'capture': c, 'efficiency': p}

def monte_carlo(max_n_reps, constants, ranges, intermediate=None,
                block="", name="", verbose=False):
    """Runs Monte Carlo simulation of the PCR amplification until
    p_value threshold for the Mann Whitney test is reached or
    the number of repetition reaches the maximum.

    Scheme:
    1) Run the first n repetitions.
    2) Calculate the highest p value for all the comparisons between
    set with the highest log likelihood median and others.
    3) Check if the highest p value with Bonferroni corrections is
    lower than the threshold. If true, stop. If false, go to step 1
    with new n now equal to 4*n if the it's the even round of
    repetitions or 2*n it's the odd. That way the sop will be made
    every 100, 500, 1000, 5000 etc. repetitions.

    Modification in intermediate mode:
    Instead of running the repetitions
    """
    pvalue_threshold = constants["pvalue_threshold"]
    """successful - parameter that helps distinguish between
    the simulations needing more repetitions and the ones
    that are beyond the abilities of SONiCS"""
    successful = True
    results = list()
    run_reps = 0
    reps_round = 1
    reps = 100 if max_n_reps > 100 else max_n_reps
    while run_reps < max_n_reps:

        results.extend(list(map(
            one_repeat,
            repeat("two_alleles", reps),
            repeat(constants),
            repeat(ranges),
            repeat(intermediate),
            repeat(reps)
        )))

        """even out the comparisons between homozygous and each
        of the heterozygous genotypes group by genotype"""
        results_pd_tmp = pd.DataFrame.from_records(results).groupby(3)
        #get the median of number of simulations per each genotype
        homo_reps = int(results_pd_tmp.size().median())

        results.extend(list(map(
            one_repeat,
            repeat("one_allele", homo_reps),
            repeat(constants),
            repeat(ranges),
            repeat(intermediate),
            repeat(reps)
        )))

        run_reps += reps

        results_pd = pd.DataFrame.from_records(results)
        #get the allele for which the median log likelihood is the highest
        highest_loglike = results_pd.groupby(3)[2].median().sort_values(ascending=False).head(n=1).index[0]
        #get the set of other alleles
        other_alleles = set(results_pd[3]) - set([highest_loglike])
        high_pval = 0
        n_tests = 0
        for b in other_alleles:
            if successful:
                try:
                    stat, pval = mannwhitneyu(
                        results_pd.groupby(3)[2].get_group(highest_loglike),
                        results_pd.groupby(3)[2].get_group(b),
                        alternative="greater"
                    )
                    n_tests += 1

                    high_pval = max(pval, high_pval)
                except ValueError:
                    successful = False
                    logging.warning(("SONiCS cannot resolve genotype for "
                                     "sample: %s block: %s"), name, block)
                    break

        #check if p_value threshold is satisfied
        high_pval *= n_tests
        if high_pval < pvalue_threshold:
            logging.debug("Will break! P-value: {}".format(high_pval))
            break
        #calculate additional repetitions
        reps = 4 * run_reps if reps_round % 2 == 1 else run_reps

        """make sure that the number of reps does not exceed
        the maximum number of reps"""
        reps = max_n_reps - run_reps if reps + run_reps > max_n_reps else reps

        reps_round += 1

    if verbose:
        results_pd.to_csv("{}_{}.txt".format(block, name),
                          sep="\t", header=False, index=False)

    if successful:
        results_groupedby_allele = results_pd[results_pd[2] > -999999].groupby(3)
        best_allele = results_groupedby_allele.get_group(highest_loglike)
        best_guess = best_allele.sort_values(0, ascending=False).head(n=1)
        ret = "{}\t{}\t{}\t{}\t{}\t{}".format(
            best_guess[0].item(), #ident
            best_guess[1].item(), #r2
            best_guess[2].item(), #log_likelihood
            high_pval, #highest p_value
            run_reps, #repetitions
            best_guess[3].item() #genotype n_reps/n_reps
        )
    else:
        ret = "{}\t{}\t{}\t{}\t{}\t{}".format(0, 0, 0, 1, run_reps, "./.")
    return ret


# FUNCTIONS RUN EVERY SIMULATION
def rsq(np.ndarray x, np.ndarray y):
    """Calculates the coefficient of determination.
    """
    a = np.array(sorted(x), dtype=int)
    b = np.array(sorted(y), dtype=int)
    x_ = a.mean()
    ss_tot = sum([(i - x_)**2 for i in a])
    ss_tot = 1e-16 if ss_tot == 0 else ss_tot
    ss_res = sum([i**2 for i in a - b])
    return 1 - ss_res / ss_tot


def one_repeat(str simulation_type, dict constants, tuple ranges,
               intermediate=None, how_many_reps=100):
    """Calls PCR simulation function, based on PCR products generates 
    genotype and calculates model statistics"""
    cdef int total_molecule, first, second, genotype_total, max_allele
    cdef dict parameters
    cdef str initial
    cdef float identified, r_squared, prob_a, noise_coef, noise_threshold
    cdef np.ndarray[DTYPE_t, ndim=1] alleles, alleles_nonzero, noise
    genotype_total = constants['genotype_total']
    noise_coef = constants['noise_coef']
    noise_threshold = constants['noise_threshold']
    max_allele = constants['max_allele']
    PCR_products = np.zeros(constants['max_allele'], dtype=DTYPE)
    parameters = generate_params(ranges, constants['up_preference'])
    alleles = constants['alleles']
    total_molecules = 0

    #Select initial molecules based on simulation parameters
    if simulation_type == "one_allele":
        if constants['random']:
            first = np.random.choice(alleles.nonzero()[0])
        else:
            # if not random, simulate PCR for an allele with max reads in the genotype
            first = np.argmax(alleles)
        # logging.info(PCR.max_allele)
        PCR_products[first] = constants['start_copies']
        initial = "{}/{}".format(first, first)
    else:
        if len(alleles.nonzero()[0]) == 1:
            raise Exception("Less then two alleles as starting conditions! Aborting.")
        else:
            if constants['random']:
                first, second = tuple(np.random.choice(alleles.nonzero()[0], 2))
            elif constants['half_random']:
                first = alleles.argmax()
                second = np.random.choice(alleles.nonzero()[0][alleles.nonzero()[0] != first])
            else:
                # if not random, simulate PCR for the alleles with max reads in the genotype
                first = alleles.argmax()
                second = np.argsort(-alleles)[1]
        PCR_products[first] = constants['start_copies'] / 2
        PCR_products[second] = constants['start_copies'] / 2
        #logging.debug("Starting PCR with alleles: %d, %d" %(first, second))
        if first < second:
            initial = "{}/{}".format(first, second)
        else:
            initial = "{}/{}".format(second, first)

    if noise_coef < noise_threshold:
        noise = np.copy(alleles)
        noise[noise > 0] = noise_coef * sum(PCR_products)
        PCR_products += noise

    else:
        if intermediate != None:
            """TODO: check if range ok => if starting conditions of 
            the pool included the alleles from present input"""
            dir_path = os.path.join(intermediate, "_".join(initial.split("/")))
            try:
                os.stat(dir_path)
            except:
                os.mkdir(dir_path)

            saved_pools = os.listdir(dir_path)

            if len(saved_pools) > how_many_reps:
                to_load = np.random.choice(saved_pools)
                PCR_products = np.load(os.path.join(dir_path, to_load))
            else:
                to_save = str(max([int(f.strip(".npy")) for f in saved_pools]) + 1) if len(saved_pools) > 0 else '1'
                # PCR simulation
                PCR_products = simulate(PCR_products, constants, parameters)
                np.save(os.path.join(dir_path, to_save), PCR_products)
        else:
            # PCR simulation
            PCR_products = simulate(PCR_products, constants, parameters)

    PCR_total_molecules = np.sum(PCR_products)

    # genotype generation
    mid = map(
        lambda al: list(repeat(
            al,
            np.random.binomial(genotype_total, PCR_products[al] / PCR_total_molecules)
        )),
        range(max_allele)
    )

    mid = list(chain.from_iterable(mid))
    loglike_a = mhl(alleles, PCR_products) if sum(PCR_products < alleles) == 0 else -999999

    # pick from PCR pool genotype
    try:
        """Simulate readout from PCR pool and compare it to the readout
        from the input."""
        readout = np.bincount(np.random.choice(mid, genotype_total))
        readout.resize(max_allele)
        # model statistics
        alleles_nonzero = alleles.nonzero()[0]
        identified = (sum([min(alleles[index], readout[index]) for index in alleles_nonzero])) / genotype_total
        r_squared = rsq(alleles, readout)
    except Exception:
        loglike_a = -999999
        identified = 0
        r_squared = 0

    report = [identified, r_squared, loglike_a, initial, noise_coef]
    return report


def simulate(np.ndarray products, dict constants, dict parameters):
    """Simulates PCR run, includes capture step if specified by PCR parameters
    """
    cdef int floor, ct, ct_up, al, n, namp, nslip, nup, ndown, ncorrect, cc
    cdef float up, down, efficiency, capture, pu, pd, pslip, pu_norm, floor_cap, cap_set, hit
    cdef np.ndarray[DTYPE_t, ndim=1] nzp, cs
    cc = constants['capture_cycle']
    floor = products.nonzero()[0][0] - constants['floor']
    floor = floor if floor > 1 else 1
    up = parameters['up']
    down = parameters['down']
    capture = parameters['capture']
    efficiency = parameters['efficiency']
    for cycle in range(1, constants['n_cycles']+1):
        # capture step
        if cycle == cc and len(products) > 1:
            nzp = products.nonzero()[0]
            cap_set = capture / (max(products) - min(nzp))
            floor_cap = 1 - (cap_set * len(nzp) / 2)
            n = 1
            for al in nzp:
                ct = products[al]
                hit = (floor_cap + cap_set * n) * ct
                ct_up = np.random.poisson(hit)
                ct_up = ct_up if ct_up < ct else ct
                products[al] = ct_up
                n += 1
        nzp = products.nonzero()[0]
        # cycle simulation for each allele
        for al in nzp:
            if al > floor:
                ct = products[al]
                prob_up = 1 - (1 - up) ** al
                prob_down = 1 - (1 - down) ** al
                prob_slip = 1 - (1 - prob_down) * (1 - prob_up)

                try:
                    prob_up_norm = prob_up / (prob_up + prob_down) 
                except ZeroDivisionError:
                    logging.warning("Encountered precision error!")
                    logging.warning("""allele: {}
                                    products: {}
                                    parameters: {}
                                    constants: {}""".format(al, 
                                                            products, 
                                                            parameters, 
                                                            constants))

                np.random.seed(np.random.randint(1, 4294967295))
                """number of molecules to which the polymerase bound,
                i.e. number of successes where number of trials is ct
                and probability is PCR efficiency"""
                mol_amp = np.random.binomial(ct, efficiency)
                """number of slips, where number of trials is the number
                of times the polymerase bound the molecule and
                probability is the probability of slippage"""
                mol_slip = np.random.binomial(mol_amp, prob_slip)
                """number of up stutters where number of trails is
                number of slips and the probability is the normalized
                probability of stutter up"""
                mol_up = np.random.binomial(mol_slip, prob_up_norm)
                mol_down = mol_slip - mol_up

                #number of molecules with correct number of repetitions
                mol_allele = mol_amp - mol_slip
                products[al] += mol_allele

                if mol_down > 0:
                    """polymerase slipped producing mol_down molecules
                    with one less repetition of the motif"""
                    products[al - 1] += mol_down
                if mol_down > 0:
                    """polymerase slipped producing mol_dup molecules
                    with one more repetition of the motif"""
                    products[al + 1] += mol_up
    return products
