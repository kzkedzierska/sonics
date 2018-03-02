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

__author__ = "Katarzyna Kedzierska"
__email__ = "kzk5f@virginia.edu"

# FUNCTIONS RUN ONLY ONCE/TWICE PER GENOTYPE
def get_alleles(genot_input):
    """get a dictionary with alleles as keys and number_of_reads as values
    from the input readout ('allele1|#;allele2|#').
    """
    max_allele = 500
    alleles = np.zeros(max_allele, dtype=DTYPE)
    for f in genot_input.split(';'):
        if f != "":
            pair = f.split("|")
            alleles[int(pair[0])] = int(pair[1])
    #alleles[alleles < max(alleles) * 0.01] = 0
    n_alleles = len(alleles.nonzero()[0])
    #logging.info(alleles)
    return alleles, max_allele, n_alleles

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

    return {'down': d, 'up': u, 'capture': c, 'efficiency': p}

def monte_carlo(max_n_reps, constants, ranges, all_simulation_params):
    """Runs Monte Carlo simulation of the PCR amplification until
    p_value threshold for the Mann Whitney test is reached or
    the number of repetition reaches the maximum.
    
    Arguments:
    max_n_reps -- upper limit for number of repetitions run
    constants -- constants throughout the simulations
    ranges -- ranges for generation PCR simulation specific parameters
    all_simulation_params -- parameters shared throughout all simulations
    Scheme:
    1) Run the first n repetitions.
    2) Calculate the highest p value for all the comparisons between
    set with the highest log likelihood median and others.
    3) Check if the highest p value with Bonferroni corrections is
    lower than the threshold. If true, stop. If false, go to step 1
    with new n now equal to 4*n if the it's the even round of
    repetitions or 2*n it's the odd. That way the sop will be made
    every 100, 500, 1000, 5000 etc. repetitions.
    """
    padjust = constants["padjust"]
    loglike_threshold = constants["loglike"]
    block = all_simulation_params['block']
    name = all_simulation_params['name']
    successful = False
    """successful - parameter that helps distinguish between
    the simulations needing more repetitions and the ones
    that are beyond the abilities of SONiCS"""
    results = list()
    run_reps = 0
    reps_round = 1

    if all_simulation_params['monte_carlo']:
        reps = max_n_reps
    else: 
        reps = 100 if max_n_reps > 100 else max_n_reps

    while run_reps < max_n_reps:

        results.extend(list(map(
            one_repeat,
            repeat(constants, reps),
            repeat(ranges),
            repeat(reps)
        )))

        run_reps += reps
        #group by the initial genotype
        results_colnames = [
            'ident', 
            'r_squared', 
            'log_like', 
            'genotype', 
            'noise_coef', 
            'down', 
            'up', 
            'capture', 
            'efficiency'
        ]
        results_pd = pd.DataFrame.from_records(results, 
                                               columns=results_colnames)

        if all_simulation_params['monte_carlo']:
            if all_simulation_params['save_report']:
                report_path = os.path.join(all_simulation_params['out_path'],
                                           "{}_{}.txt".format(block, name))
                results_pd.to_csv(report_path, index=False, sep="\t")

            best_guess = results_pd.sort_values(by="log_like", 
                                                ascending=False).head(n=1)
            genotype = best_guess["genotype"].item()

            genotype_pd = results_pd.groupby("genotype", as_index=False).get_group(genotype)
            quantiles = genotype_pd.quantile(0.75)

            ret = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            genotype, #genotype n_reps/n_reps
            best_guess["ident"].item(), #identity
            quantiles["ident"].item(), #identity quantile
            best_guess["r_squared"].item(), #r^2
            quantiles["r_squared"].item(), #r^2 quantile
            best_guess["log_like"].item(), #log_likelihood 
            quantiles["log_like"].item(), #log_like quantile
            run_reps #repetitions
            )

            return ret

        results_pd_succ = results_pd[results_pd.log_like > -999999].groupby("genotype")

        results_pd = results_pd.groupby("genotype", as_index=False)
        # print(results_pd.head(n=2))

        #check what's the minimum of simulations per genotype
        min_sim = results_pd_succ.size().sort_values().iloc[0]
        #check for minimum number of simulations
        if min_sim > 25:
            #get the medians for log likelihoods in groups
            results_medians = results_pd.median().sort_values(by="log_like", 
                                                              ascending=False)

            #get the allele for which the median log likelihood is the highest
            highest_loglike_median = results_medians['genotype'].iloc[0]
            second_highest_median = results_medians['genotype'].iloc[1]
            best_allele = results_pd.get_group(highest_loglike_median)
            # print(highest_loglike_median)
            # print(best_allele.head(n=2))
            second_best = results_pd.get_group(second_highest_median)
            # compare best likelihoods
            best_loglike_first = best_allele.sort_values(by="log_like", 
                                                         ascending=False)['log_like'].iloc[0]
            best_loglike_second = second_best.sort_values(by="log_like", 
                                                          ascending=False)['log_like'].iloc[0]
            best_loglike_ratio = best_loglike_first - best_loglike_second
            #compare median likelihoods 
            median_loglike_first = results_medians['log_like'].iloc[0]
            median_loglike_second = results_medians['log_like'].iloc[1]
            median_loglike_ratio = median_loglike_first - median_loglike_second
            #get the set of other alleles
            other_alleles = set(results_medians.genotype) - set([highest_loglike_median])
            high_pval = 0
            n_tests = 0
            for b in other_alleles:
                try:
                    stat, pval = mannwhitneyu(
                        best_allele.iloc[:,2],
                        results_pd.get_group(b).iloc[:,2],
                        alternative="greater"
                    )

                    high_pval = max(pval, high_pval)
                except ValueError:
                    high_pval = 1

                n_tests += 1

            #Bonferroni correction in its essence
            high_pval *= n_tests
            #check if p_value threshold is satisfied
            #get the medians for log likelihoods in groups
            results_maxs = results_pd.max().sort_values(by="log_like", 
                                                        ascending=False)
            highest_loglike_max = results_maxs['genotype'].iloc[0]

            if highest_loglike_max == highest_loglike_median:
                if (
                    high_pval < padjust and 
                    median_loglike_ratio > loglike_threshold and
                    best_loglike_ratio > loglike_threshold
                ):
                    successful = True
                    logging.debug("Will break! P-value: {}".format(high_pval))
                    break

            else:
                #TEST BEST NOW
                second_highest_max = results_maxs['genotype'].iloc[1]
                best_allele_max = results_pd.get_group(highest_loglike_max)
                # print(best_allele_max.head(n=2))
                second_best_max = results_pd.get_group(second_highest_max)
                # compare best likelihoods
                best_loglike_first_max = results_maxs['log_like'].iloc[0]
                best_loglike_second_max = results_maxs['log_like'].iloc[1]
                best_loglike_ratio_max = best_loglike_first_max - best_loglike_second_max
                #compare median likelihoods 
                median_loglike_first_max = best_allele_max.median()['log_like']
                median_loglike_second_max = second_best_max.median()['log_like']
                median_loglike_ratio_max = median_loglike_first_max - median_loglike_second_max
                #get the set of other alleles
                other_alleles = set(results_maxs.genotype) - set([highest_loglike_max])
                high_pval_max = 0
                for b in other_alleles:
                    try:
                        stat, pval = mannwhitneyu(
                            best_allele.iloc[:,2],
                            results_pd.get_group(b).iloc[:,2],
                            alternative="greater"
                        )

                        high_pval = max(pval, high_pval_max)
                    except ValueError:
                        high_pval_max = 1

                #Bonferroni correction in its essence
                high_pval_max *= n_tests

                if (
                    high_pval < padjust and 
                    median_loglike_ratio > loglike_threshold and
                    best_loglike_ratio > loglike_threshold
                ):
                    if (
                        high_pval_max > padjust or 
                        median_loglike_ratio_max < loglike_threshold or
                        best_loglike_ratio_max < loglike_threshold
                    ):
                        successful = True
                        logging.debug("Will break! P-value: {}".format(high_pval))
                        break
                else:
                    high_pval = high_pval_max
                    median_loglike_ratio = median_loglike_ratio_max
                    best_loglike_ratio = best_loglike_ratio_max
                    best_allele = best_allele_max

                    if (
                        high_pval_max < padjust or 
                        median_loglike_ratio_max > loglike_threshold or
                        best_loglike_ratio_max > loglike_threshold
                    ):
                        successful = True
                        logging.debug("Will break! P-value: {}".format(high_pval))
                        break

        #calculate additional repetitions
        reps = 4 * run_reps if reps_round % 2 == 1 else run_reps

        # make sure that the number of reps does not exceed the maximum 
        reps = max_n_reps - run_reps if reps + run_reps > max_n_reps else reps

        reps_round += 1

    if all_simulation_params['save_report']:
        results_pd_csv = pd.DataFrame.from_records(results, 
                                                   columns=results_colnames)
        report_path = os.path.join(all_simulation_params['out_path'],
                                   "{}_{}.txt".format(block, name))
        results_pd_csv.to_csv(report_path, index=False, sep="\t")

    high_pval = high_pval if high_pval < 1 else 1
    not_zero = best_allele.log_like > -999999
    best_guess = best_allele[not_zero].sort_values("log_like", 
                                                   ascending=False).head(n=1)

    if best_guess.empty:
        filt = "no_success"
        #this can happen if there is noise from very distant alleles
        ret = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            "./.", #genotype
            ".", #identity
            ".", #r^2
            ".", #log_like
            filt, #FILTER
            ".", #Mann-Whitney U test, p_val
            ".", #best loglike
            ".", #median loglike
            run_reps #reps
        )
        return ret

    if not successful:
        conditions = [
            "MWU_test" if high_pval > padjust else "",
            "best_ratio" if best_loglike_ratio < constants["loglike"] else "",
            "median_ratio" if median_loglike_ratio < constants["loglike"] else ""
        ]
        filt = ",".join([cond for cond in conditions if cond != ""])
    else:
        filt = "PASS"

    ret = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            best_guess["genotype"].item(), #genotype n_reps/n_reps
            best_guess["ident"].item(), #identity
            best_guess["r_squared"].item(), #r^2
            best_guess["log_like"].item(), #median log_likelihood 
            filt, #FILTER
            high_pval, #highest p_value
            best_loglike_ratio, #best loglikelihood ratio
            median_loglike_ratio, #median loglikelihood ratio
            run_reps #repetitions
        )

    return ret


# FUNCTIONS RUN EVERY SIMULATION
def rsq(np.ndarray true_values, np.ndarray pred_values):
    """Calculates the coefficient of determination between the truth (x) and
    prediction (y).
    """
    fr = min(
        true_values.nonzero()[0][0], 
        pred_values.nonzero()[0][0]
    )
    to = 1 + max(
        true_values.nonzero()[0][-1], 
        pred_values.nonzero()[0][-1]
    )
    true_values = true_values[fr:to]
    pred_values = pred_values[fr:to]
    true_mean = true_values.mean()
    ss_tot = sum([(i - true_mean) ** 2 for i in true_values])
    ss_tot = 1e-16 if ss_tot == 0 else ss_tot
    ss_res = sum([i ** 2 for i in true_values - pred_values])
    return 1 - ss_res / ss_tot

def one_repeat(dict constants, tuple ranges,
               int how_many_reps=100):
    """Calls PCR simulation function, based on PCR products generates 
    genotype and calculates model statistics"""
    cdef int total_molecule, first, second, genotype_total, max_allele, floor
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
    if constants['floor'] == -1:
        floor = 1
    else:
        floor = constants['alleles'].nonzero()[0][0] - constants['floor']
        floor = floor if floor > 1 else 1

    if len(alleles.nonzero()[0]) == 1:
            raise Exception(("Less then two alleles as starting conditions!"
                             " Aborting."))

    if constants['random']:
        first, second = tuple(np.random.choice(alleles.nonzero()[0], 2))
    else:
        first = alleles.argmax()
        second = np.random.choice(alleles.nonzero()[0])

    PCR_products[first] += constants['start_copies'] / 2
    PCR_products[second] += constants['start_copies'] / 2

    if first < second:
        initial = "{}/{}".format(first, second)
    else:
        initial = "{}/{}".format(second, first)

    if noise_coef < noise_threshold:
        noise = np.copy(alleles)
        noise[noise > 0] = noise_coef * sum(PCR_products)
        PCR_products += noise
        PCR_products = simulate(PCR_products, constants, parameters, floor)

    else:
        # PCR simulation
        PCR_products = simulate(PCR_products, constants,
                                parameters, floor)

    PCR_total_molecules = np.sum(PCR_products)

    #genotype generation
    mid = []
    for allele in range(max_allele):
        """The binomial distribution is frequently used to model the number of 
        successes in a sample of size n drawn with replacement from 
        a population of size N. If the sampling is carried out without
        replacement, the draws are not independent and so the resulting 
        distribution is a hypergeometric distribution, not a binomial one.
        However, for N much larger than n, the binomial distribution 
        remains a good approximation, and is widely used. [Wikipedia]"""

        n_times = np.random.binomial(genotype_total, 
                                     PCR_products[allele] / PCR_total_molecules)
        allele_molecules = list(repeat(allele, n_times))
        mid.extend(allele_molecules)

    
    if sum(alleles > PCR_products) != 0:
        loglike_a = -999999
    else:
        loglike_a = mhl(alleles, PCR_products) 

    """Simulate readout from PCR pool and compare it to the readout
    from the input."""
    try:
        readout = np.bincount(np.random.choice(mid, genotype_total))
        readout.resize(max_allele)
    except ValueError:
        #
        readout = np.zeros(max_allele, dtype=DTYPE)
    
    # model statistics
    alleles_nonzero = alleles.nonzero()[0]
    if readout.nonzero()[0].size == 0:
        #this happens if the input genotype is very small
        #basically the fragments did not get sequenced
        identity = 0
        r_squared = -999999
    else:
        identity = ((sum([min(alleles[i], readout[i]) for i in alleles_nonzero]))
                      / genotype_total)
        r_squared = rsq(alleles, readout)

    report = [identity, r_squared, loglike_a, initial, noise_coef]
    prmtrs = [
        parameters['down'], 
        parameters['up'], 
        parameters['capture'], 
        parameters['efficiency']
    ]
    report.extend(prmtrs)
    return report


def simulate(np.ndarray products, dict constants, dict parameters, int floor):
    """Simulates PCR run, includes capture step if specified by PCR parameters
    """
    cdef int ct, ct_up, al, n, namp, nslip, nup, ndown, ncorrect, cc
    cdef float efficiency, capture, floor_cap, cap_set, hit
    cdef double up, down, prob_up, prob_down, prob_slip, pu_norm
    cdef long seed_n
    cdef np.ndarray[DTYPE_t, ndim=1] nzp, cs
    cc = constants['capture_cycle']
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
            if al >= floor:
                seed_n = np.random.randint(1, 4294967295)
                ct = products[al]
                prob_up = 1 - (1 - up) ** al
                prob_down = 1 - (1 - down) ** al
                prob_slip = 1 - (1 - prob_down) * (1 - prob_up)

                try:
                    prob_up_norm = prob_up / (prob_up + prob_down) 
                except ZeroDivisionError:
                    logging.warning(("Encountered precision error!\n"
                                     "allele: %s\n"
                                     "parameters: %s\n"
                                     "constants: %s\n"
                                     "prob_up: %.s\n"
                                     "prob_down: %.s\n"
                                     "prob_down + prob_up: %s\n"
                                     "up: %s\n"
                                     "down: %s\n"
                                     "seed: %s\n"), al, 
                                    parameters, constants, prob_up, prob_down, 
                                    prob_up + prob_down, up, down, seed_n) 
                
                np.random.seed(seed_n)
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
            else:
                mol_amp = np.random.binomial(ct, efficiency)
                products[al] += mol_amp
    return products
