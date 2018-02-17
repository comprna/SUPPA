# -*- coding: utf-8 -*-
"""
Created on Wed May 25 04:20:00 CEST 2016

@authors: Juan C Entizne
@email: juancarlos.entizne01[at]estudiant.upf.edu

Modified by Juan L. Trincado
@email: juanluis.trincado[at].upf.edu

"""

import os
import sys
import math
import logging
import warnings
import numpy as np
import pandas as pd
from functools import reduce
from bisect import bisect_left
from collections import defaultdict
from itertools import combinations, islice
from scipy.stats import wilcoxon, mannwhitneyu
from statsmodels.sandbox.stats.multicomp import multipletests
from statsmodels.distributions.empirical_distribution import ECDF


def progressbar(prefix, i, lst_len):
    print(prefix, " ", "%d / %d. " % (i+1, lst_len), "%.2f%% completed." % ((i/lst_len)*100), end="\r", flush=True)


def flatten(d):

    try:
        fd = {k: sum(v, []) for k, v in d.items()}
    except Exception as e:
        pass

    try:
        fd = {k: sum(v, ()) for k, v in d.items()}
    except Exception as e:
        pass

    return fd


def create_dict(arg):

    d = defaultdict(list)
    with open(arg) as fh:
        next(fh)
        for event in fh:
            line = event.split()
            event_id = line[0]
            event_vals = []
            for val in line[1:]:
                try:
                    event_vals.append(float(val))
                except:
                    event_vals.append(float('nan'))
            d[event_id].append(event_vals)

    return flatten(d)


def get_psi_values(dict1, dict2):

    psi_values = defaultdict(list)
    for d in (dict1, dict2):
        for key, value in d.items():
            psi_values[key].append(value)

    return psi_values

def get_proportion_nans(psi_list):

    count = 0
    for x in psi_list:
        if(math.isnan(x)):
            count += 1
    size = len(psi_list)
    return float(count)/float(size)


def calculate_delta_psi(psi_values, median, nan_th):

    abs_dt, dt, discarded_events = (defaultdict(list) for _ in range(3))
    for event in psi_values:

        #Get the proportion of missing values in an event
        prop0 = get_proportion_nans(psi_values[event][0])
        prop1 = get_proportion_nans(psi_values[event][1])

        # event will be excluded if any of the proportions overtake the nan_threshold
        if nan_th < prop0 or nan_th < prop1:
            discarded_events[event].append([float("nan"), 1.0000000000])

        else:

            #if passes the threshold, remove all the nan values form each list
            psi_values_0 = [x for x in psi_values[event][0] if str(x) != 'nan']
            psi_values_1 = [x for x in psi_values[event][1] if str(x) != 'nan']

            if median:

                abs_dpsi_val = abs(np.nanmedian(psi_values_1) - np.nanmedian(psi_values_0))
                abs_dt[event].append(abs_dpsi_val)

                dpsi_val = np.nanmedian(psi_values_1) - np.nanmedian(psi_values_0)
                dt[event].append(dpsi_val)

            else:
                # Ignore empty slice warning when calculating the mean/median
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore', r'Mean of empty slice')

                    abs_dpsi_val = abs(np.nanmean(psi_values[event][1]) - np.nanmean(psi_values[event][0]))
                    abs_dt[event].append(abs_dpsi_val)

                    dpsi_val = np.nanmean(psi_values[event][1]) - np.nanmean(psi_values[event][0])
                    dt[event].append(dpsi_val)

    # Flatten the list of list of the dictionary values
    dpsi_abs_values = {k: sum(v) for k, v in abs_dt.items()}
    dpsi_values = {k: sum(v) for k, v in dt.items()}

    return dpsi_abs_values, dpsi_values, flatten(discarded_events)


def get_events_transcripts(ioe):

    td = defaultdict(list)
    with open(ioe) as fh_ioe:
        for line in fh_ioe:
            event_id_ioe = line.split()[2]
            tot_transcripts = line.split()[4].split(',')
            td[event_id_ioe].append(tot_transcripts)

    transcripts_values = flatten(td)

    return transcripts_values


def get_tpm_values(tpm1_values, tpm2_values, transcripts_values):

    discarded_transcript_events = []
    tpm_values = defaultdict(list)
    for tpm_dt in (tpm1_values, tpm2_values):
        for event in transcripts_values:
            transcript_vals = []
            for transcript in transcripts_values[event]:
                try:
                    transcript_vals.append(tpm_dt[transcript])
                except:
                    discarded_transcript_events.append(event)
            tpm_values[event].append(transcript_vals)

    return tpm_values


def calculate_transcript_abundance(tpm_values, tpm_th):

    if(tpm_th!=0):
        tpm_th_log10 = math.log10(tpm_th)
    else:
        tpm_th_log10 = -float('Inf')

    temp_between_conditions_logtpm = defaultdict(list)
    for event in tpm_values:

        conditions_average_logtpm = []
        for transcript_vals in tpm_values[event]:

            # Group the TPMs according to their replicate of origin
            replicates_transcript_values = list(zip(*transcript_vals))

            try:
                replicates_logtpms = [math.log10(sum(rep_tpms)) for rep_tpms in replicates_transcript_values]
                average_replicate_transcript_abundance = sum(replicates_logtpms)/len(replicates_logtpms)

                if average_replicate_transcript_abundance >= tpm_th_log10:
                    conditions_average_logtpm.append(average_replicate_transcript_abundance)

            except:
                pass

        # Filter out the events for which it was not possible to calculate log10(sum(TPMs)) for one of the conditions
        if len(conditions_average_logtpm) == 2.0:
            between_conditions_average_transcript_abundance = 0.5 * sum(conditions_average_logtpm)

            if between_conditions_average_transcript_abundance >= tpm_th_log10:
                temp_between_conditions_logtpm[event].append(between_conditions_average_transcript_abundance)
        else:
            pass

    # Flatten the list of list of the dictionary values
    between_conditions_avglogtpm = {k: sum(v) for k, v in temp_between_conditions_logtpm.items() if v[0]}

    return between_conditions_avglogtpm


def merge_dict(d1, d2):

    md = defaultdict(list)
    for d in (d1, d2):
        for key, value in d.items():
            md[key].append(value)

    merged_dict = defaultdict(list)
    for k in md:
        if len(md[k]) == 2.0:
            merged_dict[k].append(md[k])
        else:
            pass

    return flatten(merged_dict)


def get_closest_number(lst, n):
    """
    Assumes lst is sorted. Returns closest value to n.
    If two numbers are equally close, return the smallest number.
    Source: http://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value/
    """
    pos = bisect_left(lst, n)
    if pos == 0:
        return lst[0]
    if pos == len(lst):
        return lst[-1]
    before = lst[pos - 1]
    after = lst[pos]
    if after - n < n - before:
        return after
    else:
        return before


def slice_list(lst, index, slice_len):

    half_len = int(slice_len*0.5)
    diff = index - half_len

    if diff < 0:

        left_bound = 0
        right_bound = index + half_len + (-diff) + 1

    elif index + half_len >= len(lst):

        upper_diff = index + half_len - len(lst) + 1

        left_bound = diff - upper_diff
        right_bound = index + half_len + 1

    else:
        left_bound = diff
        right_bound = index + half_len + 1

    local_dpsi = lst[left_bound:right_bound]

    return local_dpsi


def calculate_empirical_pvalue(local_area, dpsi_abs_value):

    abs_local_area = [abs(val) for val in local_area]

    ecdf = ECDF(abs_local_area)

    # It is divided by 2 because we are using abs(deltaPSI) values and therefore it is a one-tailed test
    event_pvalue = (1.0 - ecdf(dpsi_abs_value)) * 0.5

    return event_pvalue


def calculate_between_conditions_distribution(cond1, cond2, tpm1, tpm2, ioe, save_tpm, median, tpm_th, nan_th, output):

        cond1_psi_values = create_dict(cond1)
        cond2_psi_values = create_dict(cond2)
        psi_values = get_psi_values(cond1_psi_values, cond2_psi_values)

        dpsi_abs_values, dpsi_values, discarded_events = calculate_delta_psi(psi_values, median, nan_th)

        transcripts_values = get_events_transcripts(ioe)

        tpm1_values = create_dict(tpm1)
        tpm2_values = create_dict(tpm2)
        tpm_values = get_tpm_values(tpm1_values, tpm2_values, transcripts_values)

        between_conditions_avglogtpm = calculate_transcript_abundance(tpm_values, tpm_th)

        if save_tpm:
            #Save between_conditions_avglogtpm object
            print("Saving between_conditions_avglogtpm...")
            output = output + "_avglogtpm.tab"
            outFile = open(output, 'w')
            for key in between_conditions_avglogtpm.keys():
                line = key + "\t" + str(between_conditions_avglogtpm[key]) + "\n"
                outFile.write(line)
            outFile.close()

            print("Saved "+output)

        between_conditions_absdpsi_logtpm = merge_dict(dpsi_abs_values, between_conditions_avglogtpm)

        return between_conditions_absdpsi_logtpm, psi_values, tpm_values, dpsi_abs_values, dpsi_values, discarded_events


def create_replicates_distribution(between_conditions_distribution, psi_dict, tpm_dict):

    unsorted_replicates_distribution, unsorted_rep_dist_for_plot = ([] for _ in range(2))
    for event in between_conditions_distribution.keys():
        conds_psi_rep_tpms = list(zip(psi_dict[event], tpm_dict[event]))

        for cond_psi_trans in conds_psi_rep_tpms:
            psis = cond_psi_trans[0]
            trans_tpms = cond_psi_trans[1]

            # Group the TPMs according to their replicate of origin
            rep_trans = list(zip(*trans_tpms))

            rep_psi_trans_lst = list(zip(psis, rep_trans))

            cond_psi_trans_lst = []
            for psi_trans in rep_psi_trans_lst:
                rep_psi_val = psi_trans[0]
                trans = psi_trans[1]

                try:
                    rep_logtpm = math.log10(sum(trans))
                    rep_psi_logtpm_pair = (rep_psi_val, rep_logtpm)
                    cond_psi_trans_lst.append(rep_psi_logtpm_pair)

                except Exception as e:
                    pass

            psi_trans_paired = list(combinations(cond_psi_trans_lst, r=2))

            for pair in psi_trans_paired:
                # A rep_pair contains (replicate_psi_value, replicate_avg_log10_tpm_value)
                rep1_pair = pair[0]
                rep2_pair = pair[1]

                try:
                    rep_delta_psi = rep2_pair[0] - rep1_pair[0]
                    rep_pair_avg_logtpm = (rep1_pair[1] + rep2_pair[1]) * 0.5
                    unsorted_replicates_distribution.append((rep_delta_psi, rep_pair_avg_logtpm))

                    unsorted_rep_dist_for_plot.append((event, rep_delta_psi, rep_pair_avg_logtpm))

                except Exception as e:
                    pass

    # It's important to sort because get_closest_number assume a sorted list
    replicates_distribution = sorted(unsorted_replicates_distribution, key=lambda x: x[1])

    # List converted to numpy array for better performance
    return np.array(replicates_distribution)


def get_local_distribution(ev_logtpm, replicates_distribution, replicates_logtpms, windows_len):

    close_rep_logtpm = get_closest_number(replicates_logtpms, ev_logtpm)
    local_dist = slice_list(replicates_distribution, replicates_logtpms.index(close_rep_logtpm), windows_len)

    return local_dist


def calculate_events_pvals(between_conditions_distribution,
                           replicates_distribution, area, abs_dpsi_dict, cutoff):

    replicates_logtpms = [event[1] for event in replicates_distribution]

    lst_len = len(between_conditions_distribution)

    uncorrected_pvals, event_lst = ([] for _ in range(2))
    for i, event in enumerate(between_conditions_distribution):

        progressbar("Calculating events empirical p-value:", i, lst_len)

        between_cond_obs_dpsi = abs_dpsi_dict[event]
        ev_logtpm = between_conditions_distribution[event][1]
        local_dist = get_local_distribution(ev_logtpm, replicates_distribution, replicates_logtpms, area)
        local_dpsi = [e[0] for e in local_dist]

        if -cutoff < between_cond_obs_dpsi < cutoff:
            event_pval = 1.0
            uncorrected_pvals.append(event_pval)
            event_lst.append(event)

        else:
            event_pval = calculate_empirical_pvalue(local_dpsi, between_cond_obs_dpsi)
            uncorrected_pvals.append(event_pval)
            event_lst.append(event)

    progressbar("Calculating events empirical p-value:", i+1, lst_len)
    print("\nDone!\n")

    return event_lst, uncorrected_pvals


def nan_eliminator(lst1, lst2, paired):

    if paired:
        z = list(zip(lst1, lst2))
        try:
            l1, l2 = zip(*[e for e in z if not math.isnan(e[0]) and not math.isnan(e[1])])
        except:
            l1, l2 = [], []
    else:
        l1 = [e for e in lst1 if not math.isnan(e)]
        l2 = [e for e in lst2 if not math.isnan(e)]

    return l1, l2


def pval_multiple_test_corrector(pval_dict, alpha):

    pval_lst, raw_pvals = ([] for _ in range(2))
    for event in pval_dict:
        pval_lst.append((event, pval_dict[event]))
        raw_pvals.append(pval_dict[event])

    _, pvals_corrected, _, _ = multipletests(raw_pvals, method='fdr_bh', alpha=alpha)

    unflat_corrected_pval_dict = defaultdict(list)
    for i, j in zip(pval_lst, pvals_corrected):
        unflat_corrected_pval_dict[i[0]].append(j)

    corrected_pval_dict = {k: sum(v) for k, v in unflat_corrected_pval_dict.items()}

    return corrected_pval_dict


def write_temp_output_files(dpsi_pval_dict, output, i, cond1_name, cond2_name):

    # Order alphabetically
    results_lst = sorted(dpsi_pval_dict.items(), key=lambda x: x[0])

    with open("%s.dpsi.temp.%d" % (output, i), 'w+') as fh:
        cond_id = cond1_name+"-"+cond2_name
        f_line = "Event_id\t%s_dPSI\t%s_p-val\n" % (cond_id, cond_id)
        fh.write(f_line)
        for event in results_lst:
            line = "%s\t%.10f\t%.10f\n" % (event[0], event[1][0], event[1][1])
            fh.write(line)


def merge_temp_output_files(output):

    # Set working directory
    if os.path.isabs(output):
        current_path = os.path.dirname(output)+"/"
    else:
        current_path = os.getcwd()+"/"

    dpsi_files = []
    for fl in os.listdir(current_path):
        if ".dpsi.temp." in fl:
            dpsi_files.append(current_path+fl)

    dpsi_files.sort(key=lambda x: x[-1])

    df_lst = []
    for lst in dpsi_files:
        df = pd.read_table(lst, sep='\t', index_col=0, header=0)
        df_lst.append(df)

        merged_dpsi_results = reduce(lambda left, right: pd.merge(left, right,
                                                        left_index=True, right_index=True,
                                                        how='outer'), df_lst)

        header = merged_dpsi_results.columns.values

        with open("%s.dpsi" % output, "w+") as fh:
            ln = "\t".join(header)
            fh.write(ln+"\n")

        with open("%s.dpsi" % output, "a") as fh:
            merged_dpsi_results.to_csv(fh, sep="\t", na_rep="nan", header=False)

    # Delete temp filesdis
    for fl in os.listdir(current_path):
        if ".dpsi.temp." in fl:
            os.remove(current_path+fl)

    return os.path.abspath("%s.dpsi" % output)


def write_psivec_file(psi_lst, output):

    df_lst = []
    for fl in psi_lst:
        df = pd.read_table(fl, sep='\t', skiprows=[0], index_col=0, header=None)

        old_header = df.columns.values
        new_header = [os.path.basename(fl).split(".")[0]+"_"+str(col_id) for col_id in old_header]
        df.rename(columns=dict(zip(old_header, new_header)), inplace=True)

        df_lst.append(df)

    merged_psi_results = reduce(lambda left, right: pd.merge(left, right,
                                                            left_index=True, right_index=True,
                                                            how='outer'), df_lst)

    header = merged_psi_results.columns.values

    with open("%s.psivec" % output, "w+") as fh:
            ln = "\t".join(header)
            fh.write(ln+"\n")

    with open("%s.psivec" % output, "a") as fh:
            merged_psi_results.to_csv(fh, sep="\t", na_rep="nan", header=False)

    return os.path.abspath("%s.psivec" % output)


def empirical_test(cond1, tpm1, cond2, tpm2, ioe, area, cutoff, save_tpm, median, tpm_th, nan_th, output):

    between_conditions_distribution, psi_dict, tpm_dict, abs_dpsi_dict, dpsi_vals, discarded_dt \
        = calculate_between_conditions_distribution(cond1, cond2, tpm1, tpm2, ioe, save_tpm, median, tpm_th, nan_th, output)

    replicates_distribution = create_replicates_distribution(between_conditions_distribution, psi_dict, tpm_dict)

    event_lst, uncorrected_pvals = calculate_events_pvals(between_conditions_distribution,
                                              replicates_distribution, area, abs_dpsi_dict, cutoff)

    return event_lst, uncorrected_pvals, dpsi_vals, discarded_dt


def classical_test(cond1, cond2, paired, median):

    if paired:
        test = wilcoxon
    else:
        test = mannwhitneyu

    cond1_dict = create_dict(cond1)
    cond2_dict = create_dict(cond2)
    psi_dict = get_psi_values(cond1_dict, cond2_dict)

    uncorrected_pvals, event_lst = ([] for _ in range(2))
    for i in psi_dict:
        l1, l2 = nan_eliminator(psi_dict[i][0], psi_dict[i][1], paired)
        try:
            stat_score, pval = test(l1, l2)
            event_lst.append(i)
            uncorrected_pvals.append(pval)
        except:
            # Wilcoxon fails when the two vectors given are identical, even if the number of reps is enough
            # For those cases we decided to give a pval of 1.0
            pval = 1.0
            event_lst.append(i)
            uncorrected_pvals.append(pval)

    _, dpsi_vals, discarded_dt = calculate_delta_psi(psi_dict, median)

    return event_lst, uncorrected_pvals, dpsi_vals, discarded_dt


def sliding_windows(seq, n=2):
    "Returns a sliding sliding_windows (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def pval_multiple_test_corrector_by_gene(pvals_dict, alpha):

    evid_pvals_dict, corrected_pval_dict = (defaultdict(list) for _ in range(2))
    for ev_id in pvals_dict.keys():
        gene = ev_id.split(";")[0]

        evid_pvals_dict[gene].append((ev_id, pvals_dict[ev_id]))

    for gene in evid_pvals_dict:
        events, raw_pvals = zip(*evid_pvals_dict[gene])
        _, pvals_corrected, _, _ = multipletests(raw_pvals, method='fdr_bh', alpha=alpha)

        evid_corrected_pvals_list = list(zip(events, pvals_corrected))

        for evid_pval in evid_corrected_pvals_list:
            corrected_pval_dict[evid_pval[0]] = evid_pval[1]

    return corrected_pval_dict


def multiple_test_correction(event_lst, uncorrected_pvals, alpha):

    _, corrected_pvals, _, _ = multipletests(uncorrected_pvals, alpha=alpha, method='fdr_bh', returnsorted=False)

    corrected_pvals_dict = {k: v for k, v in zip(event_lst, corrected_pvals)}

    return corrected_pvals_dict


def convert_to_log10pval(dpsi_pval_values, sig_threshold=0.05, dpsi_threshold=0.05):
    '''
    Convert p-values into -log10_pvalues
    If p-value = 0 then the -log10_pvalue is calculated using the half of the lowest p-value in the dictionary
    '''

    # values = [dpsi, pvalue]
    min_pval = min([values[1] for values in dpsi_pval_values.values() if values[1] != 0.0])*0.5
    log10_vals, events_sig = ([] for _ in range(2))
    for key, values in dpsi_pval_values.items():
        try:
            log10_vals.append((key, -math.log10(values[1])))
        except:
            log10_vals.append((key, -math.log10(min_pval)))

        if values[1] <= sig_threshold and abs(values[0]) >= dpsi_threshold:
            sig = "significant"
        else:
            sig = "not-significant"

        events_sig.append((key, sig))

    log10_pvalues = {k: v for k, v in log10_vals}
    events_significance = {k: v for k, v in events_sig}

    return log10_pvalues, events_significance


def multiple_conditions_analysis(method, psi_lst, tpm_lst, ioe, area, cutoff, paired, gene_cor, alpha,
                                 save_tpm, comb, median, tpm_th, nan_th, output):

    # Setting logging preferences
    logger = logging.getLogger(__name__)

    z_lst = list(zip(psi_lst, tpm_lst))

    # Pair the conditions in a sequential manner or in all possible combinations
    if not comb:
        seq_lst = list(sliding_windows(z_lst))
    else:
        seq_lst = list(combinations(z_lst, r=2))

    for i, paired_cond_tpm in enumerate(seq_lst):
        cond1, tpm1 = paired_cond_tpm[0][0], paired_cond_tpm[0][1]
        cond2, tpm2 = paired_cond_tpm[1][0], paired_cond_tpm[1][1]

        cond1_name, cond2_name = os.path.basename(cond1).split(".")[0], os.path.basename(cond2).split(".")[0]

        print("Calculating differential analysis between conditions: %s and %s "
              "" % (cond1_name, cond2_name))

        if method == 'empirical':
            #event_lst, uncorrected_pvals, dpsi_vals, discarded_dt = empirical_test(cond1, tpm1, cond2, tpm2,ioe, area, cutoff)
            event_lst, uncorrected_pvals, dpsi_vals, discarded_dt = empirical_test(cond1, tpm1, cond2, tpm2, ioe, area,
                                                                                   cutoff, save_tpm, median, tpm_th,
                                                                                   nan_th, output)
            corrected_pvals_dict = {k: v for k, v in zip(event_lst, uncorrected_pvals)}

        elif method == 'classical':
            event_lst, uncorrected_pvals, dpsi_vals, discarded_dt = classical_test(cond1, cond2, paired, median)
            corrected_pvals_dict = multiple_test_correction(event_lst, uncorrected_pvals, alpha)

        else:
            logger.error("Unknown error: {}".format(sys.exc_info()))
            sys.exit(1)

        pvals_dict = {k: v for k, v in zip(event_lst, uncorrected_pvals)}

        if gene_cor:
            corrected_pvals_dict = pval_multiple_test_corrector_by_gene(pvals_dict, alpha)
        else:
            corrected_pvals_dict = pvals_dict

        dpsi_pval_values = merge_dict(dpsi_vals, corrected_pvals_dict)

        dpsi_pval_values.update(discarded_dt)

        # Commented out while creating plot functionality and fixing "min() empty" bug
        # log10_pvalues, events_significance = convert_to_log10pval(dpsi_pval_values, sig_threshold=0.05)

        write_temp_output_files(dpsi_pval_values, output, i, cond1_name, cond2_name)

    dpsi_fl_path = merge_temp_output_files(output)
    psivec_fl_path = write_psivec_file(psi_lst, output)
