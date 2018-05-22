#!/usr/bin/python3
import copy
import getopt
import os
import sys

import numpy as np
import pandas as pd
import scipy.stats  # For fisher exact test
import statsmodels.sandbox.stats.multicomp  # For Bejamini Hochberg FDR correction
from analytical_settings import HAP1_control_1, HAP1_control_2, HAP1_control_3, HAP1_control_4, fisher_plot_title
from global_functions import read_csv, write_csv, plot


##################################
# Function 'test_against_control'
##################################
# Arguments: data       <pandas DataFrame with normalized sense and antisense counts of replicate>
#            control    <pandas DataFrame with normalized sense and antisense counts of control>
# Returns:   data       <extended version of data with non-fdr-corrected and fdr-corrected p-values of
#                       replicate against control>
##################################
def test_against_control(data, control):  # control is tuple, [0] = name, [1] = filename
    controldata = read_csv(control[1])[['gene', 'sense', 'antisense']]
    controldata.columns = ['gene', 'sense_control', 'antisense_control']
    data = pd.merge(data, controldata, on='gene', how='inner')  # This means that if a gene isn't present in any of the
                                                        # 4 controls, we won't be able to find it back in the screen!
    fdr_corrected_pv_colname = ''.join(['fcpv_control_', control[0].split("_")[-1]])
    non_fdr_corrected_pv_colname = ''.join(['pv_control_', control[0].split("_")[-1]])
    data = data.fillna(0)
    print("Performing Fisher Exact tests...", end='')
    data[non_fdr_corrected_pv_colname] = data.apply(
        lambda r: scipy.stats.fisher_exact([[r.sense, r.antisense], [r.sense_control, r.antisense_control]])[1],
        axis=1)
    print("done\nApplying Multiple Testing correction using Benjamini Hochberg method...", end='')
    data[fdr_corrected_pv_colname] = statsmodels.sandbox.stats.multicomp.multipletests(
        data[non_fdr_corrected_pv_colname], alpha=0.05, method='fdr_bh', returnsorted=False)[1]
    data = data[['gene', fdr_corrected_pv_colname, non_fdr_corrected_pv_colname]]
    print('done')
    return data


def main(argv):
    syntax = 'test_replicate_vs_controls.py -r <directory where replicate resided>'
    try:
        opts, args = getopt.getopt(argv, "hr:", ["replicate="])
    except getopt.GetoptError:
        print(syntax)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(syntax)
            sys.exit()
        elif opt in ("-r", "--replicate"):
            sampledir = arg
        else:
            assert False, "unhandled option"

    # Create dictionary to store the names of the controls (keys) and their associated records (values)
    controls = {'HAP1_control_1': HAP1_control_1, 'HAP1_control_2': HAP1_control_2, 'HAP1_control_3': HAP1_control_3,
                'HAP1_control_4': HAP1_control_4}
    sampledir = os.path.join(sampledir, '')
    sample = read_csv("".join([sampledir, "sense_vs_antisense_counts_normalized.csv"]))[['gene', 'sense', 'antisense']]
    finaldata = copy.deepcopy(sample)
    for c in controls.items():
        print("Comparing replicate against", c[0])
        controls[c[0]] = test_against_control(sample, c)
        finaldata = pd.merge(finaldata, controls[c[0]], on='gene', how='inner')

    # Rename columns of sense and antisense with normalized data to sense_normalized and antisense_normalized to keep
    # these values in the database
    finaldata.rename(columns={'sense': 'sense_normalized', 'antisense': 'antisense_normalized'}, inplace=True)

    # Merge the FDR values with the non-normalized data for import in Django and for creating a standalone SVG plot
    # Create a results file that consists of both binomal as well as fisher extacts results
    binom_data = read_csv("".join([sampledir, "sense_vs_antisense_counts.csv"]))[['gene', 'sense', 'antisense', 'fdr']]
    binom_data.rename(columns={'fdr': 'binom_fdr'}, inplace=True)
    finaldata = pd.merge(binom_data, finaldata, on='gene', how='inner')
    write_csv("results.csv", sampledir, finaldata)  # Write to CSV file for Django

    # Create a an 'accumulated' FDR column for labeling and coloring only those genes that are significant against all
    # controls
    finaldata['fdr'] = np.where((
        (finaldata['fcpv_control_1'] <= 0.05) &
        (finaldata['fcpv_control_2'] <= 0.05) &
        (finaldata['fcpv_control_3'] <= 0.05) &
        (finaldata['fcpv_control_4'] <= 0.05)), 0, 1)
    finaldata = finaldata[['gene', 'sense', 'antisense', 'fdr']]
    plot(finaldata, "".join([sampledir, "replicate_vs_controls.svg"]), fisher_plot_title)

if __name__ == "__main__":
    main(sys.argv[1:])
