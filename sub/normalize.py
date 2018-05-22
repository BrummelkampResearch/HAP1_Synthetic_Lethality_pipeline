#!/usr/bin/python3
import sys, getopt, os
import pandas as pd
import numpy as np

from analytical_settings import HAP1_control_1, HAP1_control_2, HAP1_control_3, HAP1_control_4
from analytical_settings import grouping
from analytical_settings import cauchy_compatibility

from global_functions import read_csv, write_csv
from create_phenosaurus_reference import build_phenosaurus_reference


##################################
# Function 'normalize_to_median'
##################################
# Arguments: reference <pandas DataFrame with accumulated counts sense and antisense>
#            sample    <pandas DataFrame with sense and antisense counts for
#                       sample that requires normalization>
# Returns:   ns        <pandas DataFrame with identical structure to sample but with normalized counts>
##################################
def normalize_to_median(raw, reference):
    # Calculate sense/total ratio for reference dataset and store in 'ref_ratio' column in DataFrame
    reference = reference.assign(ref_ratio=(reference['sense_reference'] + 1)
                                           / (reference['sense_reference'] + reference['antisense_reference'] + 2))
    normalized_data = pd.DataFrame(columns=['gene'])
    if cauchy_compatibility:
        print("original Cauchy.R algorithm...", end='')
    else:
        print("latest algorithm...", end='')
    subset = ['gene', 'sense', 'antisense']
    s = pd.merge(reference, raw, on='gene', how='outer')

    # Split sample(s) into part of sample that can be normalized (sn) and the part the cannot be normalized (scn)
    sn = s[((s['sense'] + s['antisense']) > 0) & (
    (s.sense_reference + s.antisense_reference) > 0)]  # meaning counts in reference and sample
    scn = s[((s['sense'] + s['antisense']) > 0) & (
    (s.sense_reference + s.antisense_reference) == 0)]  # meaning counts in sample but not in reference

    # Calculate sense-ratio for current sample and sort DataFrame based on sense-ratio of reference set
    sn = sn.assign(ratio=(s[subset[1]] + 1) / (s[subset[1]] + s[subset[2]] + 2))
    sn = sn.sort_values(by='ref_ratio')

    # Create empty dataframe for the normalized sample (ns)
    ns = pd.DataFrame(columns=subset)

    if not cauchy_compatibility:

        # Calculate optimal grouping value based on requested grouping (to create equal groups and hence identical normalization for every gene)
        divider = np.round(sn.shape[0] / grouping)
        groupsize = (sn.shape[0] / divider)

        # Group sorted DataFrame in equally sized groups and loop over these groups
        for g, df in sn.groupby(np.arange(len(sn)) // groupsize):
            ref_median = df['ref_ratio'].median()
            sample_median = df['ratio'].median()
            df = df.assign(norm_sense=np.where(df['ratio'] <= sample_median
                                               , np.where(
                    np.round((df['ratio'] / (sample_median) * ref_median * (df[subset[1]] + df[subset[2]]))) > (
                    df[subset[1]] + df[subset[2]])  # ie. the normalized sense count is higer than the total counts
                    , (df[subset[1]] + df[subset[2]])
                    , np.round((df['ratio'] / (sample_median) * ref_median * (df[subset[1]] + df[subset[2]]))))
                                               , np.where(np.round(((1 - (1 - df['ratio']) / (1 - sample_median) * (
                1 - ref_median)) * (df[subset[1]] + df[subset[2]]))) > (df[subset[1]] + df[subset[2]])
                                                          , (df[subset[1]] + df[subset[2]])
                                                          , np.round(((1 - (1 - df['ratio']) / (1 - sample_median) * (
                    1 - ref_median)) * (df[subset[1]] + df[subset[2]]))))))

            df = df.assign(norm_antisense=(df[subset[1]] + df[subset[2]] - df['norm_sense']))
            df[['norm_sense', 'norm_antisense']] = df[['norm_sense', 'norm_antisense']].astype(np.int32)
            df = df[['gene', 'norm_sense', 'norm_antisense']]
            df.columns = subset
            ns = ns.append(df)  # Append to DataFrame with normalized counts of previous group

    else:
        # If cauchy_compatibility==True, grouping of genes will be performed according to the previous version of
        # the pipeline as in 'cauchy.R'

        # Calculate the number of groups
        divider = int(np.ceil(sn.shape[0] / grouping))  # The total number of groups
        genes = sn.shape[0]

        i = 1  # Current index of input dataframe
        for g in range(0, divider):  # Loop over the groups of genes
            df = pd.DataFrame()  # Create temporary dataframe for current group
            for gi in range(0, grouping):  # Append the genes to the current group
                i = gi + g * grouping  # DataFrame Index
                if i <= genes:  # For the last group
                    df = df.append(sn[i:i + 1])
            ref_median = df['ref_ratio'].median()
            sample_median = df['ratio'].median()
            df = df.assign(norm_sense=np.where(df['ratio'] <= sample_median
                                               , np.where(
                    np.round((df['ratio'] / (sample_median) * ref_median * (df[subset[1]] + df[subset[2]]))) > (
                    df[subset[1]] + df[subset[2]])  # ie. the normalized sense count is higer than the total counts
                    , (df[subset[1]] + df[subset[2]])
                    , np.round((df['ratio'] / (sample_median) * ref_median * (df[subset[1]] + df[subset[2]]))))
                                               , np.where(np.round(((1 - (1 - df['ratio']) / (1 - sample_median) * (
                1 - ref_median)) * (df[subset[1]] + df[subset[2]]))) > (df[subset[1]] + df[subset[2]])
                                                          , (df[subset[1]] + df[subset[2]])
                                                          , np.round(((1 - (1 - df['ratio']) / (1 - sample_median) * (
                    1 - ref_median)) * (df[subset[1]] + df[subset[2]]))))))

            df = df.assign(norm_antisense=(df[subset[1]] + df[subset[2]] - df['norm_sense']))
            df = df[['gene', 'norm_sense', 'norm_antisense']]
            df.columns = subset
            ns = ns.append(df)  # Append to DataFrame with normalized counts of previous group

    # Merge the normalized sample with the datapoint of the sample that could not be normalized because of absense
    # of these genes in the reference
    ns = ns.append(scn[subset])
    normalized_data = pd.merge(ns, normalized_data, on='gene', how='outer')
    normalized_data[['sense', 'antisense']] = normalized_data[['sense', 'antisense']].astype(np.int32)

    return normalized_data


def process_controls(controls):
    controls = controls.copy()
    reference_cols = ['gene', 'sense_reference', 'antisense_reference']
    reference = pd.DataFrame(columns=reference_cols)
    for c, f in controls.items():  # c is control, f = file
        controls[c] = pd.read_csv(f, sep='\t')[['gene', 'sense',
                                                'antisense']]  # Read the file and overwrite the value of the name of the with the actual data
        reference = pd.merge(controls[c], reference, on='gene', how='outer').fillna(0)
        reference = reference.assign(
            sense_reference=(reference[['sense', 'sense_reference']].sum(axis=1)).astype(np.int32))
        reference = reference.assign(
            antisense_reference=(reference[['antisense', 'antisense_reference']].sum(axis=1)).astype(np.int32))
        reference.drop(['sense', 'antisense'], axis=1, inplace=True)
    return controls, reference


def main(argv):
    syntax = 'normalize.py -n <controls or directory of replicate> -o <directory where normalized control will be saved>'
    try:
        opts, args = getopt.getopt(argv, "hnosr:", ["normalize=", "outdir=", "screenname=", "refname="])
    except getopt.GetoptError:
        print(syntax)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(syntax)
            sys.exit()
        elif opt in ("-n", "--normalize"):
            data = arg
        elif opt in ("-o", "--outdir"):
            outdir = arg
        elif opt in ("-s", "--screenname"):
            screenname = arg
        elif opt in ("-r", "--refname"):
            refname = arg
        else:
            assert False, "unhandled option"

    # Create dictionary to store the names of the controls (keys) and their associated records (values)
    control_files = {'HAP1_control_1': HAP1_control_1, 'HAP1_control_2': HAP1_control_2, 'HAP1_control_3': HAP1_control_3,
                'HAP1_control_4': HAP1_control_4}
    # Read in controls and calculate aggregate reference counts
    controls, reference = process_controls(control_files)

    if data == 'controls':
        for c, d in controls.items():  # c = control, d = data
            print("Normalizing control using ", end='')
            norm = normalize_to_median(d, reference).sort_values(by='gene')
            write_csv(''.join([screenname, "_replicate_", str(int(c[-1:])), "_sense_vs_antisense_counts_normalized.csv"]), os.path.join(outdir, ''), norm) # Create file for normalization of samples
            print(control_files[c])
            build_phenosaurus_reference(norm, read_csv(control_files[c]), screenname, refname, outdir, int(c[-1:]))
            print("done")
    else:
        data = os.path.join(data, '')
        raw = read_csv(''.join([data, 'sense_vs_antisense_counts.csv']))[['gene', 'sense', 'antisense']]
        print("Normalizing sample using ", end='')
        norm = normalize_to_median(raw, reference).sort_values(by='gene')
        write_csv('sense_vs_antisense_counts_normalized.csv', data, norm)
        print("done")


if __name__ == "__main__":
    main(sys.argv[1:])
