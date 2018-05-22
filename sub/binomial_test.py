#!/usr/bin/python3

import sys, getopt, os
import pandas as pd
import numpy as np
import scipy.stats  # For binomial test
import statsmodels.sandbox.stats.multicomp  # For Bejamini Hochberg FDR correction

from global_functions import write_csv, plot, testinput
from analytical_settings import binom_plot_title


# bionomtest reads the files that represent the sense and antisense insertions in genes (as generation by upstream scripts)
# the number of unique sense and antisense insertions are calculated for each gene and a bionomial test is applied on it
# FDR correct (Benajamini Hochberg) and everything is saved in a dataframe. Dataframe also includes columns for
# creating a plot (no zero-values as well as a color column based on FDR corrected p-value)
#
# Arguments: directory with trailing slash
# Returns: pd.DataFrame
#
def binomtest(datadir):
    # File locations of sense and antisense insertions of the sample population
    sensefile = '/'.join([datadir, 'replicate.gene+screen_sense.txt'])
    antisensefile = '/'.join([datadir, 'replicate.gene+screen_antisense.txt'])
    columnnames = ['chromosome', 'start', 'end', 'gene']

    print("Building table...", end='')
    # Read in files, create a pandas dataframe and collapse on genename
    sense_df = pd.read_csv(sensefile, names=columnnames, sep='\t').groupby(
        'gene').size().to_frame().reset_index().rename(columns={
        0: 'sense'})  # although usecols=['gene'] is more elegant, its much slower than reading in the whole file and discarding what is not needed...
    antisense_df = pd.read_csv(antisensefile, names=columnnames, sep='\t').groupby(
        'gene').size().to_frame().reset_index().rename(columns={0: 'antisense'})

    # Combine sense and antisense in a single dataframe
    combined = pd.merge(sense_df, antisense_df, on='gene', how='outer')

    # Replace NaN with zero's for binomial test
    combined['sense'] = combined['sense'].fillna(0).astype(int)
    combined['antisense'] = combined['antisense'].fillna(0).astype(int)

    print("done\nPerforming Binomial tests...", end='')
    # Sum sense and antisense integration to get n-trials or binom test
    combined['p'] = combined[['sense', 'antisense']].apply(scipy.stats.binom_test, axis=1,
                                                           raw=True)  # Sense is considered as a success, this is a two-sided test see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binom_test.html
    print("done\nApplying Multiple Testing correction using Benjamini Hochberg method...", end='')
    combined['fdr'] = \
    statsmodels.sandbox.stats.multicomp.multipletests(combined.p, alpha=0.05, method='fdr_bh', returnsorted=False)[1]
    print("done")
    return (combined)



def main(argv):
    syntax = 'bionomial_test.py -d <dir>'
    datadir = ''
    try:
        opts, args = getopt.getopt(argv, "hd:", ["datadir="])
    except getopt.GetoptError:
        print(syntax)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(syntax)
            sys.exit()
        elif opt in ("-d", "--datadir"):
            datadir = arg

    print("\nPerforming Fisher-Exacts tests against controls\n")
    # Validate input
    datadir = os.path.join(datadir, '')
    testinput(datadir, ['replicate.gene+screen_sense.txt', 'replicate.gene+screen_antisense.txt'], ['sense_vs_antisense_counts.csv', 'binom_plot.svg'])

    data = binomtest(datadir)
    write_csv('sense_vs_antisense_counts.csv', datadir, data)
    plot(data[['gene', 'sense', 'antisense', 'fdr']], "".join([datadir, "binom_plot.svg"]), binom_plot_title)


if __name__ == "__main__":
    main(sys.argv[1:])
