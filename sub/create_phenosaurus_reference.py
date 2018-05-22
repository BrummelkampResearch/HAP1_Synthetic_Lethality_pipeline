#/usr/bin/python3
import sys, getopt, os
import pandas as pd
import numpy as np

def main():
	print("This script is called from normalize.py, do not run it directly")

def build_phenosaurus_reference(normalized_counts, raw_counts, screenname, refname, outdir, replicate):
	raw_counts.rename(columns={'fdr':'binom_fdr', 'p': 'fdr', 'gene':'relgenename'}, inplace=True)
	normalized_counts.rename(columns={'sense': 'sense_normalized', 'antisense': 'antisense_normalized', 'gene':'relgenename'}, inplace=True)
	merged = pd.merge(normalized_counts, raw_counts, on='relgenename')
	merged['replicate'] = replicate
	merged['relscreenname'] = screenname
	merged['pv_control_1'] = merged['pv_control_2'] = merged['pv_control_3'] = merged['pv_control_4'] = 1
	merged['fcpv_control_1'] = merged['fcpv_control_2'] = merged['fcpv_control_3'] = merged['fcpv_control_4'] = 1
	merged['sc'] = np.where(merged.sense<1, 1, merged.sense)
	merged['asc'] = np.where(merged.antisense<1, 1, merged.antisense)
	merged['insertions'] = merged.sc + merged.asc
	merged['senseratio'] = merged.sense/merged.insertions
	merged['relrefname'] = refname
	merged = merged.reindex(columns = ['id','relscreenname','replicate','relgenename','sense', 'antisense' ,'binom_fdr',
									   'sense_normalized' ,'antisense_normalized' ,'fcpv_control_1', 'pv_control_1',
									   'fcpv_control_3', 'pv_control_3', 'fcpv_control_4', 'pv_control_4',
									   'fcpv_control_2' ,'pv_control_2', 'insertions', 'senseratio', 'relrefname'])
	merged.to_csv("".join([outdir, screenname, "_replicate_", str(replicate), "_for_phenosaurus.csv"]), sep=',', index=False)


if __name__ == "__main__":
	main()
