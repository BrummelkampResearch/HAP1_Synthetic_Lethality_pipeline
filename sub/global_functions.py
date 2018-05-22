from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analytical_settings import plot_width, plot_height, plot_dpi, plot_sig_color, plot_nonsig_color, plot_x_label, plot_y_label 

def testinput(datadir, infiles=None, outfiles=None):
	if outfiles is None:
		outfiles = []
	if infiles is None:
		infiles = []
	for i in infiles:
		if not Path("".join([datadir, i])).is_file(): # Check if required files are present
			print(i, "not found. \nPlease check your input. Exiting now.")
			sys.exit(2)
	for o in outfiles: # Check if screen is already analyzed
		if Path("".join([datadir, o])).is_file():
			print(o, "already present, screen already analyzed?\nPlease check your input. Exiting now.")
			sys.exit(2)

def read_csv(file):
	raw_data = pd.read_csv(file, sep='\t')
	return raw_data

def write_csv(filename, outdir, df):
	df.to_csv("".join([outdir, filename]), sep='\t', index=False)

def plot(data, filename, plot_title):
	print("Creating SVG plot", filename, "...", end='')
	data['color'] = np.where(data['fdr']<=0.05, plot_sig_color, plot_nonsig_color)

	# Replace zero's for a 1 for plotting
	data['sense_corrected'] = np.where(data['sense']==0, 1, data['sense'])
	data['antisense_corrected'] = np.where(data['antisense']==0, 1, data['antisense'])

	x = np.log10(data.sense_corrected+data.antisense_corrected) # sense_corrected + antisense_corrected is total counts
	y = (data.sense_corrected)/(data.sense_corrected+data.antisense_corrected)*100
	colors = data.color
	plt.figure(figsize=(plot_width, plot_height), dpi=plot_dpi) # Set size of figure
	plt.scatter(x, y, c=colors, alpha=0.5) # Create a scatter plot
	plt.xlabel(plot_x_label)
	plt.ylabel(plot_y_label)
	plt.ylim([0,100])
	plt.title(plot_title)

	for gene in data[data['fdr']<=0.05].iterrows():
		plt.annotate(gene[1][0], xy = (np.log10(gene[1][5]+gene[1][6]), ((gene[1][5])/(gene[1][5]+gene[1][6]))*100), xytext = (5, 5), textcoords = 'offset points', fontsize = 7)

	plt.savefig(filename, format="svg")
	print("done")

def main():
    print("This script cannot be run directly... bye\n\n")

if __name__ == "__main__":
    main()

