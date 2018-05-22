###############################
# Settings for normalization
###############################
# Location of files with HAP1 control datasets
# When running the normalize.py with the '-n controls' arguments the next four lines represent the non-normalized controls,
# when normalizing an experiment the next four lines should represent the normalized-controls
HAP1_control_1 = '/media/data/analyzed_data/synthetic_lethal_screens/ControlData-HAP1_output/replicate_1/GRCh38-NCBI_RefSeq/sense_vs_antisense_counts_normalized.csv'
HAP1_control_2 = '/media/data/analyzed_data/synthetic_lethal_screens/ControlData-HAP1_output/replicate_2/GRCh38-NCBI_RefSeq/sense_vs_antisense_counts_normalized.csv'
HAP1_control_3 = '/media/data/analyzed_data/synthetic_lethal_screens/ControlData-HAP1_output/replicate_3/GRCh38-NCBI_RefSeq/sense_vs_antisense_counts_normalized.csv'
HAP1_control_4 = '/media/data/analyzed_data/synthetic_lethal_screens/ControlData-HAP1_output/replicate_4/GRCh38-NCBI_RefSeq/sense_vs_antisense_counts_normalized.csv'


###############################
# Settings for plots
###############################
plot_width = 20
plot_height = 14
plot_dpi = 80
plot_sig_color = 'red'
plot_nonsig_color = 'black'
plot_x_label = r'$\log_10\left(total\ integrations\right)$'
plot_y_label = r'% Sense insertions'
binom_plot_title = 'Ratio sense versus antisense integrations'
fisher_plot_title = 'Normalized ratio sense versus antisense integrations'

###############################
# Advanced options
###############################
# --> ALWAYS MAKE SURE THE CONTROLS ARE NORMALIZED THE SAME WAY YOU NORMALIZE YOUR EXPERIMENT!!!! <--
# (ie. with or without cauchy.R compatibility, please refer to documentation)
cauchy_compatibility = False
grouping = 500 # Desired groupsize for normalization, 500 recommended groupsize
