# Pipeline for the analysis of synthetic lethal screens in human haploid cells

## 1. Installation

1. Either clone the GitHub repository or extract the downloaded zip archive at your favorite place
1. Ensure yourself that you have installed the following packages
   * Python 3
   * ktinker (apt-get install python3-tk or yum install python34-tkinter)
   * statsmodels
   * pandas
   * scpipy
   * numpy
1. Optional step if you want the required Python libraries to be installed in a virtual environment rather than
systemwide  
`virtualenv -p python3 venv`  
`source venv/bin/activate`
4. Install the following packages using pip  
`pip install matplotlib numpy pandas`



## Preparing the control dataset ##

* If you have installated the packages in a virtual environment, open the virtual environment by typing  
`source venv/bin/activate`
* For each of the four controls run analyze_sli.sh using the following options  
    `./analyze_sli.sh -C -R=<1-4> -S=sequence_reads.fastq.gz --name=ControlData-HAP1`
* Update the analytical_settings.py configuration file and make sure each 'HAP1_control' points to the corresponding 
sense_vs_antisense_counts.csv file that now resides in the final directory that is listed in settings.conf
* Now run the following command to normalize the counts of the controls:  
`python3 sub/normalize.py --normalize=controls --outdir=<some directory where you want to store the normalized 
controls> --screenname=<screenname> --refname=<name of the reference annotation, same as entered in 
alignment_settings.conf>`  
This will create two files for each replicate:  
  1. `ControlData-HAP1_replicate_1_sense_vs_antisense_counts_normalized.csv`
  1. `ControlData-HAP1_replicate_1_for_phenosaurus.csv`
  * Update the analytical_settings.py configuration file and make sure each 'HAP1_control' points to the corresponding 
sense_vs_antisense_normalized_counts.csv file that now resides in the temp directory.
* You are now ready to analyze your first synthetic lethal screen
