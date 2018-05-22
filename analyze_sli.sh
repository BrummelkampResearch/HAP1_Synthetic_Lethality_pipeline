#!/bin/bash

#####################################################################
# A tool to analyze synthetic lethal screens in human haploid cells
#####################################################################
# authors:
# Elmer Stickel, Vincent Blomen
# Netherlands Cancer Institute
version="1.2 (May 22th, 2018)"
# Commercial use of the software is prohibited without
# prior approval
#####################################################################

# Steps involved in this script:
# 1. Fetch commandline arguments
# 2. Fetch parameters from alignment_settings.conf file
# 3. Pre-run tests
# 4. The actual program
#       4.1. Alignment to genome using bowtie
#       4.2. Annotate insertions to genes using intersectBed
#       4.3. Compare number of sense and antisense insertions using Bionomial test
#	4.4. Normalize number of insertions towards control dataset
#	4.5. Compare sense/antisense ratio's of sample against control dataset
# 5. Create result files for Excel and Django
# 6. Compress output
# 7. Move files back to final directory

# Keep track of running time
starttime=$(date +%s.%N) #Keep track of current date in seconds (to calculate total running time)
#


#################################################
# 1. Fetch command line parameters
#################################################

for i in "$@"
do
case $i in
	-S=*|--seq=*)   # File containing sequence reads
        SEQFILE="${i#*=}"
        shift
	;;
	-N=*|--name=*)  # Name of the screen
        SCREENNAME="${i#*=}"
        shift
	;;
	-R=*|--replicate*) # Name of the replicate (1,2,3...)
	    REPLICATE="${i#*=}"
	    shift
	;;
	-C|--createcontrol) # Create control option
	    cc="yes"
	    shift
	;;
	-h|-?|--help)
        HELP="yes"
        shift
    ;;
	*)
	;; # For any unknown options
esac
done

if [[ $HELP == "yes" ]]
        then
            	cat sub/help # Just cat the help file
                exit
fi

###########################################################
# 2. Fetch parameters from settings files
###########################################################
bowtie=$(awk -F '\t' '$1 ~ /location_bowtie/ {print $2}' alignment_settings.conf)
ref_genome=$(awk -F '\t' '$1 ~ /ref_genome/ {print $2}' alignment_settings.conf)
gene_ref=$(awk -F '\t' '$1 ~ /gene_ref/ {print $2}' alignment_settings.conf)
ref_id=$(awk -F '\t' '$1 ~ /ref_id/ {print $2}' alignment_settings.conf)
trim_reads=$(awk -F '\t' '$1 ~ /trim_reads/ {print $2}' alignment_settings.conf)
if [[ $trim_reads == "yes" ]]
	then
               	keep_reads_from=$(awk -F '\t' '$1 ~ /keep_reads_from/ {print $2}' alignment_settings.conf)
               	trimsize=$(awk -F '\t' '$1 ~ /desired_read_length/ {print $2}' alignment_settings.conf)
fi
core_num=$(awk -F '\t' '$1 ~ /core_num/ {print $2}' alignment_settings.conf)
mem_limit=$(awk -F '\t' '$1 ~ /memory_limit/ {print $2}' alignment_settings.conf)
mismatch_num=$(awk -F '\t' '$1 ~ /mismatch_num/ {print $2}' alignment_settings.conf)
final_dir=$(awk -F '\t' '$1 ~ /final_dir/ {print $2}' alignment_settings.conf)
tmp_dir=$(awk -F '\t' '$1 ~ /tmp_dir/ {print $2}' alignment_settings.conf)

H1C=$(awk -F '= ' '{gsub(/'"'"'/, "", $2)}{if($1 ~ /HAP1_control_1/)print $2}' analytical_settings.conf)
H2C=$(awk -F '= ' '{gsub(/'"'"'/, "", $2)}{if($1 ~ /HAP1_control_2/)print $2}' analytical_settings.conf)
H3C=$(awk -F '= ' '{gsub(/'"'"'/, "", $2)}{if($1 ~ /HAP1_control_3/)print $2}' analytical_settings.conf)
H4C=$(awk -F '= ' '{gsub(/'"'"'/, "", $2)}{if($1 ~ /HAP1_control_4/)print $2}' analytical_settings.conf)


##########################################################
# 3. Pre-run tests (parameters and settings file)
##########################################################
# A bunch of variables we need
ws="\nWrong syntax given\n Run --help or see readme.txt for commandline arguments\n\n"
ff="FOUND\n"
fnf="NOT FOUND\n\tExiting now\n"
tl="\n\tp.s. To change the location of the bowtie/reference and other settings edit settings.conf"
mme="ERROR \n\tnumber of mismatches can only range from 0 - 3\n\tExiting now\n"

## Print welcome logo
cat sub/logo.art
printf "*****\t\t\t\t\t\t\t\t\t\t   *****
*****\t\t\tAnalyis tool for synthetic lethal screens\t\t   *****
*****\t\t\tVersion: $version\t\t\t\t   *****
*****\t\t\t\t\t\t\t\t\t\t   *****
****************************************************************************************\n\n
Performing a few simple checks on the input parameters:
\t---------------------------------\n\tScreenname: $SCREENNAME - replicate: $REPLICATE\n\t---------------------------------\n"

# Test whether user has entered a screenname, replicate number and whether replicate hasn't been analyzed already.
if [[ ($SCREENNAME == '' || $REPLICATE == "" || $SEQFILE == "") && $HELP != "yes" ]]; then
    printf "$ws" && exit
elif [[ " (1 2 3 4[@]) " =~ $REPLICATE ]]; then # If so, check if given replicate number is a valid number (ie. 1, 2, 3 or 4)
    printf "\tChecking screenname $SCREENNAME..."
    if [[ -d $tmp_dir$SCREENNAME"_output" ]]; then
        printf "ERROR\n\n --> Directory for $SCREENNAME already present in temp directory <--\nPlease check the log files why it's there\nQuiting now\n\n." && exit
    elif [[ -d $final_dir$SCREENNAME"_output" ]]; then
        SP=1
        printf "found in final dir\n\tChecking replicates..."
        if [ -d $final_dir$SCREENNAME"_output/replicate_"$REPLICATE ]; then
            printf "ERROR\n\n --> replicate $REPLICATE analyzed already <--\nQuiting now\n\n." && exit
        else
            printf "OK\n"
        fi
    else
        printf "\tOK\n"
        SP=0
    fi
    OUT=$tmp_dir$SCREENNAME"_output/replicate_"$REPLICATE
    mkdir -p $OUT/$ref_id && touch $OUT/genome_align_log.txt # Create log file
else
    printf "\nSorry! Wrong syntax, replicate can only be 1, 2, 3 or 4 and not "$REPLICATE"\n\n" && exit
fi


# Write some stuff to log file and create a second log file (annotation_analysis_log.txt)
printf "\n Version of analysis pipeline used: $version" >> $OUT/genome_align_log.txt | tee -a $OUT/$ref_id/annotation_analysis_log.txt # Write verion to logfiles
printf "\n---------------- screenname: $SCREENNAME ----------------\n Date:" >> $OUT/genome_align_log.txt # Write date to logfile
date >> $OUT/genome_align_log.txt && printf "\n" >> $OUT/genome_align_log.txt
printf "\tThe logfiles will be written to: \n\t\t--> $OUT/genome_align_log.txt <-- and \n\t\t--> $OUT/$ref_id/annotation_analysis_log.txt <--\n" | tee -a $OUT/genome_align_log.txt
cp $OUT/genome_align_log.txt $OUT/$ref_id/annotation_analysis_log.txt
printf "\tTrimming reads?" | tee -a $OUT/genome_align_log.txt # Check for read trimming and write to log
if [[ $trim_reads == "yes" ]]; then printf "\tYes, new length is $trimsize bases starting from base $keep_reads_from\n" | tee -a $OUT/genome_align_log.txt
else printf "\tNo\n" | tee -a $OUT/genome_align_log.txt
fi

# Check for presence of fastq file of replicate
printf "\tFastq file of synthetic lethal replicate: ... " | tee -a $OUT/genome_align_log.txt # Check for fastq input file and write to log
if [ -f $SEQFILE ]; then printf "$ff" | tee -a $OUT/genome_align_log.txt; else printf "$fnf" && exit; fi

# Check for presence of normalized controls
if [[ $cc != "yes" ]]; then
    i=1
    for c in $H1C $H2C $H3C $H4C
    do
        printf "\tChecking control $i... " | tee -a $OUT/$ref_id/annotation_analysis_log.txt
        if [ -f $c ] ; then printf "$ff"  | tee -a $OUT/$ref_id/annotation_analysis_log.txt; else printf "$fnf" | tee -a $OUT/$ref_id/annotation_analysis_log.txt; exit; fi
        i=$(($i+1))
    done
else
    printf "\n\t*** Running in special mode to prepare controls ***\n" | tee -a $OUT/$ref_id/annotation_analysis_log.txt
fi

# Check bowtie human genome reference files and write to log
printf "\tHuman genome reference files: $ref_genome... " | tee -a $OUT/genome_align_log.txt
if [[ $(ls -l $ref_genome* | wc -l)  -gt 0 ]]; then printf "$ff" | tee -a $OUT/genome_align_log.txt; else printf "$fnf"; printf "$tl"; exit; fi

# Check of Bowtie is present and write to log
printf "\tLocation of bowtie: $bowtie... " | tee -a $OUT/genome_align_log.txt
if [ -f $bowtie ]; then printf "$ff" | tee -a $OUT/genome_align_log.txt; else printf "$fnf" && printf "$tl"; exit; fi

# Check number of allowed mismatches and write to log
printf "\tChecking number of mismatches:..." | tee -a $OUT/genome_align_log.txt
if [[ $mismatch_num>=0 && $mismatch_num<4 ]]; then printf "($mismatch_num) OK " | tee -a $OUT/genome_align_log.txt; else printf "($mismatch_num) $mme" && printf "$tl"; exit; fi
printf "\n\t-----------------------
Everyting seems fine\n\n"

# Copy the settings file to the output directory, always  useful to have
cp alignment_settings.conf $OUT/alignment_settings.conf
cp analytical_settings.conf $OUT/$ref_id/analytical_settings.conf



####################################
# 4. Here the actual program starts
####################################
# Start by creating a symlink to the fastq file and create a variable (firstinput) that can be used to open a stream of the lines in the sequence file
rep=replicate
if [[ $SEQFILE =~ .gz ]]; then	# If the low file is compresses (.gz) then create a symlink with .gz behind it
	ln -s $SEQFILE $OUT/$rep.fastq.gz
	firstinput="zcat $OUT/$rep.fastq.gz"
else
	ln -s $SEQFILE $OUT/$rep.fastq
	firstinput="cat $OUT/$rep.fastq"
fi

##########################################################
# 4.1. Preprocessing
##########################################################




# Now extract only the reads from the file, the quality info is not used
# In case the reads are trimmed as well, only the desired piece of the reads is used.
# This is the fastq format that is expected:
#
# line 1: @HWI-ST867:213:C3R0EACXX:4:1101:8930:1999 1:N:0:
# line 2: NCTGCAGCACAAAAGACGTGATGACTCTTCTCCAGTGAGTCAAGGGGCCGT
# line 3: +
# line 4: 4:BDDDEHHGGHIIIIIHIIIIIIIIIIIIIIHIDGGGEHGHEBDFHHI;
#
# In short, the reads are every 4th line, starting with the 2nd line.
#
# Extracting the reads and trimming is done using AWK.
# The info from where the reads should be kept and its length is parsed as single variable (b) and separared by a pipe character
# In the BEGIN statement of the awk script, the split command spilts the variable b into a array called trim with two two elements, the first element holds the startposition from where to keep the read and the second element holds the lengtgh of the read
# For all even-line-numbers (these contain the actual read and quality info) trim the length
 # Why check whether trimming is needed? To speed up the process if it isn't, otherwise an identical copy of the fastq file is generated without any need for it.
if [[ $trim_reads == "yes" ]]; then
	printf "Trimming and extracting reads from fastq file..."
	$firstinput | awk -v b="$keep_reads_from|$trimsize" '
		BEGIN{
			split(b,trim,"|");
		}
        {
			if(NR%4==2){
				print ">\n"substr($0, trim[1], trim[2])
			}
        }' > $OUT/$rep.fa
	printf "done\n"
else
	printf "Extracting reads (without trimming them) from fastq file..."
	$firstinput | awk '
		{
			if(NR%4==2){
				print ">\n"$0
			}
		}' > $OUT/$rep.fa
	printf "done\n"
fi


##########################################################
# 4.2. Alignment to genome using bowtie
##########################################################


# This is where bowtie is executed, and will be i-times according to the number of mismatches allowed for (max 3, more mismatches are not allowed by bowtie, nor does it make any sense)
# The rational for running it for each number of mismatches individually is that we require bowtie to align uniquely to the genome. So for
# For example, take the following sequence AATTGG
# If we allow one mismatch, it allign can allign to following sequencing on the genome AATTGG and AATTGC. Ie. bowtie cannot assign the reads to a unique genomic location and hence discards this reads because of the unique-requirment.
# However, if we now run bowtie without allowing any mismatches, the reads now can align uniquely to the genome and the read is kept. Thus, if we allow 1 mismatch, we run bowtie once with 0 mismatches and then with 1 mismatch
# Fet the read length, needed for Bowtie command. Just read the first line of the file containing the reads and determine the length of it using awk
printf "Please take a break while the program is mapping the reads to the genome ...\n"
printf "\nAligments in $rep:\n" >> $OUT/genome_align_log.txt

bwtinfile=$OUT/$rep.fa
for ((i=$mismatch_num;i>-1;i-=1)); do   # Loop in reverse order over number of mismatches
	if [[ $i -lt $mismatch_num ]]; then # If not first iteration: take previously suppressed reads as input
		bwtinfile=$bwtsuppressedfile
        fi
	bwtoutfile=$OUT/$rep.mapped_mm$i.bwt
	bwtsuppressedfile=$OUT/$rep.suppressed_mm$i.bwt

	printf "\n-------------------------------\nAlligning with $i mismatches\n-------------------------------\nIntermediate bowtie-statistics\n"
	$bowtie -p 30 --best -v $i -m 1 -f $ref_genome $bwtinfile $bwtoutfile --max $bwtsuppressedfile
	# The output file looks like (it does not have a header)
        # 469     -       chr22   30347844        GTTGCTGGTTCAAGTGATTTATTACAAGTAAAGGTTTTTTTTTTTTTTTT      IIIIIIIIIIIIIIIIIIII$
        # 462     -       chr5    54795872        ATCATATATGGTATAATTTGATTTCTCTAAAAATTTTTTTTTTTTTTTTT      IIIIIIIIIIIIIIIIIIII$
done
printf "\n-------------------------------\n"


# Concatenate the output from the multiple bowtie alignments (with n mismatches) and extract strand, chrom and pos.
printf "Merging bowtie results into single file..."
cat $OUT/$rep.mapped_mm*.bwt  > $OUT/$rep.mapped_all.bwt
printf "done\n"

# Now sort on strand, chromosome and location on chromosome and remove duplicates. We explicitely do not sort on the sequence as allowing one or more
# mismatches can yield a slightly different sequence to allign on the same location. We don't want to count them as unique insertions.
printf "Extracting distinct reads..."
# New algorithm using associative array, first chrom, then strand and then pos. because if chr then pos there can be mixup
# chr17:84782 -> chr1784782 then equals chr1:784782 --> chr1784782! Not good!
awk '!seen[$3$2$4]++' $OUT/$rep.mapped_all.bwt > $OUT/$rep.mapped_unique.bwt
printf "Number of unique reads in $f:" | tee -a $OUT/genome_align_log.txt
wc -l $OUT/$rep.mapped_unique.bwt | awk '{print $1}' | tee -a $OUT/genome_align_log.txt
printf "done\n"
# Alignment done

##########################################################
# 4.3. Annotate insertions to genes using intersectBed
##########################################################

# Alllocated to to genes using the intersect tool from BEDtools
# As a starting point, we use the combined, sorted and uniqued bowtie output. For these steps we only use column 2 (strand) 3 (chromosome), 4 (position) and 5 (sequence)
printf "Processing sample for gene annotation by intersectBed..."
awk -F"\t" '
{
    if ($2=="+"){
		# If + strand
	    print $3"\t"$4"\t"$4+1"\t"$5"\t"$2
		# This looks like (in case of low file)
		# chr1		234749136	234749137	GTTTCCAGAGCTTACTCCAGTCAAGCAGCTATTAACACACGGAAGCCCTG	+
		# In the 3th column we add 1 to the position because of we are using the USCS refseq genome browser representation which has a one-base end (see https://genome.ucsc.edu/FAQ/FAQtracks.html)
	}
	else{
		# if - strand
		print $3"\t"$4+length($5)-1"\t"$4+length($5)"\t"$5"\t"$2
		# This looks like (in case of low file, length of GTTGCTGGTTCAAGTGATTTATTACAAGTAAAGGTTTTTTTTTTTTTTTT is 50)
		# chr22		30347893	30347894	GTTGCTGGTTCAAGTGATTTATTACAAGTAAAGGTTTTTTTTTTTTTTTT	-
	}
}' $OUT/$rep.mapped_unique.bwt > $OUT/$rep.input4enrichment.bed
printf "done\n"

# Now run intersectbed with -wo option
printf "Mapping insertions to genes using intersectbed..."
intersectBed -a $OUT/$rep.input4enrichment.bed -b $gene_ref -wo >| $OUT/$ref_id/$rep.inputtocoupletogenes+gene_all.txt
# The output file looks like
# chr10 100001181       100001182       AATATAGTCATTTCTGTTGATTGATTTATTATATTGGAGAATGAAAGTTT	+       chr10   99894381        100004654       R3HCC1L 0       +       1
# chr10 100010840       100010841	ATGGCTGATTCCACAGTGGCTCGGAGTTACCTGTGTGGCAGTTGTGCAGC      -       chr10   100007443       100028007       LOXL4   0       -       1
printf "done\n"

# Now the reads have been mapped to genes we separate the sense integrations from the antisense integrations. In the remainder of this program only the sense integrations are used but for legacy purposes we still keep the antisense integrations
# If integration is sense than strand (col 4) and orientation (col 11) match
# The format of the output file is as follows:
# chr1    10001144        10001145        LZIC
printf "Extracting sense integrations..."
awk -F"\t" '{ if($5==$11) print $1"\t"$2"\t"$3"\t"$9 }' $OUT/$ref_id/$rep.inputtocoupletogenes+gene_all.txt | sort --parallel=$core_num --buffer-size=$mem_limit > $OUT/$ref_id/$rep.gene+screen_sense.txt
printf "done\n"
# If integration is antisense,  then strand (col 4) and orientation (col 11) do not match
printf "Extracting antisense integrations..."
awk -F"\t" '{ if($5!=$11) print $1"\t"$2"\t"$3"\t"$9 }' $OUT/$ref_id/$rep.inputtocoupletogenes+gene_all.txt | sort --parallel=$core_num --buffer-size=$mem_limit > $OUT/$ref_id/$rep.gene+screen_antisense.txt
printf "done\n"

##########################################################
# 5. Call Python scripts for stat. tests
##########################################################

python3 sub/binomial_test.py --data=$OUT/$ref_id/
if [[ $cc != "yes" ]]; then
    python3 sub/normalize.py --normalize=$OUT/$ref_id/
    python3 sub/test_replicate_vs_control.py --replicate=$OUT/$ref_id/

    # Create input file for Phenosaurus
    head -n 1 $OUT/$ref_id/results.csv | awk -F '\t' '{print "id,relscreenname,replicate,relgenename,"$2","$3","$4","$5","$6","$7","$8","$9","$10","$11","$12","$13","$14",insertions,senseratio,relrefname"}' > $OUT/$ref_id/for_phenosaurus_head.csv
    awk -F'\t' 'NR>1{print $0}' $OUT/$ref_id/results.csv | awk -v var="$SCREENNAME|$REPLICATE|$ref_id" -f sub/create_phenosaurus_input.awk > $OUT/$ref_id/for_phenosaurus_content.csv
    cat $OUT/$ref_id/for_phenosaurus_head.csv $OUT/$ref_id/for_phenosaurus_content.csv > $OUT/$ref_id/for_phenosaurus.csv
    rm $OUT/$ref_id/for_phenosaurus_content.csv $OUT/$ref_id/for_phenosaurus_head.csv
fi
##########################################################
# 6. Compressing and moving files to their final directory
##########################################################

printf "\nCompressing and moving files to their final directory\n"

bash sub/compress_output.sh $OUT
bash sub/compress_gene_ref_output.sh $OUT/$ref_id/

endtime=$(date +%s.%N)
printf "*****************************"
printf "\nTotal runtime hh:mm:ss: " | tee -a $OUT/genome_align_log.txt
runtime=$(python3 -c "import datetime; print(str(datetime.timedelta(seconds=(${endtime} - ${starttime}))))")
printf $runtime"\n" | tee -a $OUT/genome_align_log.txt

if [[ $SP == 0 ]]; then
    mkdir $final_dir$SCREENNAME"_output/"
fi
mv $OUT $final_dir$SCREENNAME"_output/"
rm -r $tmp_dir$SCREENNAME"_output"
