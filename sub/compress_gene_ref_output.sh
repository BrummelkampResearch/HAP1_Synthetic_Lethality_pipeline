#!/bin/bash
#
# Script to compress output results from taqcount.sh
# First argument $1 is parsed into script by taqcount.sh
#
# Create an array of the filenames to be compressed

declare -a tbc=(
replicate.inputtocoupletogenes+gene_all.txt
replicate.gene+screen_antisense.txt
replicate.gene+screen_sense.txt
)

#Check whether pigz is installed, used it if present otherwise fallback to gzip
#pigz uses multicore
if hash pigz 2>/dev/null; then
	core_num=$(awk -F '\t' '$1 ~ /core_num/ {print $2}' alignment_settings.conf) #Fetch number of cores from settings.conf
	zipper="pigz -p $core_num"
	printf "\npigz found. Compressing files on $core_num cores\n"
else
	zipper="gzip"
	prinft "\nDid not find pigz on your machine. Will used gzip instead. Consider installing pigz to facilitate compression on multiple cores\n"
fi

for f in $1/*
do
	g=${f##*/}
	if [[ " ${tbc[@]} " =~ " ${g} " ]];	then
		printf "Compressing $g...\n"
		pv $f | $zipper > $f.gz
		rm $f
	fi
done
if [[ $2 -eq 0 ]]
	then
	printf ""
else
	printf "Compressing $2\n"
	pv $2 | $zipper > $2.gz
	rm $2
fi
if [[ $2 -eq 0 ]]
	then
	printf ""
else
	printf "Compressing $3\n"
	pv $3 | $zipper > $3.gz
	rm $3
fi