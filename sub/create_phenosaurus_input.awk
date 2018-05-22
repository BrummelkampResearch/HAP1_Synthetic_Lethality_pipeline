# The columns in results.csv should be as follows:
# $1: gene    
# $2: sense   
# $3: antisense       
# $4: binom_fdr       
# $5: sense_normalized
# $6: antisense_normalized
# $7: fcpv_control_1  
# $8: pv_control_1    
# $9: fcpv_control_2  
# $10: pv_control_2    
# $11: fcpv_control_3  
# $12: pv_control_3    
# $13: fcpv_control_4  
# $14: pv_control_4

BEGIN {
	FS="\t"
	split(var,trim,"|"); #trim[1] = SCREENNAME; trim[2] = Replicate; trim[3] = RelRef
}
{
	Sense=$2
	Antisense=$3
	# Create zero-corrected counts for calculating total and ratio
	if (Sense < 1) {
		SC=1;
	}
	else {
	    SC=Sense;
	}
	if (Antisense < 1) {
		ASC=1;
	}
    else {
	    ASC=Antisense;
	}
    Insertions=SC+ASC
    SR=SC/(SC+ASC)
    print ","trim[1]","trim[2]","$1","$2","$3","$4","$5","$6","$7","$8","$9","$10","$11","$12","$13","$14","Insertions","SR","trim[3]
}
