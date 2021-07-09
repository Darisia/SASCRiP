#!/bin/bash

# This bash script contains the code for a function that will
# .. edit the 10xv1 fastq files and rearrange it so it is in the 
# .. correct format for input into < kallisto bus >

edit_10xv1_fastq(){

	IFS=","

	for val in $1;
	do

		filename=`basename $val`

		zcat $val | grep -A3 "1:N:0:0" | sed "/--/d" | bgzip -c > $2/R1_$filename
		zcat $val | grep -A3 "4:N:0:0" | sed "/--/d" | bgzip -c > $2/R2_$filename

	done
}

edit_10xv1_fastq $1 $2