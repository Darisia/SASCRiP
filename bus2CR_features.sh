#!/bin/bash

# This bash script contains the code for a function that will
# .. edit the output bustools features file and rearrange it
# .. so it is in the correct format for Seurat input

bus2CR_features(){

	bustools_features_file=$1
	output_directory=$2

	cp $bustools_features_file $output_directory/features.tsv
	gzip $output_directory/features.tsv
}

bus2CR_features $1 $2
