#!/bin/bash

# This bash script contains the code for a function that will
# .. edit the output bustools barcode file and rename it
# .. so it is in the correct format for Seurat input

bus2CR_barcodes(){

	bustools_barcodes_file=$1
	output_directory=$2

	cp $bustools_barcodes_file $output_directory/barcodes.tsv
  gzip $output_directory/barcodes.tsv
}

bus2CR_barcodes $1 $2
