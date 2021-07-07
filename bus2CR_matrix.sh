#!/bin/bash

# This bash script contains the code for a function that will
# .. rearrange the output bustools matrix so that is it in
# .. the correct format for input into Seurat (R)

bus2CR_matrix(){

	# Set the variables
	bustools_matrix_file=$1
	output_directory=$2

	# Remove the three header lines
	sed "/^%/d" $bustools_matrix_file > $output_directory/no_header_bus_matrix.txt

	# Cut out the columns and put it into separate files
	cut -d " " -f 2 $output_directory/no_header_bus_matrix.txt > $output_directory/gene_colmatrix.tsv
	cut -d " " -f 1 $output_directory/no_header_bus_matrix.txt > $output_directory/barcode_colmatrix.tsv
	cut -d " " -f 3 $output_directory/no_header_bus_matrix.txt> $output_directory/UMI_colmatrix.tsv

	# Paste everything together (in the cellranger order) into a new file
	paste  $output_directory/gene_colmatrix.tsv  \
	$output_directory/barcode_colmatrix.tsv \
	$output_directory/UMI_colmatrix.tsv > $output_directory/matrix_noheader.mtx

	# Add the header lines from the CellRanger matrix
	{ echo '%%MatrixMarket matrix coordinate integer general\n%metadata_json: {"format_version": 2, "software_version": "3.1.0"}' ; cat $output_directory/matrix_noheader.mtx ; } > $output_directory/matrix.mtx

	#bgzip the matrix file
	gzip $output_directory/matrix.mtx

	# Remove all the other unnecessary files
	rm $output_directory/no_header_bus_matrix.txt
	rm $output_directory/gene_colmatrix.tsv
	rm $output_directory/barcode_colmatrix.tsv
	rm $output_directory/UMI_colmatrix.tsv
	rm $output_directory/matrix_noheader.mtx
}

bus2CR_matrix $1 $2
