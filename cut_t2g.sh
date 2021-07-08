#!/bin/bash

# This script contains a bash function to extract and sort gene names from the
# .. transcript_to_genes file created by Kallisto

extract_genes(){
  transcripts_to_genes_file=$1
  output_directory=$2

  cut -f 2,3 $transcripts_to_genes_file | sort -u > $output_directory/ensg_gname.tsv

}

extract_genes $1 $2
