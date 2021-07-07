#!/bin/bash

# This script contains a bash function to extract and sort gene names from the
# .. transcript_to_genes file created by Kallisto

extract_genes(){
  transcripts_to_genes_file=$1
  output_directory=$2

  cut -f 2,3 $transcripts_to_genes_file | sort -u > $output_directory/genes_only_noheader.tsv
  { echo "ENSEMBL_ID  gene_name" ; cat $output_directory/genes_only_noheader.tsv ; } > $output_directory/ensg_gname.tsv

  rm $output_directory/genes_only_noheader.tsv
}

extract_genes $1 $2
