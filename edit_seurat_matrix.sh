# This bash script contains the code to edit the matrix accordingly

# Create the bash function that will convert the matrix and indices into the
# .. correct format if needed

# Set arguments
matrix_file=$1
gene_file=$2
barcode_file=$3
output_directory=$4
gene_srt_name=$5
barcode_srt_name=$6
matrix_srt_name=$7
gene_zipped=$8
barcode_zipped=$9
matrix_zipped=${10}
matrix_srt_format=${11}
gene_srt_format=${12}
t2g_file=${13}

# If the gene index file is not in the correct format (gene names) - fix it
if [ $gene_srt_format == 'False' ]
then
    while read -r ENSG
    do
      touch $output_directory/features_gene_names.tsv
      grep -w $ENSG $t2g_file | cut -f 3 | sort -u >> $output_directory/features_gene_names.tsv
    done < $gene_file
    # paste the gene file and the corresponding hgnc symbol file together
    paste $gene_file $output_directory/features_gene_names.tsv > $output_directory/features.tsv
    gzip $output_directory/features.tsv
else
  # If statement to check if the gene file is named properly
  if [ $gene_zipped == 'True' ]
  then
    if [ $gene_srt_name == 'False' ]
    then
      cp $gene_file $output_directory/features.tsv.gz
    fi
  else
    if [ $gene_srt_name == 'False' ]
    then
      cp $gene_file $output_directory/features.tsv
      gzip $output_directory/features.tsv
    else
      gzip $output_directory/features.tsv
    fi
  fi
fi

# If statements to check if the barcodes file is named properly
if [ $barcode_zipped == 'True' ]
then
  if [ $barcode_srt_name == 'False' ]
  then
    cp $barcode_file $output_directory/barcodes.tsv.gz
  fi
else
  if [ $barcode_srt_name == 'False' ]
  then
    cp $barcode_file $output_directory/barcodes.tsv
    gzip $output_directory/barcodes.tsv
  else
    gzip $output_directory/barcodes.tsv
  fi
fi

# If statemnts to check if the matrix file is in the correct arrangement and is named properly
if [ $matrix_zipped == 'False' ]
then
  if [ $matrix_srt_format == 'False' ]
  then
    sed "/^%/d" $matrix_file > $output_directory/no_header_bus_matrix.txt
    cut -d " " -f 2 $output_directory/no_header_bus_matrix.txt > $output_directory/gene_colmatrix.tsv
    cut -d " " -f 1 $output_directory/no_header_bus_matrix.txt > $output_directory/barcode_colmatrix.tsv
    cut -d " " -f 3 $output_directory/no_header_bus_matrix.txt > $output_directory/UMI_colmatrix.tsv
    paste $output_directory/gene_colmatrix.tsv $output_directory/barcode_colmatrix.tsv $output_directory/UMI_colmatrix.tsv > $output_directory/matrix_noheader.mtx
    { echo '%%MatrixMarket matrix coordinate integer general' ; cat $output_directory/matrix_noheader.mtx ; } > $output_directory/matrix.mtx
    gzip $output_directory/matrix.mtx
    # Remove all the other unnecessary files
  	rm $output_directory/no_header_bus_matrix.txt
  	rm $output_directory/gene_colmatrix.tsv
  	rm $output_directory/barcode_colmatrix.tsv
  	rm $output_directory/UMI_colmatrix.tsv
  	rm $output_directory/matrix_noheader.mtx
  else
    if [ $matrix_srt_name == 'False' ]
    then
      cp $matrix_file $output_directory/matrix.mtx
      gzip $output_directory/matrix.mtx
    else
      gzip $output_directory/matrix.mtx
    fi
  fi
else
  if [ $matrix_srt_format == 'False' ]
  then
    zcat < $matrix_file | sed "/^%/d" > $output_directory/no_header_bus_matrix.txt
    cut -d " " -f 2 $output_directory/no_header_bus_matrix.txt > $output_directory/gene_colmatrix.tsv
    cut -d " " -f 1 $output_directory/no_header_bus_matrix.txt > $output_directory/barcode_colmatrix.tsv
    cut -d " " -f 3 $output_directory/no_header_bus_matrix.txt > $output_directory/UMI_colmatrix.tsv
    paste $output_directory/gene_colmatrix.tsv $output_directory/barcode_colmatrix.tsv $output_directory/UMI_colmatrix.tsv > $output_directory/matrix_noheader.mtx
    { echo '%%MatrixMarket matrix coordinate integer general' ; cat $output_directory/matrix_noheader.mtx ; } > $output_directory/matrix.mtx
    gzip $output_directory/matrix.mtx
    # Remove all the other unnecessary files
  	rm $output_directory/no_header_bus_matrix.txt
  	rm $output_directory/gene_colmatrix.tsv
  	rm $output_directory/barcode_colmatrix.tsv
  	rm $output_directory/UMI_colmatrix.tsv
  	rm $output_directory/matrix_noheader.mtx
  else
    if [ $matrix_srt_name == 'False' ]
    then
      cp $matrix_file $output_directory/matrix.mtx.gz
    fi
  fi
fi
