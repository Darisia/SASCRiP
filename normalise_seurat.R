# This R script contains the code run sctransform
# which will be run through python

# required packages
library(tidyverse)
library(Seurat)

# Pre-defined R functions

# Convert the values in a list to their correct data type
# .. given the data type informations as a character vector
convert_data_type <- function(data_type_info_vect, list_of_values) {

	# Convert all argument values from the list_of_values list into their correct data type
	for ( number in seq(1, length(list_of_values)) ) {
		# string to correct data type
		if ( data_type_info_vect[number] == "num" ) {
			list_of_values[number] <- as.double(list_of_values[number])

		} else if ( data_type_info_vect[number] == "num_vect" ) {
			list_of_values[number] <- list(c(as.double(unlist(strsplit(unlist(list_of_values[number],
				use.names = FALSE),
			split = ",")))))

		} else if ( data_type_info_vect[number] == "char_vect" ) {
			list_of_values[number] <- list(c(unlist(strsplit(unlist(list_of_values[number],
				use.names = FALSE),
			split = ","))))

		} else if ( data_type_info_vect[number] == "logical" ) {
			list_of_values[number] <- as.logical(list_of_values[number])

		} else if ( data_type_info_vect[number] == "log_vect" ) {
			list_of_values[number] <- list(c(as.logical(unlist(strsplit(unlist(list_of_values[number],
				use.names = FALSE),
			split = ",")))))

		} else if ( data_type_info_vect[number] == "expression" ) {
			list_of_values[number] <- eval(parse(text = list_of_values[number]))

		} else if ( data_type_info_vect[number] == "expr_vect" ) {
			expr_vect <- c()
			for ( expression in parse(text = unlist(strsplit(unlist(list_of_values[number],
				use.names = FALSE),
			split = ","))) ) {
				expr_vect <- c(expr_vect, eval(expression))
			}
			list_of_values[number] <- list(c(expr_vect))

		} else if ( data_type_info_vect[number] == "list" ) {
			list_of_values[number] <- list(c(as.list(unlist(strsplit(unlist(list_of_values[number],
				use.names = FALSE),
			split = ",")))))

		} else if ( data_type_info_vect[number] == "NULL" ) {
			list_of_values[number] <- list(as.null(list_of_values[number]))
		}

	}
  return(list_of_values)
}

# General R function to load an RData file using a specified object name
loadRData <- function(filename) {
  load(filename)
  get(ls()[ls() != "filename"])
}

# R function to replace the ENSEMBL IDs with hgnc gene names within a Seurat object
ensg_to_hgnc <- function(
  ENSG_gname38_path = "working_directory",
  seurat_object,
  sample_ID
){
  # Set the path to the ENSg_gname38.tsv file
  if (ENSG_gname38_path == "working_directory") {
    ENSG_gname38_path = "./ENSG_gname38.tsv"
  } else {
    ENSG_gname38_path = ENSG_gname38_path
  }
  # read in the ENSG_gname38.tsv file
  ENSG_gname38 <- read_tsv(ENSG_gname38_path)
  # Read in the saved seurat object
  if (endsWith(seurat_object, ".rds") == TRUE | endsWith(seurat_object, "RDS") == TRUE) {
    sample_ID_srt <- readRDS(file = seurat_object)
  } else {
    sample_ID_srt <- loadRData(seurat_object)
  }
  # Update Seurat object if saved as an old seurat object
  sample_ID_srt <- UpdateSeuratObject(sample_ID_srt)
  # Generate the dataframe that will be used to match and replace the ENSEMBL IDs
  pmatch_output <- data.frame(pmatch(ENSG_gname38$ENSEMBL_ID, rownames(sample_ID_srt@assays[["RNA"]]@counts)))
  colnames(pmatch_output) <- "gene_match"
  pmatch_output$match_index <- rownames(pmatch_output)
  pmatch_output <- pmatch_output[order(pmatch_output$gene_match), ]
  pmatch_output <- filter(pmatch_output, !is.na(gene_match))
  # Match and replace the ENSEMBL IDs in the relevant slots in the seurat object with the hgnc gene symbol
  rownames(sample_ID_srt@assays[["RNA"]]@counts) <- ENSG_gname38$gene_name[as.double(pmatch_output$match_index)]
  # Return the Seurat object with hgnc symbols as rownames
  return(sample_ID_srt)
}

# Read in and check the seurat object
read_check_srt <- function(
  seurat_object,
  sample_ID,
  ENSG_gname38_path = "working_directory"
){
  # Read in the saved seurat object
  if (endsWith(seurat_object, ".rds") == TRUE | endsWith(seurat_object, "RDS") == TRUE) {
    sample_ID_srt <- readRDS(file = seurat_object)
  } else {
    sample_ID_srt <- loadRData(seurat_object)
  }
  # Update Seurat object if saved as an old seurat object
  sample_ID_srt <- UpdateSeuratObject(sample_ID_srt)
  # Check if rownames use hgnc symbols or ENSEMBL IDs - change to hgnc symbol if ENSEMBL IDs are used
  if (startsWith(rownames(sample_ID_srt@assays[["RNA"]]@counts)[1], "ENSG") == TRUE) {
    sample_ID_srt <- ensg_to_hgnc(ENSG_gname38_path, seurat_object, sample_ID)
  }
  # Remove cell barcodes with counts of 0
  sample_ID_srt <- subset(sample_ID_srt, subset = nCount_RNA > 0)
  # Check if the mito.percent column is in the meta.data slot of the Seurat object
  if ("mito.percent" %in% colnames(sample_ID_srt@meta.data)) {
    sample_ID_srt <- sample_ID_srt
  } else {
    sample_ID_srt[['mito.percent']] <- PercentageFeatureSet(sample_ID_srt, pattern = '^MT-')
  }
  return(sample_ID_srt)
}



# All arguments input in Bash to be saved to the args list
args <- commandArgs(TRUE)

# State all variables here:
# 1) path to the saved seurat object
# 2) sample_ID
# 3) ooutput_directory
# 4) ENSG_gname38_path
# 5) cell_cycle_normalisation
# 6) is_add_args = as.logical()
# 7) add_args_list
# 8) data_type_information

# Define all variables from bash arguments
seurat_object <- args[1]
sample_ID <- args[2]
output_directory <- args[3]
ENSG_gname38_path <- args[4]
is_add_args <- as.logical(args[5])
data_type_information <- args[6]
add_args_list <- args[7]

# Read in the filtered Seurat Object
sample_ID_srt <- read_check_srt(seurat_object, sample_ID)

# Check if there are any additional arguments to process
if ( is_add_args == TRUE ) {

	# Convert data type information into a character vector
	data_type_information <- unlist(strsplit(data_type_information, ";"))

	# Convert the additional argument string into a list
	# .. where each argument is a different item

	# Split the string by ; and create a character vector
	add_args_split <- unlist(strsplit(add_args_list, ";"))

	# Split each character by the first equal sign separating the argument name and the argument value
	add_args_name <- regmatches(
		add_args_split,
		regexpr("=", add_args_split),
		invert = TRUE) # Should produce the list where different arguments are separate elements within the list

	# Create new list with only the argvals as the elements within the list
	additional_argument_values <- list()
	additional_argument_names <- c()

	for ( argument in add_args_name ) {
		argument_name <- argument[1]
		argument_value <- argument[2]

		additional_argument_names <- c(additional_argument_names, argument_name) # character vector
		additional_argument_values <- c(additional_argument_values, argument_value) # list of characters (for now)

	}

	# Set the argument names as the names of each of the elements (argvals) in the additional_argument_values list
	names(additional_argument_values) <- additional_argument_names

	# Convert the values in the additional_argument_values into their
	# .. correct data type
	additional_argument_values <- convert_data_type(
		data_type_information,
		additional_argument_values)

	# Convert that list into a string for eval_parse and then assign the argname to the argval
	vector_with_args <- c()

	for ( argument_pos in additional_argument_values ) {
		argname <- names(additional_argument_values)[argument_pos]
		argval <- additional_argument_values[[argument_pos]]
		argument_string <- sprintf("%s = %s", argname, argname)
		vector_with_args <- c(vector_with_args, argument_string)
		assign(argname, argval)
	}

	string_with_args <- paste(vector_with_args, collapse = ", ")
	function_with_args <- paste0("SCTransform(sample_ID_srt, ", string_with_args, ")")
	sample_ID_srt <- eval(parse(text = function_with_args))
} else {
	sample_ID_srt <- SCTransform(sample_ID_srt)
}

# Save the normalised Seurat object to rds file
saveRDS(
  sample_ID_srt,
  file = sprintf("%s/%s_normalised_seurat.rds", output_directory, sample_ID)
)
