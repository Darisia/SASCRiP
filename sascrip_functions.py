# This python script contains all the functions I've generated for single-cell RNA sequencing
# .. data preprocessing.

# ALL required modules
import subprocess
import os
from gen_func import mkdirs
# All things that need to be important for these the kb-python functions to work
import kb_python
from kb_python import ref
from kb_python import count
import gen_func # The python 'module' I created that contains general functions


#########################################################################################
# Import code from count.py from kb_python to define the bustools_text() function       #                                                                              #
from urllib.parse import urlparse                                                       #
                                                                                        #
from kb_python.config import get_bustools_binary_path, get_kallisto_binary_path         #
from kb_python.constants import (                                                       #
    ADATA_PREFIX,                                                                       #
    BUS_CDNA_PREFIX,                                                                    #
    BUS_FILENAME,                                                                       #
    BUS_FILTERED_FILENAME,                                                              #
    BUS_INTRON_PREFIX,                                                                  #
    BUS_S_FILENAME,                                                                     #
    BUS_SC_FILENAME,                                                                    #
    BUS_UNFILTERED_FILENAME,                                                            #
    COUNTS_PREFIX,                                                                      #
    ECMAP_FILENAME,                                                                     #
    FILTER_WHITELIST_FILENAME,                                                          #
    FILTERED_COUNTS_DIR,                                                                #
    INSPECT_FILENAME,                                                                   #
    TXNAMES_FILENAME,                                                                   #
    UNFILTERED_COUNTS_DIR,                                                              #
    WHITELIST_FILENAME,                                                                 #
)                                                                                       #
from kb_python.utils import (                                                           #
    copy_whitelist,                                                                     #
    stream_file,                                                                        #
    import_matrix_as_anndata,                                                           #
    overlay_anndatas,                                                                   #
    run_executable,                                                                     #
    whitelist_provided,                                                                 #
)                                                    #
#########################################################################################

# Set up the logging file
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s:%(levelname)s:\n%(message)s")
file_handler = logging.FileHandler("sascrip.log")
file_handler = setLevel(logging.DEBUG)
file_handler = setFormatter(formatter)
stream_handler = logging.StreamHandler()
logger.addHandler(file_handler)

############################################################################
# Function to install required R packages
############################################################################

def install_R_packages():

    '''

    Installs the R packages (tidyverse and Seurat) that are required to run SASCRiP functions

    '''

    command = 'Rscript ./install_packages.R'
    check_process = subprocess.run(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

    # Check to see if there was an error. If so, print the standard error - This should be edited to print readible errors
    if check_prcoess.returncode != 0:
        print(check_process.stderr.decode())

    # Add the standard output to the logfile
    logger.info(check_process.stdout.decode())


#####################################################################################################################################

############################################################################
# All functions for checking the quality of fastq sequencing (only)
############################################################################

### Define the function -- qc -- that will run the fastqc program
def qc(fastq_path, fastqc_output, fastqc_options_list):

    '''

    Runs fastqc on input fastq file through bash

    Parameters:
    fastq_path(str): the path for the fastq file
    fastqc_output(str): The path for the directory where the output files will be saved
    fastqc_options_list(list,str): A list containing the options you want to use in the fastqc program

    '''
    command = 'fastqc {} {} --outdir {}'.format(' '.join(fastqc_options_list), fastq_path, fastqc_output)
    subprocess.run(command, shell = True)
    check_fastqc = subprocess.run(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

    # Check to see if there was an error. If so, print the standard error - This should be edited to print readible errors
    if check_fastqc.returncode != 0:
        print(check_fastqc.stderr)

### Define the function -- loop -- that will loop the fastqc program
def run_fastqc(fastq_input_list, fastqc_output_list, fastqc_options_list):

    '''

    Loops the qc function so fastqc can be run on multiple fastq files from a directory containing the files

    Additionally - A list of folders containing the fastq files can be input here and qc will be run on all fastq files within those folders

    Parameters:
    fastq_input_list(list,str): A list containing the paths for the directories with the fastq.qz files
    fastqc_output_list(list,str): A list containing the paths for the directories you want the output fastqc reports to be saved
    fastqc_options_list(list,str): A list containing the options you want to use in the fastqc program

    '''
    for output_directory in fastqc_output_list:
        try:
            os.listdir(output_directory)
            break
        except FileNotFoundError:
            os.mkdir(output_directory)
    import fastqcfunc

    # Combine the input directory list and output directory list into a dictionary
    fastq_input_list.extend(fastqc_output_list)
    for i in range(0, len(fastqc_output_list), 1):
        x = i + len(fastqc_output_list)
        if x < len(fastq_input_list):
            list_dict = {fastq_input_list[i] : fastq_input_list[i + len(fastqc_output_list)] for i in range(0, len(fastqc_output_list), 1)}

    # Run the qc function on all files in all input directories
    for fastq_input, fastqc_output in list_dict.items():
        for file in os.listdir(fastq_input):
            if file.endswith('.fastq.gz'):
                fastqcfunc.qc(fastq_input + file, fastqc_output, fastqc_options_list)
    # End of function print statement
    print("Done")

#####################################################################################################################################

####################################################################################################
# Function to edit the 10xv1 fastq files into kallisto format for input into kallisto_bustools_count
####################################################################################################

def edit_10xv1_fastq(
	input_folder,
	output_folder):

	'''
	Separates the RA fastq 10xv1 file into two separate files where one file contains the UMI sequences and the other contains the transcript sequences

	Parameters:

	input_folder(str): Path to the directory containing the RA 10xv1 fastq files
	output_folder(str): Path to the output directory where the new separated fastq files will be saved

	'''

	# First we use the input folder - search for all files
	# .. that start with read-RA and put those filenames into a list
	# .. with the containing folder

	# Before the loop - let's create an empty list for all the RA files
	RA_10xv1_files = []

	for fastq_file in os.listdir(input_folder):
		if fastq_file.startswith("read-RA_"):
			RA_file_name = fastq_file
			RA_file_path = os.path.join(input_folder, RA_file_name)
			RA_10xv1_files.append(RA_file_path)

	# Once the list is created - need to combine the list into a string
	# .. separated by a comma (,)
	RA_10xv1_files_string = ",".join(RA_10xv1_files)

	# Edit folder paths
	if output_folder.endseith("/"):
		output_folder = os.path.dirname(output_folder)
	else:
		output_folder = output_folder

	# Let's create the bash command
	command = "./edit_fastq_function.sh {} {}".format(RA_10xv1_files_string, output_folder)

	# The Bash command needs to be run and checked
	subprocess.run(command, shell = True)

	check_process = subprocess.run(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

	if check_process.returncode != 0:
   		return(check_process.stderr.decode())

    logger.info(check_process.stdout.decode())

#####################################################################################################################################

####################################################################
# All functions for alignment and quantification - using kb-python
####################################################################

# Converts a bus file to a text file using the same layout as kb_python
def bustools_text(bus_path, out_path):
    logger.info(
        'Generating text file {} from BUS file {}'.format(out_path, bus_path)
    )
    command = [
        get_bustools_binary_path(), 'text', '-o', out_path, bus_path
    ]
    run_executable(command)
    return {'bus_text': out_path}


# Converts barcode, UMI, and transcript length into the bc:umi:seq kallisto bus format
def convert_bc_umi_seq(barcode_bp, umi_bp, transcript_bp):
    logger.info('Converting single-cell technology into bc:umi:seq format')
    barcode_bp = '14'
    tenxv1_tech_edit = '2,0,' + barcode_bp + ':1,0,' + umi_bp + ':0,0,' + transcript_bp
    return tenxv1_tech_edit

# Defines the function to filter a busfile using the barcodes and then generating the filtered count matrix
def filter_busfile(sorted_corrected_busfile, species_t2g, ecmap_path, transcripts_path, output, memory):

    logger.info('Filtering BUS file')

    # Set paths for all files and folders used and generated within this function
    generated_whitelist = os.path.join(output, 'generated_whitelist.txt')
    captured_busfile = os.path.join(output, 'captured_barcode.bus')
    filtered_busfile = os.path.join(output, 'filtered_barcode.bus')
    sorted_captured_busfile = os.path.join(output, 'sorted_captured.bus')
    path_to_prefix_count_files = 'filtered_counts'
    filtered_count_folder_files = os.path.join(output, path_to_prefix_count_files)

    # Make any directories if they don't exist already
    mkdirs(output)

    # Run all the functions from kb_python.count() that will filter and then count the busfile
    count.bustools_whitelist(sorted_corrected_busfile, generated_whitelist)
    count.bustools_capture(sorted_corrected_busfile, captured_busfile, generated_whitelist, ecmap_path, transcripts_path, capture_type='barcode')
    #count.bustools_sort(captured_busfile, sorted_captured_busfile, memory=memory)
    count.bustools_count(captured_busfile, filtered_count_folder_files, species_t2g, ecmap_path, transcripts_path)

    # Move files with a particular prefix to a specific folder
    move_with_prefix(output, 'filtered_counts.', filtered_count_folder_files)

# Defines the function to check if ERCCs are present within the dataset
def check_ercc(
            ERCC_fasta,
            all_out_path,
            list_of_fastqs,
            single_cell_technology,
            UMI_bp='0',
            barcode_bp='0',
            transcript_bp='0'
):
    # Make the output directory if it does not already exist
    mkdirs(all_out_path)

    ercc_log_file = os.path.join(all_out_path, 'check_ercc.log')

    # Create paths for all files and directories generated within this function
    all_ERCC_out_path = os.path.join(all_out_path, 'ERCC_analysis')
    ERCC_index_path = os.path.join(all_out_path, 'ERCC_index.idx')
    output_ERCC_busfile = os.path.join(all_ERCC_out_path, 'output.bus')
    ERCC_busfile_to_text = os.path.join(all_ERCC_out_path, 'ERCC_bus.txt')

    # Generate a kallisto index for the ERCC sequences using functions from the kb_python package
    ref.kallisto_index(ERCC_fasta, ERCC_index_path)

    # Generate BUSfiles - taking into consideration the technology
    if single_cell_technology=='10xv1':
        converted_tech = convert_bc_umi_seq(barcode_bp, UMI_bp, transcript_bp)
        count.kallisto_bus(list_of_fastqs, ERCC_index_path, converted_tech, all_ERCC_out_path, threads=8)
    else:
        count.kallisto_bus(list_of_fastqs, ERCC_index_path, single_cell_technology, all_ERCC_out_path, threads=8)

    # Convert the generated ERCC busfile to a human readable text file
    bustools_text(output_ERCC_busfile, ERCC_busfile_to_text)

    # Count the number of lines in the ERCC busfile
    busfile_lines = count_lines(ERCC_busfile_to_text)

    if busfile_lines >= 2:
        ERCC_included = True
    else:
        ERCC_included = False
    return ERCC_included

# Defines the function to generate unfiltered and filtered count matrices for the single-cell dataset of interest
def kallisto_bustools_count(
                        list_of_fastqs,
                        single_cell_technology,
                        all_out_path, # change to output_directory
                        generate_index = False,
                        species_index = None,
                        species_t2g = None,
                        species_fasta = None,
                        species_gtf = None,
                        k_mer_length = 31,
                        intron = False,
                        filter=True,
                        UMI_bp='0',
                        barcode_bp='0',
                        transcript_bp='0',
                        whitelist_path='None',
                        path_to_prefix_count_files='unfiltered_counts',
                        memory = '4G'
):
    # Make the output directories if it does not exist
    mkdirs(all_out_path)

    count_log_file = os.path.join(all_out_path, 'kb_count.log')

    # Generate the Kallisto index and the transcript_to_index mapping file if generate_index is True
    if generate_index==True:
        # Create paths for kallisto index and gene mapping file
        species_index = os.path.join(all_out_path, 'Kallisto_index.idx')
        species_t2g = os.path.join(all_out_path, 'transcripts_to_genes.txt')

        # Generate the kallisto_index
        if k_mer_length!=31:
            ref.kallisto_index(species_fasta, species_index, k = k_mer_length)
        else:
            ref.kallisto_index(species_fasta, species_index)

        # Generate the transcript_to_genes mapping file
        if intron!=False:
            ref.create_t2g_from_gtf(species_gtf, species_t2g, intron = intron)
        else:
            ref.create_t2g_from_gtf(species_gtf, species_t2g)

    # If no index/fasta is given - as default, use GRCh kallisto index that I created w hile ago
    if generate_index==False and species_index is None:
        # set species_index to the saved index and species_t2g to the saved gene mapping file
        species_index = './kallisto_index.idx'
        species_t2g = './transcript_to_genes.txt'

    # Create paths for all the other files and directories generated within this function
    all_count_out_path = os.path.join(all_out_path, 'Count_analysis')
    unfiltered_count_folder_files = os.path.join(all_count_out_path, path_to_prefix_count_files)
    output_busfile = os.path.join(all_count_out_path, 'output.bus')
    sorted_busfile = os.path.join(all_count_out_path, 'sorted.bus')
    inspect_json_file = os.path.join(all_count_out_path, 'inspect.json')
    ecmap = os.path.join(all_count_out_path, 'matrix.ec')
    corrected_busfile = os.path.join(all_count_out_path, 'corrected.bus')
    sorted_corrected_busfile = os.path.join(all_count_out_path, 'sorted_corrected.bus')
    transcripts = os.path.join(all_count_out_path, 'transcripts.txt')

    # Generate BUSfiles - taking into consideration the technology
    if single_cell_technology=='10xv1':
        converted_tech = convert_bc_umi_seq(barcode_bp, UMI_bp, transcript_bp)
        count.kallisto_bus(list_of_fastqs, species_index, converted_tech, all_count_out_path, threads=8)
    else:
        count.kallisto_bus(list_of_fastqs, species_index, single_cell_technology, all_count_out_path, threads=8)

    # Run bustools_sort and copy_or_create_whitelist from the kb-python package
    count.bustools_sort(output_busfile, sorted_busfile, threads=8, memory=memory)
    count.copy_or_create_whitelist(single_cell_technology, sorted_busfile, all_count_out_path)

    # Set the whitelist path for the different single-cell technologies
    if single_cell_technology=='10xv1' or \
    single_cell_technology=='10xv2' or \
    single_cell_technology=='10xv3':
        whitelist_path = all_count_out_path + '/' + single_cell_technology + '_whitelist.txt'
    elif single_cell_technology=='inDrops':
        whitelist_path = all_count_out_path + '/' + single_cell_technology + 'v3_whitelist.txt'
    else:
        whitelist_path = all_count_out_path + '/whitelist.txt'

    # Run functions from the kb-python package
    count.bustools_inspect(sorted_busfile, inspect_json_file, whitelist_path, ecmap)
    count.bustools_correct(sorted_busfile, corrected_busfile, whitelist_path)
    count.bustools_sort(corrected_busfile, sorted_corrected_busfile, threads=8, memory=memory)
    count.bustools_count(sorted_corrected_busfile, unfiltered_count_folder_files, species_t2g, ecmap, transcripts)

    # Move the generated unfiltered files into the unfiltered directory
    move_with_prefix(all_count_out_path, path_to_prefix_count_files+'.', unfiltered_count_folder_files)

    # Run filter_busfile if specified
    if filter==True:
        filter_busfile(sorted_corrected_busfile, species_t2g, ecmap, transcripts, all_count_out_path, memory=memory)

    # Return a dictionary containing generated files from this function
    gen_files = {
    "transcript_to_genes": species_t2g,
    "ec_mapping_file": ecmap
    }

    return(gen_files)

# Defines the function to check first if ERCCs are included within the dataset and proceeds accordingly
def include_ERCC_bus_count(
                        list_of_fastqs,
                        single_cell_technology,
                        all_out_path,
                        ERCC_fasta,
                        species_fasta,
                        generate_index = False,
                        species_index = None,
                        species_t2g = None,
                        species_gtf = None,
                        k_mer_length = 31,
                        intron = False,
                        filter = True,
                        UMI_bp='0',
                        barcode_bp='0',
                        transcript_bp='0',
                        whitelist_path='None',
                        path_to_prefix_count_files='unfiltered_counts',
                        memory = '4G'
):
    mkdirs(all_out_path)

    # Run the check_ercc function which returns a True or False
    is_ercc_included = check_ercc(
                            ERCC_fasta=ERCC_fasta,
                            all_out_path=all_out_path,
                            list_of_fastqs=list_of_fastqs,
                            single_cell_technology=single_cell_technology,
                            UMI_bp=UMI_bp,
                            barcode_bp=barcode_bp,
                            transcript_bp=transcript_bp
                        )
    # Do the following if ERCCs are included in the dataset
    if is_ercc_included == True:

        # Create paths for all the files and directories generated within this function
        all_ERCC_out_path = os.path.join(all_out_path, 'ERCC_analysis')
        species_ERCC_fasta = os.path.join(all_ERCC_out_path, 'species_ERCC_fasta.fa.gz')
        ERCC_combined_index = os.path.join(all_ERCC_out_path, 'ERCC_combined_index.idx')
        ERCC_annotations = os.path.join(all_ERCC_out_path, 'ERCC_annotations.txt')
        ERCC_transcripts = os.path.join(all_ERCC_out_path, 'transcripts.txt')
        ERCC_combined_t2g = os.path.join(all_ERCC_out_path, 'ERCC_combined_t2g.txt')

        # Combine the ERCC fasta file and the species fasta files
        concat_gz(species_ERCC_fasta, ERCC_fasta, species_fasta)

        # Generate the kallisto_index using functions from the kb_python package
        ref.kallisto_index(species_ERCC_fasta, ERCC_combined_index)

        # Generate a file in the same t2g format for the ERCC sequences
        with open(ERCC_annotations, 'w') as outfile:
            with open(ERCC_transcripts, 'r') as infile:
                infile_lines = infile.readlines()
                for line in infile_lines:
                    new_line = line.replace('\n', '\tERCC\tERCC\n')
                    outfile.write(new_line)

        # Combine the ERCC t2g and species t2g
        concat_files(ERCC_combined_t2g, ERCC_annotations, species_t2g)

        species_index = ERCC_combined_index
        species_t2g = ERCC_combined_t2g

    else:
        species_index = species_index
        species_t2g = species_t2g

    # Run the kallisto_bustools_count function to generate unfiltered and filtered count matrices for the dataset of interest
    kallisto_bustools_count(
                        list_of_fastqs = list_of_fastqs,
                        single_cell_technology = single_cell_technology,
                        all_out_path = all_out_path,
                        generate_index = generate_index,
                        species_index=species_index,
                        species_t2g=species_t2g,
                        species_fasta = species_fasta,
                        species_gtf = species_gtf,
                        k_mer_length = k_mer_length,
                        intron = intron,
                        filter = filter,
                        UMI_bp=UMI_bp,
                        barcode_bp=barcode_bp,
                        transcript_bp=transcript_bp,
                        whitelist_path=whitelist_path,
                        path_to_prefix_count_files=path_to_prefix_count_files,
                        memory = memory
    )

####################################################################################################################

#######################################################################################################
# Function for editing BUStools output for input into Seurat (cell_cqc)
#######################################################################################################
def seurat_matrix(
    bustools_mtx_matrix,
    bustools_gene_index,
    bustools_barcode_index,
    output_directory
):

    '''
	Runs a bash function that edits the bustools feature file so that is resembles a CellRangers features file so it can be input into Seurat

	Parameters:
    bustools_mtx_matrix(str): Path to the bustools mtx matrix file
	bustools_gene_index(str): Path to the output bustools gene index file
    bustools_barcode_index(str): Path to the bustools barcodes file
	transcripts_to_genes_file(str): Path to the transcripts to genes file
	output_directory(str): Path to the output directory where the matrix files will be saved

	'''

    logger.info("Rearranging bustools output matrix for input into Seurat")

    # Let's start by checking and/or creating the output directory
	# .. That means we need the gen_func script
	if output_directory == "working_directory":
		output_directory = "./"
	else:
		gen_func.mkdirs(output_directory)

    ## First edit the gene index
    # Let's put together the bash command
	command = "sh ./bus2CR_features.sh {} {}".format(
		bustools_gene_index,
		output_directory)

	# Run the bash command using subproces.run and save the standard output to check_process
	check_process = subprocess.run(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

	# Check the returncode - if there were errors - print the standard errors
	if check_process.returncode != 0:
		print(check_process.stderr.decode())

    logger.info(check_process.stdout.decode())

    # Edit the barcode file
    # put the bash command together
	command = "sh ./bus2CR_barcodes.sh {} {}".format(
		bustools_barcode_index,
		output_directory)

	# Run the bash command using subprocess.run and save the standard output to check_process
	check_process = subprocess.run(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

	# check the returncode - if there are errors - print the standard errors
	if check_process.returncode != 0:
		print(check_process.stderr.decode())

    logger.info(check_process.stdout.decode())

    # Edit the mtx output matrix
    # put the bash command together
	command = "sh ./bus2CR_matrix.sh {} {}".format(
		bustools_mtx_matrix,
		output_directory)

	# Run the bash command using subprocess.run and save the standard output to check_process
	check_process = subprocess.run(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

	# check the returncode - if there are errors - print the standard errors
	if check_process.returncode != 0:
		print(check_process.stderr.decode())

    logger.info(check_process.stdout.decode())


####################################################################################################################

#######################################################################################################
# All functions for CellQC preprocessing steps and the generating of graphs - still needs more editing
#######################################################################################################

## Create the python function
def run_cqc(input_file_or_folder,
            sample_ID,
            output_folder = "working_directory",
            generate_seurat_object = True,
            subset_seurat_object = True,
            generate_default_plots = True,
            input_seurat_object = False,
            transcripts_to_genes_file = None,
            gene_lower = 200,
            gene_higher_method = "MAD",
            gene_higher = "to_be_calculated",
            mitochondria_percent = 10,
            nMADs = 6,
            nSD = 6,
            extract_cell_metrics = False,
            output_matrix = False
):

    '''
    Performs cell quality control using Seurats quality control functions (wrapper function)

    Parameters:
    input_file_or_folder(str): Path to the folder containing the output matrix or path to the hdf5 file
    sample_ID(str): The name of the sample
    output_folder(str): The path to the output directory where all output will be saved (default = working directory)
    generate_seurat_object(bool): State whether to generate the Seurat object from output matrix (default = True)
    subset_seurat_object(bool): State whether to subset the Seurat object, removing cells identified as damaged/Multiplets (default = True)
    generate_default_plots(bool): State whether to generate all default visualisations (default = True)
    input_seurat_object(bool): State whether the input file is a saved seurat object (default = False)
    ENSG_gname38_path(str): Path to the file called ENSG_gname38.tsv that will allow us to convert ENSG name into hgnc symbols
    gene_lower(int/None): Minimum number of genes required (default = 200)
    gene_higher_method(str): One of three methods - MAD, SD, or Manual to identify outliers within the dataset
    gene_higher(int/None): If method selected as Manual - the maximum number of genes required should be given here
    mitochondria_percent(int/None): The maximum percentage of mitochondrial genes within a cell (default = 10)
    nMADs(int): If using MAD method state the threshold calculation to determine outliers (default = 6)
    nSD(int): If using SD method state the threshold calculation to determine outliers (default = 6)
    extract_cell_metrics(bool): State whether the cell metrics should be extracted and saved in a tsv file (default = False)
    output_matrix(bool): State whether to generate an mtx matrix file from the Seurat filtered data (default = False)

    '''
    # Check if output folder is specified
    if output_folder == "working_directory":
        output_folder = "./"
    else:
        output_folder = output_folder

    # Edit directory paths
    if output_folder.endswith('/'):
        output_folder = os.path.dirname(output_folder)
    else:
        output_folder = output_folder

    # Generate the directories if it does not already exist
    if output_folder != "working_directory":
        mkdirs(output_folder)

    # Convert the None inputs to a string to be input into bash
    if gene_lower is None:
        str(gene_lower)
    if mitochondria_percent is None:
        str(mitochondria_percent)
    if gene_higher is None:
        str(gene_higher)

    # Generate the directory where the output mtx matrix would be saved
    if output_matrix == True:
        output_matrix_dir = os.path.join(output_folder, sample_ID + '_QC_matrix_output')
        mkdirs(output_matrix_dir)



    # Set Seurat object if the seurat object was given instead of the output matrix
    if input_seurat_object is True:
        seurat_object = input_file_or_folder

    # State logging statements for the thresholds
    if generate_default_plots is True:
        if (gene_lower == "None" or mitochondria_percent == "None"):
            logger.error("Lower gene cutoff or mitochondria percent cutoff missing: Default plots showing damaged cells cannot be generated")
        if (gene_higher == "None"):
            logger.error("Upper gene cutoff missing: Default plots showing cell multiplets cannot be generated")

    # Convert pythons False to 'FALSE' and True to 'TRUE'
    if generate_seurat_object is False:
        generate_seurat_object = 'FALSE'
    elif generate_seurat_object is True:
        generate_seurat_object = 'TRUE'

    if subset_seurat_object is False:
        subset_seurat_object = 'FALSE'
    elif subset_seurat_object is True:
        subset_seurat_object = 'TRUE'

    if generate_default_plots is False:
        generate_default_plots = 'FALSE'
    elif generate_default_plots is True:
        generate_default_plots = 'TRUE'

    if extract_cell_metrics is False:
        extract_cell_metrics = 'FALSE'
    elif extract_cell_metrics is True:
        extract_cell_metrics = 'TRUE'

    if output_matrix is False:
        output_matrix = 'FALSE'
    elif output_matrix is True:
        output_matrix = 'TRUE'

    # if transcripts_to_genes supplied - create ensg_gname
    if transcripts_to_genes_file != None:
        command = 'sh cut_t2g_sh {} {}'.format(
        transcripts_to_genes_file,
        output_folder
        )

        check_process = subprocess.run(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

        # check to see if that variable has a returncode above 0 - if it does, print the error message
        if check_command.returncode != 0:
            print(check_command.stderr.decode())
            print(check_command.stdout.decode())

        # Print all standard output to log file
        logger.info(check_command.stdout.decode())

        ensg_gname_path = os.path.join(output_folder, "ensg_gname.tsv")
    else:
        logger.warning("No transcripts_to_genes_file given. If ESEMBL gene names are used, a Seurat object cannot be properly generated and an error will be returned")
        ensg_gname_path = "None"


    # The bash commands to run the cell quality control R script
    command = 'Rscript cqc_srt.R {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
    input_file_or_folder,
    sample_ID,
    output_folder,
    generate_seurat_object,
    subset_seurat_object,
    generate_default_plots,
    input_seurat_object,
    ensg_gname_path,
    gene_lower,
    gene_higher_method,
    gene_higher,
    mitochondria_percent,
    nMADs,
    nSD,
    extract_cell_metrics,
    output_matrix
    )

    # Python function (subprocess) to run the bash commands and store standard output in check_command
    logger.info('Running cqc_srt.R through the shell')
    check_command = subprocess.run(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

    # check to see if that variable has a returncode above 0 - if it does, print the error message
    if check_command.returncode != 0:
        print(check_command.stderr.decode())
        print(check_command.stdout.decode())

    # Print all standard output to log file
    logger.info(check_command.stdout.decode())
    logger.info("run_cqc has finished running")

#########################################################################################################

##################################################################################################
# All functions to run the normalisation preprocessing steps on the filtered single-cell UMI data
##################################################################################################

# what parameters will I need for the function?
# 1) Path to the seurat object
# 2) the sample_ID - so I know what to save the output as
# 3) Path to the output directory - optional - otherwise created in the working directory
# 4) **kwargs - additional arguments to be passed to the actual SCTransform function
# 5) ENSG_gname38_path - the path to the file containing the ENSG names and gene symbols
# 6) cell_cycle_normalisation - boolean - determine if cell-cycle-normalisation should be performed aswell
# 7) considered parameter: separate_files - boolean - choose one output or two outputs if cell_cycle_normalisation is True ##

# sctransform function that Runs SCTransform from within Seurat to normalise the single-cell data for seuqencing depth
def sctransform_normalize(
	seurat_object,
	sample_ID,
	output_directory = "working_directory",
    ouput_log_matrix = False,
    output_count_matrix = False,
	ENSG_gname38_path = "working_directory",
	**additional_sctransform_arguments): # Don't forget that there were some issues with this

	'''
	Wraps the SCTransform function within Seurat to perform normalisation on the raw counts

	Parameters:
	seurat_object(str): Path to the saved seurat object
	sample_ID(str): Name of the sample
	output_directory(str): Path to the output directory where all generated files and folders will be saved
	ENSG_gname38_path(str): Path to the file called ENSG_gname38.tsv that will allow us to convert ENSG name into hgnc symbols
	**additional_sctransform_arguments: additional arguements (with key words) will be passed to Seurat's SCTransform function

	'''

	# Sort out the output directory
	if output_directory == "working_directory":
		output_directory = "./"
	else:
		gen_func.mkdirs(output_directory)

	# Sort out the ENSG_gname38_path path
	if ENSG_gname38_path == "working_directory":
		ENSG_gname38_path = "./"
	else:
		ENSG_gname38_path = ENSG_gname38_path

	# Create list of values from the additional argument dictionary
	# .. if there is additional arguments - should test first
	if len(additional_sctransform_arguments) == 0:
		is_add_args = False
	else:
		is_add_args = True

	if is_add_args is True:

		# put the values from additional argument dictionary into a list
		additional_argument_values = list(additional_sctransform_arguments.values())

		# Use that list in the get_data_type function from the gen_func "module"
		data_type_information_both = gen_func.get_data_type(
			additional_argument_values,
			return_type = "both",
			sep_string_by = ";")

		# Extract the data type information string
		data_type_information = data_type_information_both[1]

		# Extract the data type information list
		data_type_information_list = data_type_information_both[0]

		# Check if any input parameters is classified as an unknown data type
		if "unknown" in data_type_information_list:
			data_type_error = True
		else:
			data_type_error = False

		# Convert the additional arguments into a long string using function from gen_func "module"
		add_args_list = gen_func.argument_to_string(
			additional_sctransform_arguments,
			return_type = "string",
			sep_string_by = ";")

		# All additional argument information is given

	else:

		data_type_information = ""
		add_args_list = ""

	# all variables for the bash code
	# 1) path to seurat object
	# 2) sample_ID
	# 3) output_directory
	# 4) ENSG_gname_path
	# 5) cell_cycle_normalisation
	# 6) is_add_args
	# 7) data_type_information
	# 8) add_args_list

	# let's put together the bash code
	command = "Rscript ./normalise_seurat.R {} {} {} {} {} {} {}".format(
		seurat_object,
		sample_ID,
		output_directory,
		ENSG_gname38_path,
		is_add_args,
		data_type_information,
		add_args_list)

	# Use subprocess to run the bash command through python
	check_process = subprocess.run(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

	# Check to see if there was an error. If so, print the standard error
	if check_process.returncode != 0:
		print(check_process.stderr)

    # Print all standard output to the log file
    logger.info(check_process.stdout.decode())
    logger.info("sctransform_normalize has finished running")
