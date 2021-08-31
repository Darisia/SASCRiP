# SASCRiP Documentation
This repository contains all the code and documentation for the SASCRiP Python package

## Table of Contents
1. [Installation](https://github.com/Darisia/SASCRiP/blob/main/README.md#installation)
2. <details>
     <summary><a href="https://github.com/Darisia/SASCRiP/blob/main/README.md#sascrip-functions:-user-guide">SASCRiP functions: User guide</a></summary>
     <br>
     <a href="https://github.com/Darisia/SASCRiP/blob/main/README.md#edit_10xv1_fastq">edit_10v1_fastq</a>
     <br>
     <a href="https://github.com/Darisia/SASCRiP/blob/main/README.md#check_ercc">check_ercc</a>
     <br>
     <a href="https://github.com/Darisia/SASCRiP/blob/main/README.md#kallisto_bustools_count">kallisto_bustools_count</a>
     <br>
     <a href="https://github.com/Darisia/SASCRiP/blob/main/README.md#include_ercc_bus_count">include_ERCC_bus_count</a>
     <br>
     <a href="https://github.com/Darisia/SASCRiP/blob/main/README.md#run_cqc">run_cqc</a>
     <br>
     <a href="https://github.com/Darisia/SASCRiP/blob/main/README.md#stransform_normalize">sctransform_normalize</a>
     <br>
     <a href="https://github.com/Darisia/SASCRiP/blob/main/README.md#sascrip_preprocess">sascrip_preprocess</a>
     </details>
3. [Authors](https://github.com/Darisia/SASCRiP#authors)
4. [License](https://github.com/Darisia/SASCRiP#license)
5. [Reference](https://github.com/Darisia/SASCRiP#reference)
6. [Read next](https://github.com/Darisia/SASCRiP#read-next)


## Installation

### Requirements

SASCRiP uses multiple single-cell analysis packages such as Seurat and kb-python. Since SASCRiP makes use of the R packages such as Seurat and Tidyverse for plotting, these packages are required. A full list of the requirements is shown below

1. Python (>v3.7) is required to run SASCRiP functions
2. R (>v3.6) is required to be installed 
3. Seurat R package (can be installed through SASCRiP)
4. Tidyverse R package (can be installed through SASCRiP)
5. kb-python (can be installed when SASCRiP is installed)

### Installation code

The SASCRiP package can be installed using pip from the terminal

```bash

pip install sascrip

```

## SASCRiP functions: User guide

## `install_R_packages`

This function allows the user to install missing R packages that are required for SASCRiP to work. `install_R_packages` first checks if the package is installed and if not, install.packages() is run.

#### Usage

````python
import sascrip
from sascrip import sascrip_functions

sascrip_functions.install_R_packages()

````  

## `edit_10xv1_fastq`

`edit_10xv1_fastq` prepares FastQ files, obtained using the 10xv1 sequencing chemistry, for input into Kallisto for pseudoalignment. The directory/s containing the FastQ files are input into `edit_10xv1_fastq` which searches for the RA FastQ file that contains both the UMI and transcript sequence. The UMI and Transcript sequences are separated into their own FastQ files. The directory/s containing all these FastQ files can be input into kallisto_bustools_count for further processing. The produced FastQ files will be saved in the same directory containing all the input FastQ files.

#### Usage

````python
import sascrip
from sascrip import sascrip_functions

sascrip_functions.edit_10xv1_fastq(input_directories)

````

#### Required parameters

| Parameter | Description |
| --- | --- |
| `input_directories` (str-list) | Path to the directory/s containing the RA 10xv1 FastQ files. If there are more than one directories. A list should be given with all input directories |


## `check_ercc`

`check_ercc` allows the user to check if the single-cell dataset may contain RNA spike-ins that can be used as an additional cell-quality control metric. If `check_ercc is run`, True or False - depending on whether ERCCs are included or not, will be returned as standard output 

#### Usage

````python
import sascrip
from sascrip import sascrip_functions

sascrip_functions.check_ercc(
     ERCC_fasta,
     output_directory_path,
     list_of_fastqs,
     single_cell_technology,
     input_directory = False,
     read_separator = None,
     UMI_bp = '0',
     barcode_bp = '0',
     transcript_bp = '0'
)

````
#### Required Parameters
___________________


| Parameter | Description |
| --- | --- |
| `ERCC_fasta` (str) | Path to the ERCC FASTA file |
| `output_directory_path` (str) | Path to the output directory where output files will be saved |
| `list_of_fastqs` (str-list) | Python list of the paths to input FastQ files in the *order specified by Kallisto. The folder containing all the input files can be given instead and the relevant FastQ files will be captured and sorted in the correct order. To use this feature - input the path to the directory here (If multiple directories are used, input all the directories as a list), set input_directory = True and provide the strings used to separate the reads in the read_separator parameter. * See example below |
| `single_cell_technology` (str) | The single-cell sequencing technology that was used as *specified by Kallisto. If 10xv1 technology was used, the UMI_bp and barcode_bp parameters are required |

#### Optional parameters
___________________

| Parameter | Description |
| --- | --- |
|`input_directory` (bool) | Indicate whether the list_of_fastqs parameter is given the path to the directory containing all input FastQ files * See example below |
| `read_separator` (str-list) | The strings used to separate the input reads. Within a list - the first element should contain the string used to identify read 1, the second element should contain the string used to identify read 2 and if "10xv1" FastQ files are used, a third element is required that contains the string used to identify read 3. * See example below |
| `UMI_bp` (str) | The number of base pairs sequenced for the UMI sequence. If 10xv1 technology is used, this parameter is required |
| `barcode_bp` (str) | The number of base pairs sequenced for the barcode sequence. If 10xv1 technology is used, this parameter is required |
| `transcript_bp` (str) | The number of base pairs sequenced for the transcript sequence |


* **Kallisto specified FastQ file order and single-cell technologies**

| Single-cell tech |     FastQ file order     |
| ---------------- | ------------------------ |
| 10xv1            | Transcript, UMI, barcode |
| 10xv2            | barcode-UMI, Transcript  |
| 10xv3            | barcode-UMI, Transcript  |
| CELSeq           | barcode-UMI, Transcript  |
| CELSeq2          | barcode-UMI, Transcript  |
| DropSeq          | barcode-UMI, Transcript  |
| inDrops          | barcode-UMI, Transcript  |
| SCRBSeq          | barcode-UMI, Transcript  |
| SureCell         | barcode-UMI, Transcript  |

* **Working with more than one set of fastq files**

Include the additional sets in the fastq list, keeping the sets together. For example: if you have 2 sets of FastQ files from the 10xv2 chemistry

```python

list_of_fastqs = ["barcode_UMI_R1_1.fastq.gz", "Transcript_R2_1.fastq.gz", "barcode_UMI_R1_2.fastq.gz", "Transcript_R2_2.fastq.gz"]

```
Additionally, the directory/s containing all the relevant FastQ files can be given to the `list_of_fastqs` parameter and the FastQ files will be captured and sorted within the function, provided that `input_directory` is set to True and the `read_separator` is given. An example showing how these parameters should be set is shown below. 

Two sets of FastQ files (as shown below) are contained within a directory called "FastQ_directory":
barcode_UMI_R1_1.fastq.gz
Transcript_R2_1.fastq.gz
barcode_UMI_R1_2.fastq.gz
Transcript_R2_2.fastq.gz

```python

list_of_fastqs = "path/to/FastQ_directory"
input_directory = True
read_separator = ["R1", "R2"]

```


## `kallisto_bustools_count`

Runs both Kallisto and Bustools through kb-python to perform pseudoalignment and gene quantification respectively. `kallisto_bustools_count` generates unfiltered and filtered count matrices for the single-cell dataset of interest

#### Usage

````python
import sascrip
from sascrip import sascrip_functions

sascrip_functions.kallisto_bustools_count(
     list_of_fastqs,
     single_cell_technology,
     output_directory_path,
     species_index,
     species_t2g,
     input_directory = False,
     read_separator = None,
     generate_index = False,
     species_fasta = None,
     species_gtf = None,
     k_mer_length = 31,
     intron = False,
     filter = True,
     UMI_bp = '0',
     barcode_bp = '0',
     transcript_bp = '0',
     whitelist_path = None,
     path_to_prefix_count_files = 'unfiltered_counts',
     memory = '4G'
)

````
#### Required parameters
___________________

| Parameter | Description |
| --- | --- |
| `list_of_fastqs` (str-list) | Python list of the paths to input FastQ files in the *order specified by Kallisto. The folder containing all the input files can be given instead and the relevant FastQ files will be captured and sorted in the correct order. To use this feature - input the path to the directory here (If multiple directories are used, input all the directories as a list), set input_directory = True and provide the strings used to separate the reads in the read_separator parameter. * See example above (working with more than one set of FastQ files) |
| `single_cell_technology` (str) | The single-cell sequencing technology that was used as *specified by Kallisto. If 10xv1 technology was used, the UMI_bp and barcode_bp parameters are required |
| `output_directory_path` (str) | Path to the main output directory where all files and new directories will be created and saved |
| `species_index` (str) | Path to the kallisto_index for the species of interest. If the Kallisto index needs to be generated, set species_index = None, generate_index = True and species_fasta = "path/to/cDNA FASTA". |
| `species_t2g` (str) | Path to the transcript-to-genes mapping file for the species of interest. If the transcripts-to-genes mapping file needs to be generated, set species_t2g = None, generate_index = True and species_gtf = "path/to/GTF file". |


#### Optional parameters
___________________

| Parameter | Description |
| --- | --- |
| `input_directory` (bool) | Indicate whether the list_of_fastqs parameter is given the path to the directory containing all input FastQ files * See example above (working with more than one set of FastQ files) |
| `read_separator` (str-list) | The strings used to separate the input reads. Within a list - the first element should contain the string used to identify read 1, the second element should contain the string used to identify read 2 and if "10xv1" FastQ files are used, a third element is required that contains the string used to identify read 3. * See example above (working with more than one set of FastQ files) |
| `generate_index` (bool) | Indicate whether the Kallisto Index should be generated within the function. If set to True: the species_fasta and the species_gtf are required |
| `species_fasta` (str) | Path to the transcriptome (cDNA) FASTA file for the species of interest. This will be used to generate the Kallisto index. If generate_index is set to True, this parameter is required |
| `species_gtf` (str) | Path to the GTF file for the species of interest. This will be used to create the transcripts-to-genes mapping file. If generate_index is set to True, this parameter is required. |
| `k_mer_length` (int) | The length of the K-mers that should be generated when creating the Kallisto index |
| `intron` (bool) | Indicate whether or not to include intron transcript ids |
| `filter` (bool) | Indicate whether or not to filter the BUS file prior to generating the gene-count matrix in mtx format |
| `UMI_bp` (str) | The number of base pairs sequenced for the UMI sequence. If 10xv1 technology is used, this parameter is required |
| `barcode_bp` (str) | The number of base pairs sequenced for the barcode sequence. If 10xv1 technology is used, this parameter is required |
| `transcript_bp` (str) | The number of base pairs sequenced for the transcript sequence |
| `whitelist_path` (str) | Path to the barcode whitelist that will be used for barcode correction |
| `path_to_prefix_count_files` (str) | Prefix of the output matrix files and indices |
| `memory` (str) | Amount of memory to use |


* **Kallisto specified FastQ file order and single-cell technologies**

Please see table above

* **Working with more than one set of fastq files**

Please see description above


## `include_ERCC_bus_count`

`include_ERCC_bus_count` first checks if ERCC spike-ins are included in the dataset and generates counts with `kallisto_bustools_count` accordingly. To check if spike-ins are present, a Kallisto index is generated from the ERCC FASTA and the FastQ sequences from the dataset are aligned to the ERCC index. if no sequences align to the ERCC index, the gene-count matrix will be generated without including ERCC sequences. However, if there are seqeunces that align to the ERCC index, the ERCC FASTA seqeunce and the species FASTA sequence will be combined to create a new kallisto index. Therefore, the final gene-count matrix will include ERCC sequences.

#### Usage

````python
import sascrip
from sascrip import sascrip_functions

sascrip_functions.include_ERCC_bus_count(
     list_of_fastqs,
     single_cell_technology,
     output_directory_path,
     ERCC_fasta,
     species_index,
     species_t2g,
     species_fasta,
     input_directory = False,
     read_separator = None,
     generate_index = False,
     species_gtf = None,
     k_mer_length = 31,
     intron = False,
     filter = True,
     UMI_bp = '0',
     barcode_bp = '0',
     transcript_bp = '0',
     whitelist_path = None,
     path_to_prefix_count_files = 'unfiltered_counts',
     memory = '4G'
)

````
#### Required parameters
___________________

| Parameter | Description |
| --- | --- |
|`list_of_fastqs` (str-list) | Python list of the paths to input FastQ files in the *order specified by Kallisto. The folder containing all the input files can be given instead and the relevant FastQ files will be captured and sorted in the correct order. To use this feature - input the path to the directory here (If multiple directories are used, input all the directories as a list), set input_directory = True and provide the strings used to separate the reads in the read_separator parameter. * See example above (working with more than one set of FastQ files) |
| `single_cell_technology` (str) | The single-cell sequencing technology that was used as *specified by Kallisto. If 10xv1 technology was used, the UMI_bp and barcode_bp parameters are required |
| `output_directory_path` (str) | Path to the main output directory where all files and new directories will be created and saved |
| `ERCC_fasta` (str) | Path to the ERCC FASTA file
| `species_index` (str) | Path to the kallisto_index for the species of interest. If the Kallisto index needs to be generated, set species_index = None, generate_index = True and species_fasta = "path/to/cDNA FASTA". |
| `species_t2g` (str) | Path to the transcript-to-genes mapping file for the species of interest. If the transcripts-to-genes mapping file needs to be generated, set species_t2g = None, generate_index = True and species_gtf = "path/to/GTF file". |
| `species_fasta` (str) | Path to the transcriptome (cDNA) FASTA file for the species of interest. This will be required if ERCCs are included in the dataset |


#### Optional parameters
___________________

| Parameter | Description |
| --- | --- |
| `input_directory` (bool) | Indicate whether the list_of_fastqs parameter is given the path to the directory containing all input FastQ files * See example above (working with more than one set of FastQ files) |
| `read_separator` (str-list) | The strings used to separate the input reads. Within a list - the first element should contain the string used to identify read 1, the second element should contain the string used to identify read 2 and if "10xv1" FastQ files are used, a third element is required that contains the string used to identify read 3. * See example above (working with more than one set of FastQ files) |
| `generate_index` (bool) | Indicate whether the Kallisto Index should be generated within the function. If set to True: the species_fasta and the species_gtf are required |
| `species_gtf` (str) | Path to the GTF file for the species of interest. This will be used to create the transcripts-to-genes mapping file. If generate_index is set to True, this parameter is required. |
| `k_mer_length` (int) | The length of the K-mers that should be generated when creating the Kallisto index |
| `intron` (bool) | Indicate whether or not to include intron transcript ids |
| `filter` (bool) | Indicate whether or not to filter the BUS file prior to generating the gene-count matrix in mtx format |
| `UMI_bp` (str) | The number of base pairs sequenced for the UMI sequence. If 10xv1 technology is used, this parameter is required |
| `barcode_bp` (str) | The number of base pairs sequenced for the barcode sequence. If 10xv1 technology is used, this parameter is required |
| `transcript_bp` (str) | The number of base pairs sequenced for the transcript sequence |
| `whitelist_path` (str) | Path to the barcode whitelist that will be used for barcode correction |
| `path_to_prefix_count_files` (str) | Prefix of the output matrix files and indices |
| `memory` (str) | Amount of memory to use |

* **Kallisto specified FastQ file order and single-cell technologies**

Please see table above

* **Working with more than one set of fastq files**

Please see description above


## `seurat_matrix`

This SASCRiP function is designed to check the input matrix and index files to make sure it is in the correct format for use with Seurat. If the input files are not in the correct format `seurat_matrix` rearranges the files and saves the new files in the specified output directory.  

#### Usage

````python
import sascrip
from sascrip import sascrip_functions

sascrip_functions.seurat_matrix(
     matrix_file,
     gene_index,
     barcode_index,
     output_directory
)

````
#### Required parameters
___________________

| Parameters | Description |
| --- | --- |
| `matrix_file` (str) | Path to the mtx matrix file |
| `gene_index` (str) | Path to the gene index file |
| `barcode_index` (str) | Path to the barcode index file |
| `output_directory` (str) | Path to the output directory where the new matrix files will be saved, if new files are generated |



## `run_cqc`

This is the main cell quality control function. `run_cqc` runs Seurat functions on the input Seurat-compatible mtx matrix and performs cell quality control through the use of filtering thresholds. Multiple parameters are included within this function to allow for many different types of inputs and outputs for integration with many other single-cell analysis tools

#### Usage

````python
import sascrip
from sascrip import sascrip_functions

sascrip_functions.run_cqc(
     input_file_or_folder,
     sample_ID,
     output_directory_path = "working_directory,
     generate_seurat_object = True,
     subset_seurat_object = True,
     generate_default_plots = True,
     gene_column = 1,
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
)

````
#### Required parameters
___________________

| Parameter | Description |
| --- | --- |
| `input_file_or_folder` (str) | Path to the folder containing the Seurat-compatible matrices or path to the hdf5 file. Additionally, saved seurat objects in rds or rdata format can be input here. If a Seurat object is input, the input_seurat_object parameter must be set to True |
| `sample_ID` (str) | The name of the sample |

#### Optional parameters
___________________

| Parameter | Description |
| --- | --- |
| `output_directory_path` (str) | The path to the output directory where all output files and directories will be created and saved |
| `generate_seurat_object` (bool) | Indicate whether a seurat object should be generated from the input mtx matrix. If a Seurat object is input in the input_file_or_folder parameter, this parameter is required to be set to False. If a gene count matrix is input, this parameters is required to be set to True |
| `subset_seurat_object` (bool) | Indicate whether the Seurat object should be subset, removing low-quality cells identified by given thresholds |
| `generate_default_plots` (bool) | Indicate whether deafault plots should be generated to visualise the single-cell data and identified low-quality cells |
| `gene_column` (int) | The column number in the genes index file that should be used for the row names (gene names) in the Seurat object. 1-based |
| `input_seurat_object` (bool) | Indicate whether the input_file_or_folder parameter contains the path to a saved Seurat object. If so, this parameter should be set to True |
| `transcripts_to_genes_files` (str) | Path to the transcripts-to-genes mapping file that will allow ENSG gene names to be converted into corresponding HGNC gene symbols (within the seurat object) if required |
| `gene_lower` (int/None) | Minimum number of genes that should be detected in healthy cells. If this cell metric should not be used to identify low-quality cells then gene_lower should be set to None. However, if this parameter is set to None and generate_default_plots is set to True, a warning will be returned as the visualisations that use this threshold cannot be generated |
| `gene_higher_method` (str) | One of three methods - "MAD" (Median Absolute Deviation), "SD" (Standard Deviation), or "Manual" - that should be used to identify outliers using total gene count per cell. If Manual is selected, a value must be given for the gene_higher parameter |
| `gene_higher` (int/None) | Maximum number of genes that should be detected in single cells. If the gene_higher_method is set to "Manual", a value for gene_higher must be given. If this cell metric should not be used to identify cell doublets then gene_higher should be set to None. However, if this parameter is set to None and generate_default_plots is set to True, a warning will be returned as the visualisations that use this threshold cannot be generated |
| `mitochondria_percent` (int/None) | The maximum percentage of mitochondrial genes within a cell. If this cell metric should not be used to identify damaged cells then mitochondria_percent should be set to None. However, if this parameter is set to None and generate_default_plots is set to True, a warning will be returned as the visualisations that use this threshold cannot be generated |
| `nMADs` (int) | If the "MAD" method is selected in gene_higher_method, the value given to nMADs is used as the threshold to classify cells as outliers |
| `nSD` (int) | If the "SD" method is selected in gene_higher_method, the value given to nSD is used as the threshold to classify cells as outliers |
| `extract_cell_metrics` (bool) | Indicate whether to extract the calculated cell metrics as a .tsv file |
| `output_matrix (bool)` | Indicate whether to generate an mtx matrix file from the filtered Seurat data. This matrix will contain all healthy single cells characterised by the given thresholds |



## `stransform_normalize`

`sctransform_normalize` takes the raw UMI counts from healthy single cells (stored within a seurat object) and generates gene expression values using sctransform through Seurat. A saved Seurat object containing log normalised expression values as well as corrected counts. Additionally, the top 2000 highly variable genes are returned. This information is required for downstream analysis such as clustering or data integration. 

#### Usage

````python
import sascrip
from sascrip import sascrip_functions

sascrip_functions.sctransform_normalize(
     seurat_object,
     sample_ID,
     output_directory_path = "working_directory,
     output_log_matrix = False,
     output_count_matrix = False,
     transcripts_to_genes_file = None,
     **additional_sctransform_arguments
)

````
#### Required parameters
_________________________

| Parameter | Description |
| --- | --- |
| `seurat_object` (str) | Path to the saved filtered Seurat object |
| `sample_ID` (str) | Name of sample |


#### Optional parameters
________________________

| Parameter | Description |
| --- | --- |
| `output_directory_path` (str) | Path to the output directory where all generated files and dircetories will be saved |
| `output_log_matrix` (bool) | Indicate whether to additionally store the gene expression values per cell in an mtx matrix |
| `output_count_matrix` (bool) | Indicate whether to additionally store the corrected UMI counts in an mtx matrix |
| `transcripts_to_genes_files` (str) | Path to the transcripts-to-genes mapping file that will allow ENSG gene names to be converted into corresponding HGNC gene symbols (within the seurat object) if required |
| **additional_sctransform_arguments (dict) | Additional parameters (with key words) that should be passed to Seurat's SCTransform function - which additionally passes the parameters to the original stransform::vst function. In order to use parameter - the additional parameters to be passed should be in the form of a python dictionary where the key word is the parameter name (str) and the value is the given parameter value (in it's correct data type) - *see example below |

* Example when using additional_sctransform_arguments - To access the SCTransform parameter "conserve.memory" and set it to True:

```python

import sascrip
from sascrip import sascrip_functions

# Set up the dictionary for the additional parameters
additional_parameters = {"conserve.memory": True} # Note: If the parameter requires a boolean type, please use Pythons boolean datatype (True/False)

# Run sctransform_normalize with additional parameters
sascrip_functions.sctransform_normalize(
     seurat_object,
     sample_ID,
     **additional_parameters # Indicate the keyword dictionary using "**"
)

```


## `sascrip_preprocess`

`sascrip_preprocess` allows the user to run the entire single-cell RNA sequencing data pre-processing steps with one function. The sascrip_preprocess parameters can be customised to adjust the default settings. 

#### Usage

````python
import sascrip
from sascrip import sascrip_functions

sascrip_functions.sascrip_preprocess(
     output_directory_path,
     sample_ID,
     list_of_fastqs,
     single_cell_technology,
     species_index,
     species_t2g,
     input_directory = False,
     read_separator = None,
     filter = True,
     include_checkpoints = False,
     kallisto_bustools_count_parameters = None,
     run_cqc_parameters = None,
     additional_sctransform_arguments = None
)

````
#### Required parameters
___________________

| Parameter | Description |
| --- | --- |
| `output_directory` (str) | Path to the output directory where output files will be saved |
| `sample_ID` (str) | Name of the sample |
| `list_of_fastqs` (str-list) | Python list of the paths to input FastQ files in the *order specified by Kallisto. The folder containing all the input files can be given instead and the relevant FastQ files will be captured and sorted in the correct order. To use this feature - input the path to the directory here (If multiple directories are used, input all the directories as a list), set input_directory = True and provide the strings used to separate the reads in the read_separator parameter. * See example above (working with more than one set of FastQ files) |
| `single_cell_technology` (str) | The single-cell sequencing technology that was used as *specified by Kallisto. If 10xv1 technology was used, the UMI_bp and barcode_bp parameters are required |
| `species_index` (str) | Path to the kallisto_index for the species of interest. If no index is given (species_index = None), the default Kallisto index created using the GRCh 38 transcriptome assembly will be used. If the kallisto_index needs to be generated, set kallisto_index = None and create a keyword dictionary (as described for **additional_sctransform_arguments) for the following kallisto_bustools_count parameters; "generate_index", "species_fasta", "species_gtf". This dictionary should be created for kallisto_bustools_count_parameters |
| `species_t2g` (str) |  Path to the transcript-to-genes mapping file for the species of interest. If no mapping file is given, the default file created using the GRCh 38 GTF file will be used. If the transcripts-to-genes mapping file needs to be generated, set species_t2g = None. The dictionary supplied to kallisto_bustools_count_parameters will allow the kallisto_bustools_count function to create both the kallisto index and the transcripts-to-genes mapping file |        

#### Optional parameters
___________________

| Parameter | Description |
| --- | --- |
| `input_directory` (bool) | Indicate whether the list_of_fastqs parameter is given the path to the directory containing all input FastQ files * See example above (working with more than one set of FastQ files) |
| `read_separator` (str-list) | The strings used to separate the input reads. Within a list - the first element should contain the string used to identify read 1, the second element should contain the string used to identify read 2 and if "10xv1" FastQ files are used, a third element is required that contains the string used to identify read 3. * See example above (working with more than one set of FastQ files) |
| `filter` (bool) | Indicate whether to filter the BUS file prior to generating the count matrix |
| `include_checkpoints` (bool) | Indicate whether to print all log statements to standard output |
| `kallisto_bustools_count_parameters` (dict) | Any additional parameters to be adjusted for the kallisto_bustools_count function in the same format as described above for **additional_sctransform_arguments. In this case the ** is not needed |
| `run_cqc_parameters` (dict) | Any additional parameters to be adjusted for the run_cqc functions in the same format as described above for **additional_sctransform_arguments. In this case the ** is not needed |
| `additional_sctransform_arguments` (dict) | Additional parameters (with key words) that should be passed to Seurat's SCTransform function - which additionally passes the parameters to the original stransform::vst function. This dictionary should be created in the same format as previously described for this parameter. In this case the ** is not needed |


* **Kallisto specified FastQ file order and single-cell technologies**

Please see table above

* **Working with more than one set of fastq files**

Please see description above

## Authors

#### Darisia Moonsamy
email: darisia@outlook.com <br>
Institution: University of the Witwatersrand, Johannesburg, South Africa <br>
Department: School of Molecular and Cell Biology <br>
Research areas: Bioinformatics (genetics and immunology) <br>

#### Dr Nikki Gentle


## License

[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)


## Reference


## Read next

#### [Seurat vignette: A guided clustering workflow](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) 

#### [Kallisto documentation](https://pachterlab.github.io/kallisto/manual)

#### [BUStools documentation](https://bustools.github.io/manual)
