# SASCRiP_draft_repo
This repo contains all the code and documentation for the SASCRiP package

## Table of Contents
1. [Overview](https://github.com/Darisia/SASCRiP_draft_repo/blob/main/README.md#overview)
2. [Installation](https://github.com/Darisia/SASCRiP_draft_repo/blob/main/README.md#installation)
3. [SASCRiP's workflow](https://github.com/Darisia/SASCRiP_draft_repo/blob/main/README.md#sascrips-workflow)
4. <details>
     <summary><a href="https://github.com/Darisia/SASCRiP_draft_repo/blob/main/README.md#sascrip-functions---user-guide">SASCRiP functions - User guide</a></summary>
     <br>
     <a href="https://github.com/Darisia/SASCRiP_draft_repo/blob/main/README.md#edit_10xv1_fastq">Edit_10v1_fastq</a>
     </details>


## Overview


## Installation


## SASCRiP's workflow


## SASCRiP functions - User guide

### run_fastqc

## edit_10xv1_fastq
edit_10xv1_fastq prepares FastQ files, obtained using the 10xv1 sequencing chemistry, for input into Kallisto for pseudoalignment and quantification. The directory containing the FastQ files are input into edit_10xv1_fastq which searches for the RA FastQ file that contains both the UMI and transcript sequence. The UMI and Transcript sequences are separated into their own FastQ files. The directory containing all these FastQ files can be input into kallisto_bustools_count for further processing

#### Usage

````python
import sascrip
from sascrip import sascrip_functions

sascrip_functions.edit_10xv1_fastq(input_directory, output_directory)

````
#### Parameters

````
Required parameters
___________________

input_directory (str): Path to the directory containing the RA 10xv1 FastQ files

output-directory (str): Path to the output directory where the new separated FastQ files will be saved

````    

## check_ercc

check_ercc allows the user to check if the single-cell dataset may contain RNA spike-ins that can be used as an additional cell-quality control metric. If check_ercc is run, True or False - depending on whether ERCCs are included or not, will be returned as standard output 

#### Usage

````python
import sascrip
from sascrip import sascrip_functions

sascrip_functions.check_ercc(
     ERCC_fasta,
     output_directory,
     list_of_fastqs,
     single_cell_technology,
     UMI_bp = '0',
     barcode_bp = '0',
     transcript_bp = '0'
)

````
#### Parameters

````
Required parameters
___________________

ERCC_fasta (str):             Path to the ERCC FASTA file

output_directory (str):       Path to the output directory where output files will be saved

list_of_fastqs (str-list):    Python list of the paths to input FastQ files in the *order specified by Kallisto

single_cell_technology (str): The single-cell sequencing technology that was used as *specified by Kallisto. If 10xv1 technology was used, the UMI_bp and barcode_bp parameters are required

Optional parameters
___________________

UMI_bp (str):                 The number of base pairs sequenced for the UMI sequence. If 10xv1 technology is used, this parameter is required

barcode_bp (str):             The number of base pairs sequenced for the barcode sequence. If 10xv1 technology is used, this parameter is required

transcript_bp (str):          The number of base pairs sequenced for the transcript sequence

```` 

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

Include the additional sets in the fastq list, keeping the sets together. For example: if you have 2 sets of fastq files from the 10xv2 chemistry
```python

list_of_fastqs = ["barcode-UMI_R1_1.fastq.gz", "Transcript_R2_1.fastq.gz", "barcode_UMI_R1_2.fastq.gz", "Transcript_R2_2.fastq.gz"]

```

## kallisto_bustools_count

## include_ERCC_bus_count

## run_cqc

## stransform_normalize

## sctransform_cell_cycle

## sascrip-preprocess


## Authors


## Licence


## Reference


## Read next
