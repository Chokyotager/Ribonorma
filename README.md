# Ribonorma
As part of the paper **Transcripts Per Million Ratio: a novel batch and sample control methodover an established paradigm** by H. Lam and R. Pincket.


## Installation
Installation can be done with Pip. Python 3.6+.

`pip install ribonorma`
`pip3 install ribonorma`

## How to use

### Command line usage
![Pipeline](/images/Ribonorma%20analysis%20pipeline.png)

1. Suppose you have two files, one file with the raw RNASeq count data ([example](/test/maqc_count_with_lengths.tsv)) and one file with the phenotype file data ([example](/test/maqc_phenotypes.tsv)).

2. Run `ribonorma-normalise`

### Python code usage
```py
import csv
from ribonorma import ribonorma

reads = list(csv.read( ... )) # Reads
gene_length = [] # List of gene lengths

conditions = ["Standard media", "Standard media", "Standard media", "Super media", "Super media", "Super media"] # Example conditions

normalised_counts = ribonorma.tpmr(reads, gene_length, conditions, percent_housekeep=10)
```


`ribonorma.tpm(reads, gene_length)`
* `reads`: 1D list of read counts
* `gene_length`: 1D list of individual gene lengths


`ribonorma.tpmm(samples, gene_length)`
* `samples`: 2D list of sample read counts - [sample x count]
* `gene_length`: 1D list of individual gene lengths


`ribonorma.tpmr(samples, gene_length, experimental_conditions, percent_housekeep=10)`
* `samples`: 2D list of sample read counts - [sample x count]
* `gene_length`: 1D list of individual gene lengths
* `experimental_conditions`: 1D list of individual experimental conditions, same size as `samples`
* `alpha`: floating point value of percent housekeep as stated in the paper


`ribonorma.tpmr_2(samples, gene_length, experimental_conditions, percent_housekeep=10)`
* `samples`: 2D list of sample read counts - [sample x count]
* `gene_length`: 1D list of individual gene lengths
* `experimental_conditions`: 1D list of individual experimental conditions, same size as `samples`
* `alpha`: floating point value of percent housekeep as stated in the paper
