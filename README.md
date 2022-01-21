# Ribonorma
As part of the paper


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
from ribonorma import ribonorma
```
