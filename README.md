# isoform-coevolution
Accompanying code for "Sites of transcription initiation drive mRNA isoform selection" by Alfonso-Gonzalez et al. 

Given a set of dominant and non.dominant GFF files specifying genes in `input/gff/`, this pipeline extracts MAFs for each gene and the corresponding FASTA files, and using prody computes the covariance at each nucleotide position across the multiple alignment. 

A set of figures are output with optional boxes around the set of regions promoter regions specified in `config/all_box_genes.csv`, as well as statistics determining the covariance and the likelihood that this covariance is different from the background. 

## Before running

Download MAFs from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/dm6/multiz27way/maf/) to `input/mafs`. 

## Run 

Edit run.sh to suit your computational resources, then execute run.sh.

