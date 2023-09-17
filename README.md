# TE-compartmentalization

This repository contains the code associated with Gozashti and Corbett-Detig's **Universal signatures of transposable element compartmentalization across eukaryotic genes**. 

## TE-compartmentalization pipeline

We define "TE-compartmentalized genes" as genes whose flanking regions show TE density within the top 90th percentile when all genes are compared for a given species. This pipeline only requires a genome sequence file in fasta format and cooresponding gene annotation in bed format. You may also want to include a fasta file of TE-associated proteins to filter out potential spurious TE families. A basic library of TE proteins can be obtained from RepeatMasker (RepeatMasker/Libraries/RepeatPeps.lib).

### Tool Dependencies

Each of these should be visible in your path. You can install most required tools using conda.

* RepeatModeler2
* RepeatMasker
* Blast
* fastaqual_select.pl from (https://github.com/sujaikumar/assemblage/tree/master)
* bedtools
* goatools

## Usage

```
sh TE-compartmentalization_pipeline.sh {species_name} {genome.fasta} {genes.bed}
```

## Wrapper to run pipeline on all refseq and genbank species



