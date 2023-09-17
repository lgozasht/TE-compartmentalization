# TE-compartmentalization

This repository contains the code associated with Gozashti and Corbett-Detig's **Universal signatures of transposable element compartmentalization across eukaryotic genes**. 

## TE-compartmentalization pipeline

We define "TE-compartmentalized genes" as genes whose flanking regions show TE density within the top 90th percentile when all genes are compared for a given species. This pipeline requires a genome sequence file in fasta format and cooresponding gene annotations in bed format, coding trancripts in fasta format and proteins in fasta format, as well as gene ontology annotations for each gene. All required input files can be obtained for any genome on genbank or refseq. Go annotations can be obtained from uniprot. You may also want to include a fasta file of TE-associated proteins to filter out potential spurious TE families. A basic library of TE proteins can be obtained from RepeatMasker (RepeatMasker/Libraries/RepeatPeps.lib). We ran this pipeline systematically across all species on genbank/refseq, so we created a wrapper to produce all necessary inputs given a genbank or refseq genome (see more on this below). However, you can also run it a standalone.

### Dependencies

Each of these should be visible in your path. You can install most required tools using conda.

* RepeatModeler2
* RepeatMasker
* RepeatMasker utils (PathToRepeatMakser/util)
* Blast
* fastaqual_select.pl from (https://github.com/sujaikumar/assemblage/tree/master)
* bedtools
* goatools
* go-obo file (can be obtained from https://geneontology.org/docs/download-ontology/#go_basic)

### Standalone usage
Make a new directory and copy TE-compartmentalization_pipeline.sh to that directory. Then do the following in that directory.
1. Create a copy of your genome fasta file
2. Move the corresponding transcript fasta file to that directory and name it transcripts.fa
3. Move the corresponding protein fasta file to that directory and name it proteins.fa
4. Create a bed file for all genes in the genome named "geneFile_final.tsv". This file should have six columns. For each gene the first four columns should contain the scaffold, start, stop and gene name. The last two columns are not used in the analysis and can simply be ".".
5. 

```
sh TE-compartmentalization_pipeline.sh {species_name} {genome.fasta}  {TE_proteins.fa} {go-obo file} {percentile}
```

## Wrapper to run pipeline on all refseq and genbank species



