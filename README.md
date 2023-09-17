# TE-compartmentalization

This repository contains the code associated with Gozashti and Corbett-Detig's **Universal signatures of transposable element compartmentalization across eukaryotic genes**. 

## TE-compartmentalization pipeline

We define "TE-compartmentalized genes" as genes whose flanking regions show TE density within the top 90th percentile when all genes are compared for a given species. This pipeline requires a genome sequence file in fasta format and cooresponding gene annotations in bed format, coding trancripts in fasta format and proteins in fasta format, as well as gene ontology annotations for each gene. All required input files can be obtained for any genome on genbank or refseq. Go annotations can be obtained from uniprot. You may also want to include a fasta file of TE-associated proteins to filter out potential spurious TE families. A basic library of TE proteins can be obtained from RepeatMasker (RepeatMasker/Libraries/RepeatPeps.lib). We ran this pipeline systematically across all species on genbank/refseq, so we created a wrapper to produce all necessary inputs given a genbank or refseq genome (see more on this below). However, you can also run it a standalone.

### Tool dependencies

Each of these should be visible in your path. 

* RepeatModeler2
* RepeatMasker
* RepeatMasker utils (PathToRepeatMakser/util)
* Blast
* fastaqual_select.pl from (https://github.com/sujaikumar/assemblage/tree/master)
* bedtools
* goatools
* bioawk

### Data dependencies

* go-obo file (can be obtained from https://geneontology.org/docs/download-ontology/#go_basic)
* Library of TE proteins (such as RepeatMasker/Libraries/RepeatPeps.lib)

### Standalone usage

Make a new directory and copy TE-compartmentalization_pipeline.sh to that directory. Then do the following in that directory.
1. Create a copy of your genome fasta file
2. Move the corresponding transcript fasta file to that directory and name it transcripts.fa
3. Move the corresponding protein fasta file to that directory and name it proteins.fa
4. Create a bed file for all genes in the genome named "geneFile_final.tsv". This file should have six columns. For each gene the first four columns should contain the scaffold, start, stop and gene name. The last two columns are not used in the analysis and can simply be "."
5. Identical files (only in the standalone case) named "all_genes.list_final.filtered" and "genes_with_terms_final.txt", which can be generated with the commands ``` cut -f 4 geneFile_final.tsv > all_genes.list_final.filtered ``` and ``` cut -f 4 geneFile_final.tsv > genes_with_terms_final.txt ```
6. A file named "gene2GO_final.tsv" with GO annotations for each gene. This file contains two columns. For each gene, the first column gives the gene name (must be consistent with other files) and the second column contains GO annotations for that gene. If a gene has multiple GO annotations, they should be separated by ";" (e.g. GO:111111;GO:222222;GO:333333).
7. A file detailing scaffold lengths for scaffolds longer than 50kb, called "scaffolds.bed". This file can be generated with the command: ``` bioawk -c fastx '{print $name"\t"length($seq)}' genome.fasta | awk '$2>50000' > scaffolds.bed ```

Then, run the script, specifying the species name (no spaces), genome fasta file name, full path to the TE protein database, full path to the GO obo file, and the percentile cutoff you want to use for calling TE-compartmentalized genes. We tried a variety of different cutoffs but settled on 90.

```
sh TE-compartmentalization_pipeline.sh {species_name} {genome.fasta} {PATH/TO/TE_proteins.fasta} {PATH/TO/gofile.obo} {percentile}
```

### Running the pipeline on all refseq and genbank species

We created a wrapper to run the pipeline systematically on all (or a subset of) genomes on RefSeq/GenBank. The wrapper 


