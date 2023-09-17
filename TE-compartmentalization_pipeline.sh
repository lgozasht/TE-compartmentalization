########
#Usage:
# sh TE-compartmentalization_pipeline.sh {species_name} {genome.fasta} {TE_proteins.fa} {go-obo file} {percentile}
#######
species=$1
genomeFasta=$2
TE_proteins=$3
goObo=$4
percentile=$5



#Run Repeatmodeler2
mine_repeats () {
  #$1 = species name (no spaces)
  #$2 = genome fasta file
  if test -f "$1-families.fa"; then
      echo "$1-families.fa exists."

  else 
      perl BuildDatabase  -name $1 $2
      perl RepeatModeler -database $1
      echo "$1-families.fa does not exist."
  fi
}

#Filter spurious TE families from repeatmodeler
#Produces final TE library: consensi.fa.classified.filtered_for_CDS_repeats.fa
#Also produces file containing non-TE-associated transcripts: transcripts.no_tes.fa
filter_repeats () {
  #$1 = species name (no spaces)
  #$2 = path to TE protein fasta file

  if test -f "transcripts.fa"; then
    echo "no transcript file, skipping filtering step"

  else
    echo "filtering"
    #make blast database for TE protein library
    makeblastdb -in $2 -dbtype prot
  
    #blast species proteins against TE database
    blastp -query proteins.fa \
         -db $2 \
         -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
         -max_target_seqs 25 \
         -culling_limit 2 \
         -evalue 1e-5 \
         -out proteins.fa.vs.RepeatPeps.25cul2.1e5.blastp.out
  
    #Remove TE-associated 
    python ../util/FilterTEsFromTranscripts.py
  
    #make blast db for species transcripts minus those associated with TEs
    makeblastdb -in transcripts.no_tes.fa -dbtype nucl
  
    #blast TE library at transcripts
    blastn -task megablast \
         -query $1-families.fa \
         -db transcripts.no_tes.fa \
         -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
         -max_target_seqs 25 \
         -culling_limit 2 \
         -evalue 1e-25 \
         -out repeatmodeller_lib.vs.transcripts.no_tes.25cul2.1e10.megablast.out
  
    #Remove unknown TE families with hits to non-TE-associated species proteins 
    fastaqual_select.pl -f $1-families.fa \
         -e <(awk '{print $1}' repeatmodeller_lib.vs.transcripts.no_tes.25cul2.1e10.megablast.out | grep -i "Unknown" | sort | uniq) > consensi.fa.classified.filtered_for_CDS_repeats.fa
  fi
}



###Annotate repeats
annotate_repeats () {
  #$1 = genome.fasta
  if test -f "repeats/$1.tbl"; then
      echo "$1 exists."
  else
      perl RepeatMasker -lib consensi.fa.classified.filtered_for_CDS_repeats.fa -excln -s -no_is -dir repeats -u -noisy -html -xm -a -xsmall -source $1

  fi
}

get_repeat_divergence () {
  #$1 = species_name (no spaces)
  #$2 = genome.fasta
  calcDivergenceFromAlign.pl -s repeats/$1.divsum repeats/$2.align
}



#
prepare_analysis_all () {
  python ../util/clean_RM_out.py repeats/$1.out 
  python RM2Bed.py -o higher_score repeats/$1.out2  allRepeats.bed

  #Filter repeat related genes from gene file
  python ../util/getGenesWithoutTEs.py


  ###Make 100kb windows for genome
  bedtools makewindows -g scaffolds.bed -w 50000 > windows.bed
  bedtools intersect -wa -u -a windows.bed -b geneFile_final.tsv > windows_filtered.bed 

  ###Analyze each repeat class independently
  grep "LTR" allRepeats.bed > LTRs.bed
  grep "DNA" allRepeats.bed > DNA_TEs.bed
  grep -i "Simple\|sat\|Low_complexity" allRepeats.bed > Simple.bed
  grep "LINE" allRepeats.bed > LINE.bed
  grep "SINE" allRepeats.bed > SINE.bed
  grep -v -i "Simple\|sat\|Low_complexity\|rRNA\|snRNA\|tRNA" allRepeats.bed > onlyTEs.bed

  #Make sure to filter out scaffolds shorter than 100000
  bedtools coverage -a windows_filtered.bed -b allRepeats.bed | cut -f 1,2,3,7 | awk '$3-$2 >=50000' > repeat_coverage.bed
  bedtools coverage -a windows_filtered.bed -b LTRs.bed | cut -f 1,2,3,7 | awk '$3-$2 >=50000' > LTR_repeat_coverage.bed
  bedtools coverage -a windows_filtered.bed -b DNA_TEs.bed | cut -f 1,2,3,7 | awk '$3-$2 >=50000' > DNA_repeat_coverage.bed
  bedtools coverage -a windows_filtered.bed -b onlyTEs.bed | cut -f 1,2,3,7 | awk '$3-$2 >=50000' > TE_repeat_coverage.bed
  bedtools coverage -a windows_filtered.bed -b Simple.bed | cut -f 1,2,3,7 | awk '$3-$2 >=50000' > Simple_repeat_coverage.bed
  bedtools coverage -a windows_filtered.bed -b LINE.bed | cut -f 1,2,3,7 | awk '$3-$2 >=50000' > LINE_repeat_coverage.bed
  bedtools coverage -a windows_filtered.bed -b SINE.bed | cut -f 1,2,3,7 | awk '$3-$2 >=50000' > SINE_repeat_coverage.bed



}





###Convert to bed
prepare_analysis_genes_and_flanks () {


  ###Extend boundaries of each gene
  mkdir gene_centric_analysis
  cd gene_centric_analysis

  awk '{print $1"\t"0"\t"$2}' ../scaffolds.bed > scaffolds_with_starts.bed
  bedtools intersect -wa -u -a ../geneFile_final.tsv -b scaffolds_with_starts.bed > geneFile_final_primary_scaffolds.bed
  bedtools slop -b 50000 -i geneFile_final_primary_scaffolds.bed -g ../scaffolds.bed > geneFile_extended.bed


  #For now lets simplify
  #find repeat coverage across windows
  #Make sure to filter out scaffolds shorter than 100000
  bedtools coverage -a geneFile_extended.bed -b ../allRepeats.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}'  > repeat_coverage.bed
  bedtools coverage -a geneFile_extended.bed -b ../LTRs.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}'  > LTR_repeat_coverage.bed
  bedtools coverage -a geneFile_extended.bed -b ../DNA_TEs.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}' > DNA_repeat_coverage.bed
  bedtools coverage -a geneFile_extended.bed -b ../onlyTEs.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}' > TE_repeat_coverage.bed
  bedtools coverage -a geneFile_extended.bed -b ../Simple.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}'  > Simple_repeat_coverage.bed
  bedtools coverage -a geneFile_extended.bed -b ../LINE.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}'  > LINE_repeat_coverage.bed
  bedtools coverage -a geneFile_extended.bed -b ../SINE.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}'  > SINE_repeat_coverage.bed


}



#### Find genes enriched in repeat hotspots within 95, 97, and 99 percentile
run_analysis_genes_and_flanks () {
#Gets 95 percentile regions
  python ../../util/repeatDensity_gene_specific.py repeat_coverage.bed all $3
  python ../../util/repeatDensity_gene_specific.py LTR_repeat_coverage.bed LTR $3
  python ../../util/repeatDensity_gene_specific.py DNA_repeat_coverage.bed DNA $3
  python ../../util/repeatDensity_gene_specific.py TE_repeat_coverage.bed TE $3
  python ../../util/repeatDensity_gene_specific.py Simple_repeat_coverage.bed simple $3 
  python ../../util/repeatDensity_gene_specific.py LINE_repeat_coverage.bed LINE $3
  python ../../util/repeatDensity_gene_specific.py SINE_repeat_coverage.bed SINE $3


  cut -f 5 all_hotspots_${3}.bed > all_enriched_${3}.list
  cut -f 5 LTR_hotspots_${3}.bed > LTR_enriched_${3}.list
  cut -f 5 DNA_hotspots_${3}.bed > DNA_enriched_${3}.list
  cut -f 5 TE_hotspots_${3}.bed > TE_enriched_${3}.list
  cut -f 5 simple_hotspots_${3}.bed > simple_enriched_${3}.list
  cut -f 5 LINE_hotspots_${3}.bed > LINE_enriched_${3}.list
  cut -f 5 SINE_hotspots_${3}.bed > SINE_enriched_${3}.list
  


  comm -12 <(sort all_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > all_enriched.list_${3}.filtered
  comm -12 <(sort DNA_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > DNA_enriched.list_${3}.filtered
  comm -12 <(sort LTR_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > LTR_enriched.list_${3}.filtered
  comm -12 <(sort TE_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > TE_enriched.list_${3}.filtered
  comm -12 <(sort simple_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > simple_enriched.list_${3}.filtered
  comm -12 <(sort LINE_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > LINE_enriched.list_${3}.filtered
  comm -12 <(sort SINE_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > SINE_enriched.list_${3}.filtered
  

  if [[ -s all_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
         	echo "enriched.list.filtered	is populated"
          
          python find_enrichment.py  all_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0 --obo $2  --outfile all_geneRepeatEnrichment_${3}.out
          
  else
          # The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  
  if [[ -s DNA_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
          
  	echo "enriched.list.filtered    is populated"
  
          python find_enrichment.py  DNA_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0  --obo $2   --outfile DNA_geneRepeatEnrichment_${3}.out
          
  else
      	# The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  
  if [[ -s LTR_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
          
  	echo "enriched.list.filtered    is populated"
  
          python find_enrichment.py  LTR_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0 --obo $2   --outfile LTR_geneRepeatEnrichment_${3}.out
          
  else
      	# The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  
  if [[ -s TE_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
          
  	echo "enriched.list.filtered    is populated"
  
          python find_enrichment.py  TE_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0 --obo $2   --outfile TEs_geneRepeatEnrichment_${3}.out
          
  else
      	# The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  
  
  if [[ -s simple_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
          
  	echo "enriched.list.filtered    is populated"
  
          python find_enrichment.py  simple_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0 --obo $2   --outfile simple_geneRepeatEnrichment_${3}.out
          
  else
      	# The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  
  
  if [[ -s LINE_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
          
  	echo "enriched.list.filtered    is populated"
  
          python find_enrichment.py  LINE_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0 --obo $2   --outfile LINE_geneRepeatEnrichment_${3}.out
          
  else
      	# The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  
  if [[ -s SINE_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
          
          echo "enriched.list.filtered    is populated"
  
          python find_enrichment.py  SINE_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0 --obo $2   --outfile SINE_geneRepeatEnrichment_${3}.out
          
  else
          # The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  
cd ..
}



prepare_analysis_only_flanks () {


  ###Extend boundaries of each gene
  mkdir gene_centric_flank_analysis
  cd gene_centric_flank_analysis
  
  
  awk '{print $1"\t"0"\t"$2}' ../scaffolds.bed > scaffolds_with_starts.bed
  bedtools intersect -wa -u -a ../geneFile_final.tsv -b scaffolds_with_starts.bed | grep -v "NT_\|NW_" > geneFile_final_primary_scaffolds.bed
  bedtools flank -b 50000 -i geneFile_final_primary_scaffolds.bed -g ../scaffolds.bed > geneFile_extended.bed
  
  bedtools sort -i geneFile_final_primary_scaffolds.bed > geneFile_final_primary_scaffolds.sorted.bed
  bedtools merge -i geneFile_final_primary_scaffolds.sorted.bed > geneFile_final_primary_scaffolds.merged.bed
  bedtools coverage -wo -a geneFile_extended.bed -b geneFile_final_primary_scaffolds.merged.bed > geneFile_coding_flanks.bed
  
  
  
  
  
  #For now lets simplify
  #find repeat coverage across windows
  #Make sure to filter out scaffolds shorter than 100000
  bedtools coverage -a geneFile_extended.bed -b ../allRepeats.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}'  > repeat_coverage.bed.pre
  bedtools coverage -a geneFile_extended.bed -b ../LTRs.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}'  > LTR_repeat_coverage.bed.pre
  bedtools coverage -a geneFile_extended.bed -b ../DNA_TEs.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}' > DNA_repeat_coverage.bed.pre
  bedtools coverage -a geneFile_extended.bed -b ../onlyTEs.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}' > TE_repeat_coverage.bed.pre
  bedtools coverage -a geneFile_extended.bed -b ../Simple.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}'  > Simple_repeat_coverage.bed.pre
  bedtools coverage -a geneFile_extended.bed -b ../LINE.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}' > LINE_repeat_coverage.bed.pre
  bedtools coverage -a geneFile_extended.bed -b ../SINE.bed | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$4}' > SINE_repeat_coverage.bed.pre
  
  ###Intersect flanks with other genes and remove coding basepairs
  
  python ../../util/combine_flanks.py repeat_coverage.bed.pre
  python ../../util/combine_flanks.py LTR_repeat_coverage.bed.pre
  python ../../util/combine_flanks.py DNA_repeat_coverage.bed.pre
  python ../../util/combine_flanks.py TE_repeat_coverage.bed.pre
  python ../../util/combine_flanks.py Simple_repeat_coverage.bed.pre
  python ../../util/combine_flanks.py LINE_repeat_coverage.bed.pre
  python ../../util/combine_flanks.py SINE_repeat_coverage.bed.pre
  
}



#### Find genes enriched in repeat hotspots within 95, 97, and 99 percentile
run_analysis_only_flanks () {
  #Gets 95 percentile regions
  python ../../util/repeatDensity_gene_specific.py repeat_coverage.bed all $3
  python ../../util/repeatDensity_gene_specific.py LTR_repeat_coverage.bed LTR $3
  python ../../util/repeatDensity_gene_specific.py DNA_repeat_coverage.bed DNA $3
  python ../../util/repeatDensity_gene_specific.py TE_repeat_coverage.bed TE $3
  python ../../util/repeatDensity_gene_specific.py Simple_repeat_coverage.bed simple $3 
  python ../../util/repeatDensity_gene_specific.py LINE_repeat_coverage.bed LINE $3
  python ../../util/repeatDensity_gene_specific.py SINE_repeat_coverage.bed SINE $3
  
  
  cut -f 5 all_hotspots_${3}.bed > all_enriched_${3}.list
  cut -f 5 LTR_hotspots_${3}.bed > LTR_enriched_${3}.list
  cut -f 5 DNA_hotspots_${3}.bed > DNA_enriched_${3}.list
  cut -f 5 TE_hotspots_${3}.bed > TE_enriched_${3}.list
  cut -f 5 simple_hotspots_${3}.bed > simple_enriched_${3}.list
  cut -f 5 LINE_hotspots_${3}.bed > LINE_enriched_${3}.list
  cut -f 5 SINE_hotspots_${3}.bed > SINE_enriched_${3}.list
  
  
  comm -12 <(sort all_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > all_enriched.list_${3}.filtered
  comm -12 <(sort DNA_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > DNA_enriched.list_${3}.filtered
  comm -12 <(sort LTR_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > LTR_enriched.list_${3}.filtered
  comm -12 <(sort TE_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > TE_enriched.list_${3}.filtered
  comm -12 <(sort simple_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > simple_enriched.list_${3}.filtered
  comm -12 <(sort LINE_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > LINE_enriched.list_${3}.filtered
  comm -12 <(sort SINE_enriched_${3}.list) <(sort ../genes_with_terms_final.txt) > SINE_enriched.list_${3}.filtered
  

  if [[ -s all_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
          
         	echo "enriched.list.filtered	is populated"
          
          python find_enrichment.py  all_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0 --obo $2   --outfile all_geneRepeatEnrichment_${3}.out
          
  else
          # The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  
  if [[ -s DNA_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
          
  	echo "enriched.list.filtered    is populated"
  
          python find_enrichment.py  DNA_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0 --obo $2   --outfile DNA_geneRepeatEnrichment_${3}.out
          
  else
      	# The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  
  if [[ -s LTR_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
          
  	echo "enriched.list.filtered    is populated"
  
          python find_enrichment.py  LTR_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0 --obo $2   --outfile LTR_geneRepeatEnrichment_${3}.out
          
  else
      	# The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  
  if [[ -s TE_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
          
  	echo "enriched.list.filtered    is populated"
  
          python find_enrichment.py  TE_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0 --obo $2   --outfile TEs_geneRepeatEnrichment_${3}.out
          
  else
      	# The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  

  if [[ -s simple_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
          
  	echo "enriched.list.filtered    is populated"
  
          python find_enrichment.py  simple_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0 --obo $2   --outfile simple_geneRepeatEnrichment_${3}.out
          
  else
      	# The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  
  
  if [[ -s LINE_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
          
  	echo "enriched.list.filtered    is populated"
  
          python find_enrichment.py  LINE_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0 --obo $2   --outfile LINE_geneRepeatEnrichment_${3}.out
          
  else
      	# The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  

  if [[ -s SINE_enriched.list_${3}.filtered ]]; then
          # The file is not-empty.
          
          echo "enriched.list.filtered    is populated"
  
          python find_enrichment.py  SINE_enriched.list_${3}.filtered ../all_genes.list_final.filtered ../gene2GO_final.tsv --method fdr_bh --pval=1.0 --obo $2   --outfile SINE_geneRepeatEnrichment_${3}.out
          
  else
          # The file is empty.
          echo "enriched.list.filtered is empty"
  fi
  
}



mine_repeats $species $genomeFasta

filter_repeats $species $TE_proteins

annotate_repeats $genomeFasta

get_repeat_divergence $species $genomeFasta

prepare_analysis_all $genomeFasta

prepare_analysis_genes_and_flanks 
run_analysis_genes_and_flanks $species $goObo $percentile


prepare_analysis_only_flanks 
run_analysis_only_flanks $species $goObo $percentile




