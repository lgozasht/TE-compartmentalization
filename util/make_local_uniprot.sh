#!/bin/bash

pwd
IN=$1
arrIN=(${IN//_/ })
echo ${arrIN[0]}"_"${arrIN[1]}



python ../util/getData_species.py ${arrIN[0]}"_"${arrIN[1]} > uniprot_all_species.tsv

#gene name
grep "${arrIN[0]}" uniprot_all_species.tsv | grep -i -f all_genes.list | awk -F"\t" '{print $4"\t"$6}' > localUniprot.tsv

#gene ID
grep "${arrIN[0]}" uniprot_all_species.tsv | grep -f all_genesID.list | awk -F"\t" '{gsub(/;/," ", $11); print $11"\t"$6}' > localUniprotID1.1.tsv
grep "${arrIN[0]}" uniprot_all_species.tsv | grep -f all_genesID.list | awk -F"\t" '{gsub(/;/," ", $8); print $8"\t"$6}' > localUniprotID1.2.tsv
cat localUniprotID1.1.tsv localUniprotID1.2.tsv > localUniprotID.tsv

#genbank
grep "${arrIN[0]}" uniprot_all_species.tsv | grep -i -f all_genes_genbank.list | awk -F"\t" '{print $4"\t"$6}' > localUniprot_genbank.tsv
cut -d"." -f 1 all_genes_genbank.list > all_genes_genbank.list_edit
grep "${arrIN[0]}" uniprot_all_species.tsv | grep -i -f all_genes_genbank.list_edit | awk -F"\t" '{print $4"\t"$6}' >> localUniprot_genbank.tsv
cat all_genes_genbank.list_edit >> all_genes_genbank.list


python ../util/correctDatabase.py commonName
python ../util/correctDatabase.py ID
python ../util/correctDatabase.py genbank




cut -f 1 gene2GO.tsv | sort | uniq > genes_with_terms.txt
`bash -c "comm -12 <(sort all_genes.list) <(sort genes_with_terms.txt) > all_genes.list.filtered" `

cut -f 1 gene2GOID.tsv | sort | uniq > genes_with_termsID.txt
`bash -c "comm -12 <(sort all_genesID.list) <(sort genes_with_termsID.txt) > all_genes.listID.filtered" `

cut -f 1 gene2GO_genbank.tsv | sort | uniq > genes_with_terms_genbank.txt
`bash -c "comm -12 <(sort all_genes_genbank.list) <(sort genes_with_terms_genbank.txt) > all_genes.list_genbank.filtered" `
