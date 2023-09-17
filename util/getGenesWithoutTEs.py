
from sequenceAnalyzer import FastAreader
import sys
import os

#species = sys.argv[1]
thisReader = FastAreader('transcripts.no_tes.fa')
with open('passing_genes.list','w') as f:
    for head, seq in thisReader.readFasta():
        #if 'mrna' in head.lower():
        if "locus_tag=" in head: 
            f.write('{0}\n'.format(head.split('locus_tag=')[-1].split(' ')[0].strip(']')))
        if "GeneID:" in head:

            f.write('{0}\n'.format(head.split('GeneID:')[-1].split(' ')[0].strip(']')))

        if 'gene=' in head:
            f.write('{0}\n'.format(head.split('gene=')[-1].split(' ')[0].strip(']')))

os.system('sort passing_genes.list | uniq > passing_genes.list.sorted.uniq')

passDic = {}
with open('passing_genes.list.sorted.uniq','r') as passFile:
    for line in passFile:
        gene = line.strip()
        passDic[gene] = ''

with open('geneFile_final_TEs.tsv','r') as f:
    with open('geneFile_final.tsv','w') as outFile:
        for line in f:
            if line.strip().split('\t')[3] in passDic:
                outFile.write(line.strip() + '\n')

with open('all_genes.list_final_TEs.filtered','r') as f:
    with open('all_genes.list_final.filtered','w') as outFile:
       	for line in f:
       	    if line.strip() in passDic:
       	       	outFile.write(line.strip() + '\n')

