import os
import os.path
import glob
import sys
from sequenceAnalyzer import FastAreader

#check for annotation

def checkAnnotation(link, species):
    end = link.split('/')[-1].strip('\"')
    gffLink = link + '/' +  end + '_genomic.gff.gz'
    gffFile = end + '_genomic.gff.gz'
    os.system('wget {0}'.format(gffLink))
    if 'GCF' in gffFile:
        if os.path.isfile(gffFile):


            os.system('gunzip {0}'.format(gffFile))
            os.system('python ../util/makeGeneFile_RefSeq.py {0}'.format(gffFile.replace('.gz','')))
            os.system('mv geneFile.tsv {0}'.format(species))
            os.system('mv geneFileID.tsv {0}'.format(species))        
            os.system('cut -f 4 {0}/geneFile.tsv > {0}/all_genes.list'.format(species))
            os.system('cut -f 4 {0}/geneFileID.tsv > {0}/all_genesID.list'.format(species))
            os.system('rm -r {0}'.format(gffFile.replace('.gz','')))

        else:
            return False
    else:
        if os.path.isfile(gffFile):

            os.system('gunzip {0}'.format(gffFile))
            os.system('python ../util/makeGeneFile_GenBank.py {0}'.format(gffFile.replace('.gz','')))
            with open("geneFile_genbank.tsv", 'r') as fp:
                lines = len(fp.readlines())
                if lines < 50:
                    return False
            os.system('mv geneFile_genbank.tsv {0}'.format(species))
            os.system('cut -f 4 {0}/geneFile_genbank.tsv > {0}/all_genes_genbank.list'.format(species))
            os.system('rm -r {0}'.format(gffFile.replace('.gz','')))

        else:
            return False

    return True
            

def compareAnnotations(species):

    first = True
    fileList = glob.glob('all_genes.list*.filtered')
    for file in fileList:
        if file == "all_genes.list_final_TEs.filtered" or file == "all_genes.list_final.filtered":
            continue

        print(file)
        with open(file, 'r') as fp:
            lines = len(fp.readlines())
        if first == True:
            bestFile = file
            bestCount = lines
            first = False
        elif lines > bestCount:
       	    bestFile = file
            bestCount =	lines
    print("The best dataset for {0} is {1} with {2} functionally annotated genes".format(species, bestFile, bestCount)) 
    if 'genbank' in bestFile:
        database = 'genbank'
        os.system('cp genes_with_terms_genbank.txt genes_with_terms_final.txt')
        os.system('cp geneFile_genbank.tsv geneFile_final_TEs.tsv')
        os.system('cp all_genes.list_genbank.filtered all_genes.list_final_TEs.filtered')
        os.system('cp gene2GO_genbank.tsv gene2GO_final.tsv')

    elif 'ID' in bestFile:
        database = "ID"
        os.system('cp genes_with_termsID.txt genes_with_terms_final.txt')
        os.system('cp geneFileID.tsv geneFile_final_TEs.tsv')
        os.system('cp all_genes.listID.filtered all_genes.list_final_TEs.filtered')
        os.system('cp gene2GOID.tsv gene2GO_final.tsv')

    else:
        database = "refseq"
        os.system('cp genes_with_terms.txt genes_with_terms_final.txt')
        os.system('cp geneFile.tsv geneFile_final_TEs.tsv')
        os.system('cp all_genes.list.filtered all_genes.list_final_TEs.filtered')
        os.system('cp gene2GO.tsv gene2GO_final.tsv')
         
    
    return bestFile, database


def downloadData(link, species, refseq):

    end = link.split('/')[-1].strip('\"')
    fnaLink = link + '/' +  end + '_genomic.fna.gz'
    fnaFile = end + '_genomic.fna.gz'
    transcriptLink = link + '/' +  end + '_cds_from_genomic.fna.gz' #'_rna_from_genomic.fna.gz'
    transcriptFile = end + '_cds_from_genomic.fna.gz' #'_rna_from_genomic.fna.gz'

    proteinLink = link + '/' +  end + '_protein.faa.gz'
    proteinFile = end + '_protein.faa.gz'
    if os.path.isfile(fnaFile.replace('.gz','')):
        pass
    else:
     
        os.system('wget {0}'.format(fnaLink))
        os.system('gunzip {0}'.format(fnaFile))

    os.system('rm -r {0}'.format(transcriptFile))

    os.system('wget {0}'.format(transcriptLink))
    os.system('gunzip {0}'.format(transcriptFile))
    thisReader = FastAreader(transcriptFile.replace('.gz',''))
    with open('transcripts.fa','w') as f:
    
        for head, seq in thisReader.readFasta():
            f.write('>{0}\n{1}\n'.format(head,seq)) 


    os.system('rm -r {0}'.format(transcriptFile.replace('.gz','')))

    os.system('wget {0}'.format(proteinLink))
    os.system('gunzip {0}'.format(proteinFile))
    os.system('mv {0} proteins.fa'.format(proteinFile.replace('.gz','')))


    os.chdir('..')
    os.system('sh ../util/convertFasta.sh {0} {1}'.format(species,fnaFile.replace('.gz','')))
    os.chdir(species)
    
    return fnaFile

def runAnalysis(fnaFile, species, TEProteinFile, goOboFile, cutoff):
    os.system('cp ../TE-compartmentalization_pipeline.sh  .')
    os.system('sh TE-compartmentalization_pipeline.sh  {0} {1} {2} {3} {4}'.format(species,fnaFile.replace('.gz',''), TEProteinFile, goOboFile, cutoff))
    os.chdir('..')



linkFile = sys.argv[1] 
TEProteinFile = sys.argv[2]
goOboFile = sys.argv[3]
cutoff = sys.argv[4]


with open('{0}'.format(linkFile),'r') as f:
    for line in f:
        if 'GCF' not in line:
            continue
        species = line.strip().split(',')[0].strip('\"').strip('\"').replace(' ','_').replace('/','_').replace('(','_').replace(')','_').replace('[','').replace(']','').replace('\]','').replace('\[','')
        os.system('mkdir {0}'.format(species))

        
        refSeqLink = line.strip().split(',')[-1].strip('\"').strip('\"')
        genBankLink = line.strip().split(',')[-2].strip('\"').strip('\"')
        refComplete = checkAnnotation(refSeqLink, species)
       	genComplete = checkAnnotation(genBankLink, species)

        ####Landen was here
        os.chdir(species)
        os.system('sh ../util/make_local_uniprot.sh {0}'.format(species))

        bestFile, dataType = compareAnnotations(species)
        
        if 'genbank' == dataType:
            fnaFile = downloadData(genBankLink, species, refSeqLink)
        else:
            fnaFile= downloadData(refSeqLink, species, refSeqLink)
    
        
        runAnalysis(fnaFile, species)    
        
