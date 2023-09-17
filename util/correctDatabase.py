import sys

geneDic = {}
if sys.argv[1] == "ID":
    suffix = "ID"
elif sys.argv[1] == "genbank":
    suffix = "_genbank"
else:
    suffix = ''
with open('all_genes{0}.list'.format(suffix),'r') as f:
    for line in f:
        geneDic[line.strip().upper()] = line.strip()

with open('localUniprot{0}.tsv'.format(suffix),'r') as f:
    with open('gene2GO{0}.tsv'.format(suffix),'w') as finalOut:
        for line in f:
            line = line.strip()
            sp = line.split('\t')
            found = False
            if len(sp) < 2:
                continue
            genes = line.split('\t')[0]
            goTerms = line.split('\t')[1]
            if ' ' in line:
                geneNameList = genes.split()
                for gene in geneNameList:
                    found = False
                    if gene.upper() in geneDic:
                        if suffix == 'ID':
                            if gene.isdigit() == True:
                                finalOut.write('{0}\t{1}\n'.format(geneDic[gene.upper()], goTerms.replace(' ','')))
                        elif suffix == "_genbank":
                            if '.' not in gene.upper():
                                for k in range(10):
                                    if gene.upper() + '.' + str(k) in geneDic:
                                        finalOut.write('{0}\t{1}\n'.format(geneDic[gene.upper() + '.' + str(k)], goTerms.replace(' ','')))
       	       	       	                found = True
                                        break
                                if found == False:

                                    finalOut.write('{0}\t{1}\n'.format(geneDic[gene.strip().upper()], goTerms.replace(' ','')))

                            else:
                                finalOut.write('{0}\t{1}\n'.format(geneDic[gene.upper()], goTerms.replace(' ','')))


                        else:
                            finalOut.write('{0}\t{1}\n'.format(geneDic[gene.upper()], goTerms.replace(' ','')))
            elif genes.strip().upper() in geneDic:
                if suffix == 'ID':
                    if gene.isdigit() == True:
                        finalOut.write('{0}\t{1}\n'.format(geneDic[genes.strip().upper()], goTerms.replace(' ','')))
       	       	elif suffix == "_genbank":
       	       	    if '.' not in genes.strip().upper():

                        for k in range(10):
       	       	       	    if genes.strip().upper() + '.' + str(k) in geneDic:
                                finalOut.write('{0}\t{1}\n'.format(geneDic[genes.strip().upper() + '.' + str(k)], goTerms.replace(' ','')))
                                found = True
                                break
                        if found == False:
                            finalOut.write('{0}\t{1}\n'.format(geneDic[genes.strip().upper()], goTerms.replace(' ','')))

                    else:
                        finalOut.write('{0}\t{1}\n'.format(geneDic[genes.strip().upper()], goTerms.replace(' ','')))


                else:
                    finalOut.write('{0}\t{1}\n'.format(geneDic[genes.strip().upper()], goTerms.replace(' ','')))

