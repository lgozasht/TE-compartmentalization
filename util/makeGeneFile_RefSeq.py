import sys



def grabGenes(gff,geneDic,first,single,geneID):
    stack = set([]) #only take the first occurance of a gene in the case of duplicates

    with open(gff,'r') as f:
        for line in f:
            if '#' not in line and 'pseudo' not in line and 'Name=' in line and 'protein_coding' in line:
                sp = line.strip().split('\t')
                if first == True:
                    if 'gene' == sp[2]:
                         chrm = sp[0]
                         if int(sp[4]) > int(sp[3]):
                             start = sp[3]
       	       	       	     stop = sp[4]
                         else:
                             start = sp[4]
       	       	       	     stop = sp[3]
                         #geneName = sp[-1].split(';')[1].split('=')[-1]

                         if geneID == True:
                             geneName = sp[-1].split(';')[1].split(':')[-1].strip()
                         else:
                             #if sp[-1].split(',')
                             geneName = sp[-1].split('Name=')[1].split(';')[0].strip()


                         if single == False:
                             if geneName not in stack:
                                 geneDic[geneName] = [chrm,start,stop]
                                 stack.add(geneName)
                         else:
                             if geneName not in geneDic:
                                 geneDic[geneName] = [chrm,start,stop]
                             else:
                                 geneDic[geneName].append(chrm)
       	       	       	         geneDic[geneName].append(start)
       	       	       	         geneDic[geneName].append(stop)

                else:
                    if 'gene' == sp[2]:

                         chrm =	sp[0]
       	       	       	 if int(sp[4]) > int(sp[3]):
                             start = sp[3]
       	       	       	     stop = sp[4]
       	       	       	 else:
                             start = sp[4]
                             stop = sp[3]

       	       	       	 #geneName = sp[-1].split(';')[1].split('=')[-1]
                         geneName = sp[-1].split('Name=')[1].split(';')[0].strip()

       	                 if geneName in geneDic and geneName not in stack:      	    
                             geneDic[geneName].append(chrm)
                             geneDic[geneName].append(start)
                             geneDic[geneName].append(stop)
       	       	       	     stack.add(geneName)


    return geneDic

if len(sys.argv) == 4:
    geneDic = grabGenes(sys.argv[1],{},True,False)
    geneDic = grabGenes(sys.argv[2],geneDic,False,False)
    geneDic = grabGenes(sys.argv[3],geneDic,False,False)
    with open('orthologs.tsv','w') as f:
        for gene in geneDic:
            if len(geneDic[gene]) > 8:
                f.write('{0}\t{1}\n'.format(gene,'\t'.join(geneDic[gene])))
            

else:
    geneDic = grabGenes(sys.argv[1],{},True,True,True)
    with open('geneFileID.tsv','w') as f:
        for gene in geneDic:
            if len(geneDic[gene]) == 3:
       	        f.write('{1}\t{0}\t.\t.\n'.format(gene,'\t'.join(geneDic[gene])))
            else:
                for i in range(0,len(geneDic[gene]),3):
                    f.write('{1}\t{0}\t.\t.\n'.format(gene,'\t'.join(geneDic[gene][i:i+3])))

    geneDic = grabGenes(sys.argv[1],{},True,True,False)
    with open('geneFile.tsv','w') as f:
        for gene in geneDic:
            if len(geneDic[gene]) == 3:
                f.write('{1}\t{0}\t.\t.\n'.format(gene,'\t'.join(geneDic[gene])))
            else:
                for i in range(0,len(geneDic[gene]),3):
                    f.write('{1}\t{0}\t.\t.\n'.format(gene,'\t'.join(geneDic[gene][i:i+3])))


