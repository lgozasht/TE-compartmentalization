from sequenceAnalyzer import FastAreader


proteinList = []
with open('proteins.fa.vs.RepeatPeps.25cul2.1e5.blastp.out','r') as f:
    for line in f:
        sp = line.strip().split('\t')
        protein = sp[0]
        proteinList.append(protein)



proteinList = list(set(proteinList))
headerList = []
seqList = []

thisReader = FastAreader('transcripts.fa')
for head, seq in thisReader.readFasta():

    headerList.append(head)
    seqList.append(seq)

removeList = []
for protein in proteinList:
    removeList = removeList + list(filter(lambda x: protein in x, headerList))
        

removeList = set(removeList)
with open('transcripts.no_tes.fa','w') as f:
    for i in range(len(headerList)):
        if headerList[i] in removeList:
            print(headerList[i])
        else:
            f.write('>{0}\n{1}\n'.format(headerList[i], seqList[i]))
