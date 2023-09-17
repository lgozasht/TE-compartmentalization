

import sys


first = True 
inputFile = sys.argv[1]
flankCodingBases = {}
genePos = {}

with open('geneFile_final_primary_scaffolds.bed','r') as f:
    for line in f:
        sp = line.strip().split('\t') 
        genePos[sp[3]] = [sp[1],sp[2]]

with open('geneFile_coding_flanks.bed','r') as f:
    for line in f: 
        sp = line.strip().split('\t')

        gene = sp[3].strip()
        if first == True:
            flankCodingBases[gene] = ['','','','']
            leftLength = int(sp[-3])
            leftGene = sp[3].strip()
            flankCodingBases[gene] = [leftLength,'',genePos[gene][0],genePos[gene][1]] 
            first = False
        elif leftGene == gene:
            #print(gene)

            rightLength = int(sp[-3])
            flankCodingBases[gene] = [leftLength,rightLength, genePos[gene][0],genePos[gene][1]] 
            first = True
        else:
            print(leftGene)
            flankCodingBases[gene] = ['','','','']
            leftLength = int(sp[-3])
            leftGene = sp[3].strip()
            flankCodingBases[gene] = [leftLength,'',genePos[gene][0],genePos[gene][1]] 
            first = False




first = True 
with open(inputFile,'r') as f:
    with open(inputFile.replace('.pre',''),'w') as outF:
        for line in f: 
            sp = line.strip().split('\t')
            gene = sp[-1].strip()
            if first == True:
                leftLength = (int(sp[2]) - int(sp[1])) 
                leftGene = sp[-1].strip()
                leftDensity = float(sp[3])
                leftChrom = sp[0]
                first = False
            elif leftGene == gene:
                rightLength = (int(sp[2]) - int(sp[1])) 
                rightDensity =  float(sp[3])
                if (leftLength + rightLength - flankCodingBases[gene][0] - flankCodingBases[gene][1]) == 0:
                    finalDensity = 0
                else:
                    finalDensity = ((leftLength * leftDensity) + (rightLength * rightDensity))/(leftLength + rightLength - flankCodingBases[gene][0] - flankCodingBases[gene][1])
                outF.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(sp[0],flankCodingBases[gene][2],flankCodingBases[gene][3],str(finalDensity),sp[4]))
                first = True
            else:
                print(leftGene)
                if (leftLength  - flankCodingBases[leftGene][0]) == 0:
                    finalDensity = 0
                else:
                    finalDensity = (leftLength * leftDensity)/(leftLength  - flankCodingBases[leftGene][0])
                outF.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(leftChrom,flankCodingBases[leftGene][2],flankCodingBases[leftGene][3],str(finalDensity),leftGene))

                leftLength = (int(sp[2]) - int(sp[1])) 
                leftGene = sp[-1].strip()
                leftDensity = float(sp[3])

                first = False

