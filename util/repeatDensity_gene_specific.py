import numpy as np 
import pandas as pd 
import sys 
import os

def getTopPercentile(densityArray, percentile):
    p = np.percentile(densityArray, percentile)
    topPercentileIndexes = [i for i,v in enumerate(densityArray) if v > p]
    return p,topPercentileIndexes

densityDic = {}

currentDir = os.getcwd()
coverageFile = '{0}/{1}'.format(currentDir,sys.argv[1])
TEDensity = pd.read_csv(coverageFile,sep = '\t')
TEDensity.columns = ['chrm','start','stop','density','gene']

densityArray = TEDensity['density']

p,indexList = getTopPercentile(densityArray, int(sys.argv[3]))

with open('{0}/{1}_hotspots_{2}.bed'.format(currentDir,sys.argv[2],str(sys.argv[3])),'w') as f:

    for index, row in TEDensity.iterrows():
        currentDensity = row['density']

        if currentDensity > p:
            f.write("\t".join([str(i) for i in row.tolist()])+"\n")



