'''
Created on Feb 2, 2017

@author: Renat
'''
from testEmbAlignMS import loadTestEmbryosRNAi
from varLookup import AVG_CURVES_MS, nLambda_MS

def initialize():
    embs, data = loadTestEmbryosRNAi()
    print('Lambda = {0}'.format(nLambda_MS))
    for emb, d in zip(embs, data):
        print('{0} {1}'.format(emb.RNAi, emb.label))
        print('scored point of deviation = {0}'.format(d[12]))
        for curveName in AVG_CURVES_MS:  print('point of deviation for {0} = {1}'.format(curveName, emb.cutoff[curveName]))

if __name__ == '__main__':
    initialize()