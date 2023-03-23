'''
Created on Aug 10, 2016

@author: ODLab
'''
from varLookup import FOLDER_IN
from MSImageClass import MSImage
from myFunc import sort_nicely
import glob
import os

if __name__ == '__main__':
    control=False
    if control:
        folderName = FOLDER_IN + 'cropped/EMBD0000'
        strains = glob.glob('{0}/*'.format(folderName))
        for strain in strains:
            if os.path.isdir(strain):# and strain=='Z:/cropped/EMBD0000\MS':
                dates = glob.glob('{0}/*'.format(strain))
                for date in dates:
                    embs = glob.glob('{0}/*'.format(date))
                    for emb in embs:
                        try:
                            nZ, dz = 18,1
                            nT, dt = 31,30
                            imNames = glob.glob('{0}/*.tif'.format(emb))
                            if len(imNames)>0:
                                imName = glob.glob('{0}/*.tif'.format(emb))[0][:-14]+'All.tif'
                                if not os.path.exists(imName):
                                    print('converting {0}'.format(emb))
                                    im = MSImage(emb+'/', (nZ, dz), nT, 3)
                                    im.saveIms('{0}'.format(imName))
                        except:
                            print('ERROR, couldnt handle {0}'.format(emb))
        print('ALL DONE')
    else:
        folderName = FOLDER_IN + 'cropped/'
        rnais = glob.glob('{0}/*'.format(folderName))
        sort_nicely(rnais)
        for rnai in rnais:
            if os.path.isdir(rnai):
                strains = glob.glob('{0}/*'.format(rnai))
                for strain in strains:
                    if os.path.isdir(strain):
                        embs = glob.glob('{0}/*'.format(strain))
                        for emb in embs:
                            try:
                                nZ, dz = 18,1
                                nT, dt = 31,30
                                imNames = glob.glob('{0}/*.tif'.format(emb))
                                if len(imNames)>0:
                                    imName = glob.glob('{0}/*.tif'.format(emb))[0][:-14]+'All.tif'
                                    if not os.path.exists(imName):
                                        print('converting {0}'.format(emb))
                                        im = MSImage(emb+'/', (nZ, dz), nT, 3)
                                        im.saveIms('{0}'.format(imName))
                            except:
                                print('ERROR, couldnt handle {0}'.format(emb))
        print('ALL DONE')