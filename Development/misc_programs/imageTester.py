'''
Created on Mar 21, 2017

@author: Admin
'''
import numpy as np
from RNAiClass import RNAiClass
from varLookup import printLog
from skimage.filters import threshold_li
import scipy.ndimage as ndimage
from myFunc import getLargestObj

def getAvgProjInt(emb):
    emb.loadImages()
    allIms = emb.image.images[:,:,:2].astype(np.float32)
    allIms[:,:,0]-=emb.thresh[0]
    allIms[:,:,1]-=emb.thresh[1]
    allIms[np.where(allIms<0)]=0
    ims = np.mean(allIms, axis=(0,1,2))
    return np.max(ims)

def checkEmbDebris(emb):
    if getAvgProjInt(emb)>900:
        printLog('CHECK for debris {d} {r}/{s}/{e}'.format(d=emb.date, r=emb.RNAi, s=emb.strain, e=emb.label))
        return True
#         emb.showMaxProj(7)

def checkRNAi(i):
    r = RNAiClass(i)
    r.refresh_params_data()
    embs = r.GLSMoveEmbs+r.GLSNoMoveEmbs+r.MSMoveEmbs+r.MSNoMoveEmbs
    printLog('checking RNAi={0}'.format(i))
    for emb in embs:
        checkEmbDebris(emb)
        
def checkEmbFocus(emb):
    if emb.image is None: emb.loadImages()
#     if checkEmbDebris(emb):
#         print('Check debris present')
#         return False
    im = emb.image.images[0,9,0]
    th = threshold_li(im)
    radius, length, depth = [], [], []
    for slide in range(10):
        mask = (emb.image.getMaxProject(slide, 0, th)>0)
        label_im = getLargestObj(mask)
        # Now that we have only one connect component, extract it's bounding box
        slice_x, slice_y = ndimage.find_objects(label_im==1)[0]
        radius.append((slice_x.stop - slice_x.start)/2.)
        length.append((slice_y.stop - slice_y.start)/2.)
        
        ims = emb.image.images[slide,:,0]
        ims = ims[:,:,:ims[0].shape[1]/2]
        imProj = np.max(ims, axis=2)
        imProj = (imProj>th)
        label_im = getLargestObj(imProj)
#         fig = myFigure()
#         fig.imshow(zoom(label_im, (emb.dz, 1), order=0))
#         fig.show()
        slice_x, slice_y = ndimage.find_objects(label_im==1)[0]
        depth.append(emb.dz*(slice_x.stop - slice_x.start)/2.)
    return 0.9*np.mean(radius)<=np.mean(depth)

def checkSignal(emb):
    for c in range(2):
        y = emb.tInt[c]
        dy = np.roll(y,-1)-y
        np

if __name__ == '__main__':
    printLog('STARTING DEBRIS TEST')
#     for i in range(1,504):
#         checkRNAi(i)
#     embs = getMSEmbs()
#     intens = np.array([getAvgProjInt(emb) for emb in embs])
#     fig = myFigure()
#     fig.hist(intens)
#     fig.show()
        