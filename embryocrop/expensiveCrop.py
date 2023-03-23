'''
Created on Nov 9, 2015

@author: renat

performs time consuming crop
'''

date = '20140515T131046' 

loadFolder = 'Z:/'
z = 18 #number of z planes
myWells=[7]
myPVs = [0]
nWells = 1#number of wells (14)
pointVisits = 1# number of point visits (4)
RNAi, strains = [], [] #RNA condition, strain condition
trackingFile = 'Z:/Experiment_tracking_sheets/EMBD_fileNames_Tracking_Sheet.csv'

import cv2, glob, Image, time, csv, shutil
import numpy as np
import Tkinter
import tkMessageBox
from findEmbryo import *
from myFunc import *
import multiprocessing as mp
from itertools import repeat
from AttenuationCorrection import correctAttAll
from particleTracking import findTracks

debug=False
splitSym = '/'
if debug: nWells,pointVisits=1,1
folderIn = loadFolder + 'CV1000/' + date #Input folder

def removeBG(im, well):
    im = np.float32(im) - 3000
    im[np.where(im<0)]=0
    im = np.uint16(im)
    if strains[well]!='MS':
        Gsize = 41
    else:Gsize = 201
    karnel = numpy.ones([Gsize,Gsize])
    karnel = 1.0*karnel/sum(karnel)
    imBG = cv2.GaussianBlur(im, (Gsize,Gsize), 0)
    imTmp = numpy.float32(im)-numpy.float32(imBG)
    imTmp[numpy.where(imTmp<0)]=0
    if numpy.max(imTmp)==0:
        print('zero image')
        imTmp= numpy.uint16(imTmp)
    elif np.max(imTmp)>65535:
        print('large value image')
        imTmp= numpy.uint16(imTmp)
        showIm(im, 'original')
        showIm(imTmp, 'subtracted')
    imTmp= numpy.uint16(imTmp)
    return imTmp

def getConditions(date, fileName):
    ''' loads RNAi strains for a specified date from a csv file '''
    global RNAi, strains
#     csvFile = csv.reader(open(fileName, 'rb'), delimiter=',')
    csvFile = csv.reader(open(fileName, 'rU'), delimiter=',') #universal
    fileData=[]
    for row in csvFile:
        fileData.append(row[1:-1])
    myDate = [s for s in fileData if s[0]==date]
    myDate = sorted(myDate, key=lambda well: well[2])
    RNAi = [s[3] for s in myDate]
    strains =  [s[4] for s in myDate]
    return

def loadImages(folder, well, j):
    '''
    loads images from a folder and splits them in separate point visits
    
    Parameters:
    folder : folder to read images from
    
    Return:
    allImgs: list of 4 different point visits with 3 channels in each. images are numpy arrays.
    '''
    imc1, imc2, imc3 = [], [], []
    folderNames = glob.glob(folder+'/Well{0:0>3}/*F{1:0>3}*C1.tif'.format(well,j))
    folderNames.sort()
    for fileName in folderNames:
        imc1.append(cv2.imread(fileName, -1))
        
    folderNames = glob.glob(folder+'/Well{0:0>3}/*F{1:0>3}*C2.tif'.format(well,j))
    folderNames.sort()
    for fileName in folderNames:
        imc2.append(cv2.imread(fileName, -1))
        
    folderNames = glob.glob(folder+'/Well{0:0>3}/*F{1:0>3}*C3.tif'.format(well,j))
    folderNames.sort()
    for fileName in folderNames:
        imc3.append(cv2.imread(fileName, -1))
        
    allImgs=(imc1, imc2, imc3)
    return allImgs

def getPlane(imgs, z, j):
    '''
    outputs image list of only specified plane
    '''
    res = [imgs[i*z+j] for i in range(len(imgs)/z)]
    return res

def correctDrift4AllC(imgs):
    '''
    Correct drift for all images based on central plane of the third image sequence
    
    Parameters:
    imgs : tuple of three lists of images
    
    Return:
    tuple of three lists of corrected images
    '''
    im1, im2, im3 = imgs
    centralPlane = getPlane(im3, z, z/2)
    drift = findDrift(centralPlane)
    drift4All = []
    for d in drift:
        for i in range(z):
            drift4All.append(d)
    im1c = correctDrift(im1,drift4All)
    im2c = correctDrift(im2,drift4All)
    im3c = correctDrift(im3,drift4All)
    corrected = [im1c, im2c, im3c]
    print('drift corrected')
    return corrected

def cropAllC(imgs, well):
    '''
    Crops all images based on the embryos found in the first central plane of C3
    
    Parameters:
    imgs : tuple of three lists of images of different channels
    
    Return:
    list of tuples of three channels for each embryo
    '''
    im1, im2, im3 = imgs
    if True:
        imtmp=[]
        for im in im1:
            imtmp.append(removeBG(im, well))
        im1=imtmp
        imtmp=[]
        for im in im2:
            imtmp.append(removeBG(im, well))
        im2=imtmp
        del imgs, imtmp
    allParams = []
    maskPrev=None
    for tCrop in range(len(im1)/z):
        print('well {0}, scanning time point={1}'.format(well+1, tCrop+1))
        ims8b = np.uint8(255.*(im3[tCrop*z:(tCrop+1)*z] - np.min(im3[tCrop*z:(tCrop+1)*z]))/np.max(im3[tCrop*z:(tCrop+1)*z]))
        mask = getMaskStak(ims8b)
        if maskPrev is not None and np.sum(mask)<0.5*np.sum(maskPrev):
            mask=maskPrev
            print('replace mask at t={0}'.format(tCrop))
        if np.sum(mask)/255>0.5*mask.size:
            mask = getMask(a16a8(im3[tCrop*z+z/2]), True)
            print('cropAllC, too much dirt')
        eParams = findEmbsonIm(mask)
        if len(eParams)>0: allParams.append(eParams)
        else: allParams.append(allParams[-1])
        maskPrev = mask
    allParams=pairEmbs(allParams) 
    j=0
    allEmb, aspRatio = [], []
    pool = mp.Pool(processes=10)
    for params in allParams:
        print('well {0}, cropping Embryo={1}'.format(well+1, j+1))
        apR = []
        for k in [6,8,10]:
            imtmp = cropRotate([maxIntensProject(correctAttAll(im2[k*z:(k+1)*z],z,0,0.1)),params[k], False])
            apR.append(getAP(imtmp))
        if numpy.mean(apR)<0.8:flip = True
        elif numpy.mean(apR)>1.25: flip = False
        elif strains[well]=='GLS':
            apR = []
            for k in [8,10,12]:
                imtmp1 = cropRotate([maxIntensProject(correctAttAll(im1[k*z:(k+1)*z],z,0,0.1)),params[k], False])
                imtmp2 = cropRotate([maxIntensProject(correctAttAll(im2[k*z:(k+1)*z],z,0,0.1)),params[k], False])
                apR.append(getAP2(imtmp1, imtmp2))
        if numpy.mean(apR)<1: flip = True
        else: flip = False
        pAll = []
        for k in range(len(params)):
            for i in range(z):
                pAll.append(params[k])
        imAll1 = pool.map(cropRotate, zip(im1, pAll, repeat(flip)))
        imAll2 = pool.map(cropRotate, zip(im2, pAll, repeat(flip)))
        imAll3 = pool.map(cropRotate, zip(im3, pAll, repeat(flip)))
        emb = (imAll1,imAll2,imAll3)
        allEmb.append(emb)
        (a,b), center, angle = params[0]
        aspRatio.append(1.*a/b)
        j+=1
        del emb, imAll1, imAll2, imAll3
    del im1, im2, im3
    pool.close()
    pool.join()
    return allEmb, aspRatio

def pairEmbs(allParams):
    pointsList = [[p for tmp, p, a in params] for params in allParams]
    trackList = findTracks(pointsList, 50)
    aMax, bMax = np.zeros((2,len(trackList)))
    for j in range(len(trackList)):
        tr = trackList[j]
        time, points, indx = tr.getTrack()
        aMax[j]=np.mean([allParams[time[i]][indx[i]][0][0] for i in range(tr.getLength())])
        bMax[j]=np.mean([allParams[time[i]][indx[i]][0][1] for i in range(tr.getLength())])
    paramsGrouped = []
    for j in range(len(trackList)):
        tr = trackList[j]
        print('track {0} length {1}'.format(j,tr.getLength()))
        time, points, indx = tr.getTrack()
        if time[0]==0:
            params = [((aMax[j], bMax[j]),allParams[time[0]][indx[0]][1],allParams[time[0]][indx[0]][2])]
            for i in range(1,tr.getLength()):
                ang1 = params[-1][-1]
                ang2 = allParams[time[i]][indx[i]][2]
                if abs(ang1-ang2)>abs(ang1-(np.pi+ang2)): ang2=np.pi+ang2
                elif abs(ang1-ang2)>abs(ang1-(-np.pi+ang2)): ang2=-np.pi+ang2
                params.append(((aMax[j], bMax[j]),allParams[time[i]][indx[i]][1],ang2))
            if len(params)<len(allParams):
                for i in range(len(allParams)-len(params)):
                    params.append(params[-1])
            paramsGrouped.append(params)
    return paramsGrouped

def getAllEmb(folder):
    '''
    Finds and saves all embryos
    
    ParametersL
    folder: folder to load embryos
    
    Return:
    None
    '''
    global RNAi, strains
    
    print('STARTED!!!')
    getConditions(date, trackingFile)
    totalEmb = 0
    embs=[]
    for well in myWells:
#         try:
            for j in myPVs:
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#
                imgs = loadImages(folder, well+1, j+1)
                if len(imgs[0])>0:
                    print('loaded well {0} point {1}'.format(well+1, j+1))
                    es, rs = cropAllC(imgs, well)
                    embs.append([well, j, es, rs])
                    totalEmb+=len(embs[-1])-2
                else: print('no images well {0} point {1}'.format(well+1, j+1))
                del imgs
#         except Exception,e: print('Error: '+str(e))
        
    print('Done analyzing, ready to save!')
    checks = np.array([checkEmbDebris(e[2][3*z+z/2]) for well, j, es, rs in embs for e in es])
    totalEmb = checks.size
    uniqeRNAi = np.array(list(set(RNAi)))
    embN = np.ones([len(uniqeRNAi),2])
    xembN = np.ones([len(uniqeRNAi),2])
    if not debug:
        for well,j, es, rs in embs:
            if strains[well]!='MS':k=0
            else: k=1
            for l in range(len(es)):
                e=es[l]
                r=rs[l]
                print('{0} embryos left, saving...'.format(totalEmb))
                i = checks.size-totalEmb
                if checks[i]==1:
                    saveEmb(e,j,int(embN[np.where(uniqeRNAi==RNAi[well])[0],k]), well, checks[i],r)
                    embN[np.where(uniqeRNAi==RNAi[well])[0],k]+=1
                elif checks[i]==2:
                    saveEmb(e,j,int(xembN[np.where(uniqeRNAi==RNAi[well])[0],k]), well, checks[i],r)
                    xembN[np.where(uniqeRNAi==RNAi[well])[0],k]+=1
                totalEmb-=1

def checkEmbDebris(im):
    ''' lets user debug supplied image, and returns 1 for yes (save) and 0 for not and 2 for special case'''
    code = showIm(im)
    if code == ord('d') or code == ord('D'):
        result = tkMessageBox.askquestion("Delete", "Are You Sure?", icon='warning')
        if result == 'yes':
            return 0
        else:
            return checkEmbDebris(im)
    elif code == ord('x') or code == ord('X'): return 2
    else: return 1

def saveEmb(imgs, point, i, well, check, r):
    '''
    Saves embryo images according to a certain pattern
    
    Parameters:
    imgs: tuple of 3 lists of images for each channel
    f: point visit number
    i: embryo number
    r: aspect ratio
    '''
    
    im1, im2, im3 = imgs
    strain = strains[well]
    ri = RNAi[well]
    if ri!='EMBD0000': folderOut = loadFolder + 'tmp/cropped/{0}/{1}/'.format(ri,strain) #outputFolder
    else: folderOut = loadFolder + 'tmp/cropped/EMBD0000/{0}/{1}/'.format(strain,date) #outputFolder
    if check == 2 :folderOut = folderOut+'x'
    else: folderOut = folderOut
    j = 1
    folder = folderOut+ 'Emb{0}/'.format(j)
    while os.path.exists(folder):
        fileName = glob.glob(folder+'*_T01_Z01_C1.tif')
        if len(fileName)>0:
            fileName = fileName[0].split('/')[-1]
            if fileName.split('_')[2]==date:
                if i==1: break
                else: i -= 1
            j += 1
            folder = folderOut+ 'Emb{0}/'.format(j)
        else:
            j += 1
            folder = folderOut+ 'Emb{0}/'.format(j)
    fileName = '{0}_Emb{1}_{2}_W{3:0>2}F{4}_'.format(ri,j,date,well+1,point+1)
    if not debug:
        ''' correct attenuation and save local '''
        if not os.path.exists(folder): os.makedirs(folder)
        else:
            print('file exist, clearing folder', folder)
            clearFolder(folder)
            if os.path.exists(folder+'batch-output'.format(splitSym)):
                shutil.rmtree(folder+'batch-output'.format(splitSym))
        saveImgs(correctAttAll(im1,z,0,0.1), folder, fileName,1)
        saveImgs(correctAttAll(im2,z,0,0.1), folder, fileName,2)
        saveImgs(im3, folder, fileName,3)
  
        ''' populate aspect ratio file '''
        print('saveEmb aspect j=',j)
        addAspect(r,date, ri, j)
    
def saveImgs(imgs,folder, fileName, c):
    for i in range(len(imgs)):
        cv2.imwrite(folder+fileName+'T{0:0>2}_Z{1:0>2}_C{2}.tif'.format( i/z+1, i%z+1, c), imgs[i])

def addAspect(r, date, ri, j):
    '''
    Adds aspect ratio of an embryo into the csv file
    
    Parameters:
    r: aspect ratio
    date: date
    ri: RNAi conditions
    j: embryo number
    '''
    
    import operator
    fileName = loadFolder+'cropped/'+'aspects.csv'
    newData = []
    oldData = loadAspects(fileName)
    i=0
    while i<len(oldData):
        if oldData[i][0]!=ri: newData.append(oldData[i])
        elif oldData[i][1]!=date: newData.append(oldData[i])
        elif int(oldData[i][2])<j: newData.append(oldData[i])
        i+=1
    newData.append([ri, date, '{0:0>3}'.format(j), str(r)])
    newData = sorted(newData, key=operator.itemgetter(1, 2, 3))
    saveAspects(fileName, newData)

def loadAspects(fileName):
    '''
    reads aspect ratios from file. The aspect ratios are sorted by rnai condition, date, embryo number.
    fileName: name of the file to read from
    '''
    fileData = []
    try:
        csvFile = csv.reader(open(fileName, 'rU'), delimiter=',') #universal
        for row in csvFile:
            fileData.append(row)
    except: pass
    return fileData

def saveAspects(fileName, data):
    with open(fileName, 'wb') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerows(data)

if __name__ == '__main__':
    getAllEmb(folderIn)
    print('ALL DONE!!!!! :)')