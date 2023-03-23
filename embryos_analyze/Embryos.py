"""
Created on Feb 24, 2015

@author: renat

Embryo class
"""
from myMath import getSubIndex, fitSigmBump, getCenterOfMass, calcMoment, getSig
from myFigure import myFigure
from myFunc import getLargestObj, sort_nicely, maxIntensProject
import numpy as np
import csv
import glob
import os
from varLookup import printLog, PARAM_NAMES, T0_2_MOVE_AVG_GLS, T0_2_MOVE_AVG_MS, RADIUS_AVG, LENGTH_AVG, \
    AVG_CURVES_GLS, AVG_CURVES_MS, nLambda_GLS, nLambda_MS, tauScale_MS, tauScale_GLS, FILENAME_AVG_CURVES_MS_ADJTH, \
    FILENAME_AVG_CURVES_GLS, AVG_G_BASE_MS, AVG_R_BASE_MS, FOLDER_IN, GLS_MOVE_PRECISION, AG_AVG_MS, AR_AVG_MS
from MSImageClass import MSImage
import moveDetect
from skimage.filters import threshold_yen, threshold_li
import pickle
import timeAlign
import SkeletonAnalysis
import emb_handler
import db_utils_embryos
from scipy import ndimage
from copy import copy
import cv2
from lmfit import Model, report_fit, Parameters
import matplotlib.pyplot as plt


class Embryo(object):
    """
    General class that describes functions applicable to both Morph and germ strains
    """

    def __init__(self, folder=None, verb=True, check_version=True): #switch back check_version to True! 4/15/22
        """
        Initialization embryo object
        :param folder: embryo folder with images
        :param verb: print or not things along the way
        :param check_version: enforce checking that version stored is the same as the run, otherwise update embryo
        """
        folder = '/'.join(folder.split('\\'))
        if folder[-1] != '/':
            folder += '/'
        self.folder = folder
        self.fileName = folder + 'emb.pkl'
        self.verb = verb
        self.check_version = check_version
        np.random.seed(1)
        self.iniCommonVars()
        self.iniIndividualVars()
        self.iniCalulatedVars()
        if self.folder is not None and self.verb:
            printLog('loading {0}'.format(self.label))
        conn, curs = db_utils_embryos.initialize_db()
        self.id = db_utils_embryos.get_embryo_id(int(self.RNAi[4:]), self.strain, int(self.label[3:]), self.date, curs)
        if self.id is None:
            db_utils_embryos.add_embryo(int(self.RNAi[4:]), self.strain, int(self.label[3:]), self.date, self.folder,
                                        self.field, curs)
            conn.commit()
            self.id = db_utils_embryos.get_embryo_id(int(self.RNAi[4:]), self.strain, int(self.label[3:]), self.date,
                                                     curs)
        self.load()

    def iniCommonVars(self):
        """
        Initialization of variables common to all embryos (MS & GLS)
        :return:
        """
        self.splitSym = '/'
        self.z = 18
        self.channels = 3
        self.timePoints = 31
        self.dz = 8
        self.dt = 1
        self.adjThresh = True
        self.time = np.arange(1. * self.timePoints)

    def iniIndividualVars(self):
        """
        Initialization of variables that change between individual embryos
        :return:
        """
        self.setLables()  # sets self.RNAi, self.label, self.date, self.field
        self.stopInt = 31  # modified in self.fitSigm2Int()
        self.image = None  # modified in self.loadImages()
        self.loaded = True  # modified in self.loadImages()
        self.radius = None  # modified in setSize()
        self.length = None  # modified in setSize()
        self.__version__ = emb_handler.get_embryo_revision_version()

    def iniCalulatedVars(self):
        """
        Initialization of variables that are calculated for individual embryos
        :return:
        """
        self.thresh = [0, 0]  # modified in self.setThresh()
        self.t0 = 0  # modified in timeAlign.***
        self.tScale = 1.  # time scaling
        self.tMove = None  # modified in moveDetect.***
        self.movement = False  # modified in moveDetect.***
        self.dead = False  # modified in moveDetect.***
        self.startPoint = {}  # modified in self.setCurveStartPoints()
        self.scale = {}  # modified in self.setScaleCutoff()
        self.cutoff = {}  # modified in self.setScaleCutoff()
        self.tInt = np.zeros([self.channels, self.timePoints])  # modified in self.setInt()
        self.tIntTh = np.zeros([self.channels, self.timePoints])  # modified in self.setIntTh()
        self.SK = np.zeros([self.channels, self.timePoints])  # modified in self.setIntTh()
        self.dims = np.zeros([self.channels, self.timePoints, 3])  # modified in self.setIntTh()
        self.headInt = np.zeros(self.timePoints)  # modified in self.setIntTh()
        self.CMPos = np.zeros([self.channels, self.timePoints, 3])  # modified in self.getDistFromGreen()
        self.MoI = np.zeros([self.channels, 3, self.timePoints])  # modified in self.getAllMoments()
        self.params = {}  # modified in self.updateParams()
        self.embCenter = np.ones(3) * np.nan  # modified in self.getEmbCenter()
        for pN in PARAM_NAMES:
            self.params[pN] = np.nan

    def refresh(self):
        """
        Refresh of embryo: remeasures all embryo variables
        :return:
        """
        printLog('{d}/{l} refreshing...'.format(d=self.date, l=self.label))
        self.__version__ = emb_handler.get_embryo_revision_version()
        self.iniCommonVars()
        self.iniIndividualVars()
        self.iniCalulatedVars()
        self.loadImages()
        self.setInt()
        self.dead = moveDetect.getDead(self)

    def save(self):
        """
        Save all embryo calculated parameters in MySQL database
        :return:
        """
        conn, curs = db_utils_embryos.initialize_db()
        columns = copy(self.params)
        columns['version'] = self.__version__
        columns['dead'] = self.dead
        columns['t0'] = self.t0
        columns['radius'] = self.radius
        columns['length'] = self.length
        columns['emb_center_x'] = self.embCenter[0]
        columns['emb_center_y'] = self.embCenter[1]
        columns['emb_center_z'] = self.embCenter[2]
        columns['thresh_g'] = self.thresh[0]
        columns['thresh_r'] = self.thresh[1]
        columns['tInt'] = self.tInt
        columns['tIntTh'] = self.tIntTh
        columns['headInt'] = self.headInt
        columns['dims_g_x'] = self.dims[0, :, 0]
        columns['dims_g_y'] = self.dims[0, :, 1]
        columns['dims_g_z'] = self.dims[0, :, 2]
        columns['dims_r_x'] = self.dims[1, :, 0]
        columns['dims_r_y'] = self.dims[1, :, 1]
        columns['dims_r_z'] = self.dims[1, :, 2]
        columns['dims_y_x'] = self.dims[2, :, 0]
        columns['dims_y_y'] = self.dims[2, :, 1]
        columns['dims_y_z'] = self.dims[2, :, 2]
        columns['cmpos_g_x'] = self.CMPos[0, :, 0]
        columns['cmpos_g_y'] = self.CMPos[0, :, 1]
        columns['cmpos_g_z'] = self.CMPos[0, :, 2]
        columns['cmpos_r_x'] = self.CMPos[1, :, 0]
        columns['cmpos_r_y'] = self.CMPos[1, :, 1]
        columns['cmpos_r_z'] = self.CMPos[1, :, 2]
        columns['cmpos_y_x'] = self.CMPos[2, :, 0]
        columns['cmpos_y_y'] = self.CMPos[2, :, 1]
        columns['cmpos_y_z'] = self.CMPos[2, :, 2]
        columns['moi_g_x'] = self.MoI[0, 0]
        columns['moi_g_y'] = self.MoI[0, 1]
        columns['moi_g_z'] = self.MoI[0, 2]
        columns['moi_r_x'] = self.MoI[1, 0]
        columns['moi_r_y'] = self.MoI[1, 1]
        columns['moi_r_z'] = self.MoI[1, 2]
        columns['moi_y_x'] = self.MoI[2, 0]
        columns['moi_y_y'] = self.MoI[2, 1]
        columns['moi_y_z'] = self.MoI[2, 2]
        columns['curve_start_point'] = self.startPoint
        columns['curve_cutoff'] = self.cutoff
        columns['curve_scale'] = self.scale
        columns['skeleton_r'] = self.SK
        db_utils_embryos.update_row(self.id, columns, curs, 'embryos')
        conn.commit()
        conn.close()

    def load(self):
        """
        Load embryo parameters from the database to avoid calculation. If version in database different from
        the run version then refresh the embryo
        :return:
        """
        if self.verb:
            printLog('loading from MySQL...')
        conn, curs = db_utils_embryos.initialize_db()
        row = db_utils_embryos.get_row(self.id, curs, 'embryos')
        ver = row['version']
        if not isinstance(ver, str):  # check that version is not nan otherwise load from pickle
            self.load_from_pickle()
            self.save()
        # elif self.check_version and self.__version__ != ver:  # check version, if different from run refresh
        #     printLog('old version, updating embryo')
        #     self.refresh() #Toggle these 3 lines back on when done fixing head params#
            # print("SWITCH BACK EMBS 202-204")
        else:
            print('using mysql data. Version in MySQL = {}, current version = {}'.format(ver, self.__version__))
            self.__version__ = ver
            self.params = {}
            for par in PARAM_NAMES:
                self.params[par] = row[par]
            self.dead = row['dead']
            self.t0 = row['t0']
            self.tScale = row['tScale']
            self.radius = row['radius']
            self.length = row['length']
            self.embCenter[0] = row['emb_center_x']
            self.embCenter[1] = row['emb_center_y']
            self.embCenter[2] = row['emb_center_z']
            self.thresh[0] = row['thresh_g']
            self.thresh[1] = row['thresh_r']
            self.tInt = row['tInt']
            self.tIntTh = row['tIntTh']
            self.headInt = row['headInt']
            self.dims[0, :, 0] = row['dims_g_x']
            self.dims[0, :, 1] = row['dims_g_y']
            self.dims[0, :, 2] = row['dims_g_z']
            self.dims[1, :, 0] = row['dims_r_x']
            self.dims[1, :, 1] = row['dims_r_y']
            self.dims[1, :, 2] = row['dims_r_z']
            self.dims[2, :, 0] = row['dims_y_x']
            self.dims[2, :, 1] = row['dims_y_y']
            self.dims[2, :, 2] = row['dims_y_z']
            self.CMPos[0, :, 0] = row['cmpos_g_x']
            self.CMPos[0, :, 1] = row['cmpos_g_y']
            self.CMPos[0, :, 2] = row['cmpos_g_z']
            self.CMPos[1, :, 0] = row['cmpos_r_x']
            self.CMPos[1, :, 1] = row['cmpos_r_y']
            self.CMPos[1, :, 2] = row['cmpos_r_z']
            self.CMPos[2, :, 0] = row['cmpos_y_x']
            self.CMPos[2, :, 1] = row['cmpos_y_y']
            self.CMPos[2, :, 2] = row['cmpos_y_z']
            self.MoI[0, 0] = row['moi_g_x']
            self.MoI[0, 1] = row['moi_g_y']
            self.MoI[0, 2] = row['moi_g_z']
            self.MoI[1, 0] = row['moi_r_x']
            self.MoI[1, 1] = row['moi_r_y']
            self.MoI[1, 2] = row['moi_r_z']
            self.MoI[2, 0] = row['moi_y_x']
            self.MoI[2, 1] = row['moi_y_y']
            self.MoI[2, 2] = row['moi_y_z']
            self.startPoint = row['curve_start_point']
            self.cutoff = row['curve_cutoff']
            self.scale = row['curve_scale']
            self.SK = row['skeleton_r']
            self.tScale = self.params['tScale']
            self.movement = self.params['movement']
            if np.isnan(self.params['tMove']):
                self.tMove = np.nan
            else:
                self.tMove = np.round(self.params['tMove'] / self.tScale + self.t0).astype(int)

    def load_from_pickle(self):
        """
        (DEPRECIATED) If embryo is not in the database, load it from pickle.
        :return:
        """
        if not os.path.exists(self.fileName):
            printLog('no pickle found {0}...'.format(self.fileName))
            self.refresh()
        else:
            try:
                printLog('loading from {0}...'.format(self.fileName))
                with open(self.fileName, 'rb') as input:
                    self.__version__ = pickle.load(input)
                    self.thresh = pickle.load(input)
                    self.t0 = pickle.load(input)
                    self.tScale = pickle.load(input)
                    self.tMove = pickle.load(input)
                    self.movement = pickle.load(input)
                    self.dead = pickle.load(input)
                    self.startPoint = pickle.load(input)
                    self.scale = pickle.load(input)
                    self.cutoff = pickle.load(input)
                    self.tInt = pickle.load(input)
                    self.tIntTh = pickle.load(input)
                    self.SK = pickle.load(input)
                    self.dims = pickle.load(input)
                    self.headInt = pickle.load(input)
                    self.CMPos = pickle.load(input)
                    self.MoI = pickle.load(input)
                    self.params = pickle.load(input)
                    self.embCenter = pickle.load(input)
                    if self.embCenter is None:
                        self.embCenter = np.ones(3) * np.nan
                    self.radius, self.length = pickle.load(input)
                    for key in self.startPoint:
                        self.startPoint[key] = int(self.startPoint[key])
                    for pN in PARAM_NAMES:
                        if pN not in self.params: self.params[pN] = np.nan
            except Exception as e:
                printLog('Failed to load from pickle: {0} {1}'.format(self.fileName, str(e)))
                self.refresh()

    def loadImages(self):
        """
        Load embryo images
        self.image[time,z,ch,im.shape]
        """
        printLog('loading images for {0}...'.format(self.label))
        if self.folder is not None:
            fileName = glob.glob('{0}/*_All.tif'.format(self.folder))
            if len(fileName) > 0:
                try:
                    self.image = MSImage(fileName[0], (self.z, self.dz), self.timePoints, self.channels)
                    if np.max(self.image.images) == 0:
                        raise ValueError
                except:
                    print('issue with multi-tif, using singles')
                    os.unlink(fileName[0])
                    self.image = MSImage(self.folder, (self.z, self.dz), self.timePoints, self.channels)
                    imName = glob.glob('{0}/*.tif'.format(self.folder))[0][:-14] + 'All.tif'
                    self.image.saveIms(imName)
            else:
                self.image = MSImage(self.folder, (self.z, self.dz), self.timePoints, self.channels)
                imName = glob.glob('{0}/*.tif'.format(self.folder))[0][:-14] + 'All.tif'
                self.image.saveIms(imName)
                #         self.image.setIms2Edges()
                #         self.edge = copy.copy(self.image)
                #         self.edge.setIms2Edges()
        if self.strain == 'MS':
            print('cleaning image')
            self.imageClean()
        printLog('images loaded')
        self.loaded = False

    def imageClean(self):
        """
        cleans the image from debris or adjacent embryos by finding the largest object in the image
        :return:
        """
        im = self.image.images[0, 9, 0]
        th = threshold_li(im)
        im = self.image.getMaxProject(0, 0, 0)
        mask = (im > th)
        label_im = getLargestObj(mask)
        self.image.images *= label_im

    def setInt(self):
        """
        sets total intensity and its parameters
        :return:
        """
        for c in range(2):
            self.tInt[c, :] = self.image.getAllIntens(c)
        self.params['maxG'] = np.max(self.tInt[0][10:])
        self.params['maxR'] = np.max(self.tInt[1][10:])

    def setScaleCutoff(self):
        """
        Sets the cutoffs (time points of deviation) for each curve. Note: should be done after time alignment.
        :return:
        """
        if self.curvesLoaded:
            # for curveName in list(self.avgCurves.keys()):
            for curveName in self.getCurveNames():
                if np.sum(np.isnan(self.getCurve(curveName)[1])) < 28:
                    x0, nDev, yScale, xScale = timeAlign.getAlignParams(self, [curveName], show=False, shiftVary=False,
                                                                        scaleVary=True, verb=False)
                    #  determine the cutoff by using the time alignment function and allowing it to fit
                    #  only y scale and time point of deviation (no shifting and scaling in time, since already aligned)
                    self.scale[curveName] *= yScale[0]
                    self.cutoff[curveName] = int(nDev[0])
                else:
                    self.scale[curveName] = 1
                    self.cutoff[curveName] = 0

    def setCM_MoIParams(self, cn, strain):
        """
        sets parameters calculated for the center of mass and moment of inertia
        :param cn: curve name
        :param strain: strain name (GLS/MS)
        :return:
        """
        #         if not self.movement:
        x, y = self.getCurve(cn)
        if np.sum(1 - np.isnan(y)) > 2:
            if strain == 'GLS':
                self.params['{0}td_{s}'.format(cn, s=strain)] = min(x[self.cutoff[cn] - 1], T0_2_MOVE_AVG_GLS)
            else:
                self.params['{0}td_{s}'.format(cn, s=strain)] = min(x[self.cutoff[cn] - 1], T0_2_MOVE_AVG_MS)
            self.params['{0}scale_{s}'.format(cn, s=strain)] = self.scale[cn]
            if self.curvesLoaded:
                xRef, yRef, yRefErr = self.avgCurves[cn]
                yReft = np.interp(x, xRef, yRef)
            else:
                yReft = np.zeros_like(y)
            if x[-2] > self.params['{0}td_{s}'.format(cn, s=strain)]:
                self.params['{0}avgRes_{s}'.format(cn, s=strain)] = np.nanmean(
                    (y - yReft)[np.where(x > self.params['{0}td_{s}'.format(cn, s=strain)])])
                self.params['{0}stdTail_{s}'.format(cn, s=strain)] = np.nanstd(
                    (y)[np.where(x > self.params['{0}td_{s}'.format(cn, s=strain)])])
            else:
                self.params['{0}avgRes_{s}'.format(cn, s=strain)] = 0
                self.params['{0}stdTail_{s}'.format(cn, s=strain)] = 0
        else:
            self.params['{0}td_{s}'.format(cn, s=strain)] = 0
            self.params['{0}scale_{s}'.format(cn, s=strain)] = 1.
            xRef, yRef, yRefErr = self.avgCurves[cn]
            self.params['{0}avgRes_{s}'.format(cn, s=strain)] = -np.mean(yRef)
            self.params['{0}stdTail_{s}'.format(cn, s=strain)] = 0.

    def setLables(self):
        """
        Sets labels: RNAi, strain, field, date and label (embryo number)
        :return:
        """
        if self.folder is not None:
            imNames = glob.glob('{0}*.tif'.format(self.folder))
            self.batchFolder = '{1}batch-output{0}'.format(self.splitSym, self.folder)
            if len(imNames) > 0:
                name = '/'.join(imNames[2].split('\\'))
                self.RNAi, self.label, self.date, self.field = name.split(self.splitSym)[-1].split('_')[:4]
                self.label = name.split(self.splitSym)[-2]
            else:
                print('Can not find images, check connection')
                if os.path.exists(self.batchFolder):
                    ftmp = glob.glob('{0}*'.format(self.batchFolder))
                    ftmp1 = glob.glob('{0}*.*'.format(self.batchFolder))
                    ftmp = [t for t in ftmp if t not in ftmp1]
                    sort_nicely(ftmp)
                    ftmp, self.label, self.date = ftmp[0].split('_')[-9:-6]
                    ftmp = ftmp.split(self.splitSym)
                    self.RNAi = ftmp[-1]
                    if ftmp[-4] == self.date:
                        self.strain = ftmp[-5]
                    else:
                        self.strain = ftmp[-4]
                    self.field = None
                else:
                    printLog('Msg: Cannot determine lables, no expected files in the folder {0}'.format(self.folder))
                    self.date = None
                    self.label = self.folder.split(self.splitSym)[-2]
            if self.date is not None:
                ftmp = self.folder.split(self.splitSym)
                if ftmp[-3] == self.date:
                    self.strain = ftmp[-4]
                else:
                    self.strain = ftmp[-3]
        else:
            self.strain, self.RNAi, self.label, self.date = 'avg', 'avg', 'avg', 'avg'

    def setSize(self):
        self.radius = RADIUS_AVG  # The statistical average overestimates due to small tilt
        self.length = LENGTH_AVG

    def getSizeFromIm(self):
        """
        determine embryo size from image size (image has 6 pixel padding in y and 10 in x direction)
        :return:
        """
        imNames = glob.glob('{0}*.tif'.format(self.folder))
        if len(imNames) > 0:
            for name in imNames[:1]:
                im = cv2.imread(name, -1)
            y, x = im.shape
            y -= 6
            x -= 10
            y /= 2.
            x /= 2.
        return x, y

    def getMaxProj(self, t, thFlag=False):  # t is self.time (real time)
        """
        produces the maximum intensity projection image
        :param t: time point (time of development)
        :param thFlag: use automatic thresholds or not
        :return: 3 colored image ans numpy array and intensity center
        """
        slide = int(np.round(t / self.tScale + self.t0))  # convert t to slide number
        if slide > 30:
            slide = 30
        if slide < 0:
            black = True
            slide = 1
        else:
            black = False  # if time point is before imaging started return black image
        if self.image is None: self.loadImages()
        if thFlag:
            thG, thR = self.thresh
        else:
            thG, thR = 0, 0

        g = np.array(self.image.getMaxProject(slide, 0, thG), dtype=np.float32)  # /20000.
        xG, yG, zG = self.image.getIntCenter(slide, 0, thG)
        r = np.array(self.image.getMaxProject(slide, 1, thR), dtype=np.float32)  # /30000.
        xR, yR, zR = self.image.getIntCenter(slide, 1, thR)
        r = 1. * r / np.max(r)
        g = 1. * g / np.max(g)
        r[np.where(r > 1)] = 1.
        g[np.where(g > 1)] = 1.
        b = np.zeros_like(g)  # there is only green and red signal, the blue color contribution is a bunch of zeros
        if black:
            r = np.zeros_like(r)
            g = np.zeros_like(g)
        return np.dstack([r, g, b]), (xG, yG), (xR, yR)

    def checkGoodParam(self, pN):
        """
        checks if parameter value is useful based on its error being smaller than 50%
        :param pN: parameter name
        :return: True/False
        """
        pNErr = pN + 'Err'
        if pNErr in self.params:
            if np.isnan(self.params[pN]) or self.params[pN] > 0:
                return not np.isnan(self.params[pN]) and self.params[pNErr] / self.params[pN] < 0.5
            else:
                return True
        else:
            return not np.isnan(self.params[pN])

    def showInt(self, ch, fig=None):
        if fig is None: fig = myFigure()
        fig.plot(self.time - self.t0, self.tInt[ch], label=self.label)
        return fig

    def showDistCentMass(self, ax, fig=None, chs=range(3), setColor=True):
        if fig is None: fig = myFigure()
        color = ['g', 'r', 'y']
        for ch in chs:
            curveName = 'CoM{0}{1}'.format(ax, ['G', 'R', 'Y'][ch])
            if curveName in self.avgCurves: xRef, yRef, yRefErr = self.avgCurves[curveName]
            x, dist = self.getCurve(curveName)
            if setColor:
                fig.plot(x, dist, color=color[ch], label=self.label + 'scale={0}'.format(self.scale[curveName]))
                fig.plot(x[:self.cutoff[curveName]], dist[:self.cutoff[curveName]], 'k--',
                         label=self.label + 'scale={0}'.format(self.scale[curveName]))
                if curveName in self.avgCurves: fig.errorbar(xRef[::10], yRef[::10], yRefErr[::10])
            else:
                fig.plot(x, dist, label=self.label)
                fig.plot(x[:self.cutoff[curveName]], dist[:self.cutoff[curveName]], 'k--',
                         label=self.label + 'scale={0}'.format(self.scale[curveName]))
        fig.ylim((0, 1))
        fig.title(curveName)
        fig.legend()
        return fig

    def showMaxProj(self, t, bi=False):
        fig = myFigure()
        fig.markerSize = 100
        im, (xG, yG), (xR, yR) = self.getMaxProj(t)
        if bi:
            im[0] = im[0] > 0
            im[1] = im[1] > 0
        fig.imshow(im, colorbar=False)
        fig.scatter([xG], [yG], color='blue')
        fig.scatter([xR], [yR], color='pink')
        fig.noAxis()
        fig.setBGColor('k')
        fig.title('{2} t={0}, d={1}'.format(t, np.sqrt((xR - xG) ** 2 + (yR - yG) ** 2) / self.length, self.label))
        return fig

    def showSKProj(self, t):
        slide = max(0, int(np.round(t / self.tScale + self.t0)))
        if slide <= 31:
            if self.image is None: self.loadImages()
            thG, thR = np.array(self.thresh)
            rsk = self.image.getSkeleton(1, thR, slide)
            rsk[:, :, :rsk.shape[2] / 2] = 0
            gsk = self.image.getSkeleton(0, thG, slide)
            g = np.array(maxIntensProject(gsk)).astype(np.float)  # /20000.
            r = np.array(maxIntensProject(rsk)).astype(np.float)  # /30000.
            #             g = np.zeros_like(g)
            b = np.zeros_like(g)
            fig = myFigure()
            fig.imshow(np.dstack([r, g, b]), colorbar=False)
            fig.noAxis()
            fig.setBGColor('k')
            fig.title('{1} t={0}'.format((slide - self.t0) * self.tScale, self.label))
            SkeletonAnalysis.getLongestPath(rsk)
            return fig

    def showMaxProjHist(self, t):
        slide = np.round(t + self.t0)
        if slide <= 31:
            if self.image is None: self.loadImages()
            thG, thR = self.thresh
            figG = self.image.showMaxProjHist(t, 0, thG)
            figR = self.image.showMaxProjHist(t, 1, thR)
            figG.show()

    def reportParams(self):
        """
        print all parameter values
        :return:
        """
        for pN in PARAM_NAMES:
            print('{ch} {0} = {1}'.format(pN, self.params[pN], ch=self.checkGoodParam(pN)))


class MSEmbryo(Embryo):
    """ Class for MS strain """

    def __init__(self, folder=None, verb=True, check_version=True):
        Embryo.__init__(self, folder, verb, check_version)
        self.iniMSVars()

    def iniMSVars(self):
        self.nLambda = nLambda_MS
        self.tauScale = tauScale_MS
        self.setSize()  # sets self.radius, self.length
        self.loadAvgCurves()  # sets self.avgCurves

    def refresh(self):
        """
        Loads images and calculates all curves: intensities, MoIs, CoMs, etc., and runs updateParams
        :return:
        """
        Embryo.refresh(self)
        self.iniMSVars()
        self.tMove, self.movement = moveDetect.getMovementMS(self)
        self.params['tMove'] = self.tMove
        self.params['movement'] = self.movement
        self.setSigmParams()
        self.setThresh()
        self.setIntTh()
        self.setHead()
        self.updateParams()

    def updateParams(self):
        """
        Updates parameters and time aligns using curves data without loading images and saves
        :return:
        """
        printLog('{d}/{r}/{l} updating parameters...'.format(d=self.date, l=self.label, r=self.RNAi))
        for cn in AVG_CURVES_MS:
            self.cutoff[cn] = 1
            self.scale[cn] = 1.
            self.startPoint[cn] = 0
        self.setSigmParams()
        self.setCurveStartPoints()
        self.t0 = 0.
        self.tScale = 1.
        if (self.checkGoodParam('aR') and self.params['aR'] > 0) or (
                    self.checkGoodParam('aG') and self.params['aG'] > 0):
            timeAlign.alignTimeMS(self)
            if self.curvesLoaded:
                self.setScaleCutoff()
            self.setCMParams()
            self.setMIParams()
            self.setHeadParams()
            self.setLengthParams()
        self.params['tScale'] = self.tScale
        x, y = self.getSizeFromIm()
        # self.params['aspect'] = y / x
        self.params['tMove'] = (self.tMove - self.t0) * self.tScale
        self.params['mG'] = (self.params['mG'] - self.t0) * self.tScale
        self.params['mR'] = (self.params['mR'] - self.t0) * self.tScale
        self.params['mRmG'] *= self.tScale
        self.image = None
        for key in self.params:
            if not isinstance(self.params[key], int):
                self.params[key] = float(np.round(self.params[key], getSig(self.params[key])+5))
        self.save()

    def loadAvgCurves(self):
        """
        loads average curves data as a dictionary into self.avgCurves
        :return:
        """
        self.avgCurves = {}
        fileName = FILENAME_AVG_CURVES_MS_ADJTH
        if os.path.exists(fileName):
            print('loadAvgCurves: curves loaded')
            with open(fileName, 'rb') as input:
                self.avgCurves = pickle.load(input)
            self.curvesLoaded = True
        else:
            print('loadAvgCurves: NO curves loaded')
            self.curvesLoaded = False

    def setSize(self):
        """
        Calculates embryo size by using the bounding box around a signal obtained from the maximum intensity projection
        :return:
        """
        if self.radius is None or self.length is None:
            if self.image is None:
                self.loadImages()
            im = self.image.images[0, 9, 0]
            th = threshold_li(im)  # determine threshold
            radius, length = [], []
            for slide in range(10):  # check for 10 time points and use their average
                mask = (self.image.getMaxProject(slide, 0, th) > 0)  # gets binary mask of MIP for each time point
                label_im, nb_labels = ndimage.label(mask)  # labels objects on the mask
                if nb_labels > 1:  # find the largest object in the binary mask
                    sizes = ndimage.sum(mask, label_im, range(nb_labels + 1))
                    mask_size = sizes < np.max(sizes)
                    remove_pixel = mask_size[label_im]
                    label_im[remove_pixel] = 0
                    labels = np.unique(label_im)
                    label_im = np.searchsorted(labels, label_im)
                if np.max(label_im) > 0:
                    # Now that we have only one connect component, extract it's bounding box
                    slice_x, slice_y = ndimage.find_objects(label_im == 1)[0]
                    radius.append((slice_x.stop - slice_x.start) / 2.)
                    length.append((slice_y.stop - slice_y.start) / 2.)
            self.radius = np.round(np.mean(radius)).astype(int)  # calculate average size from 10 slides
            self.length = np.round(np.mean(length)).astype(int)

    def setCurveStartPoints(self):
        """
        populates self.startPoint, a time point when a curve becomes reliable, by identifying time point
        when there is significant signal
        :return:
        """
        if self.checkGoodParam('aR') and self.checkGoodParam('mR') and self.checkGoodParam('sR'):
            minR = self.params['mR'] - 2 * self.params['sR'] + 1  # time point of good red signal
        else:
            minR = 0
        if self.checkGoodParam('aG') and self.checkGoodParam('mG') and self.checkGoodParam('sG'):
            minG = self.params['mG'] - 2 * self.params['sG'] + 1  # time point of good green signal
        else:
            minG = 0
        indR = getSubIndex(self.time, minR, 31)[0]  # determine indices with good signal
        indG = getSubIndex(self.time, minG, 31)[0]
        if indG > 20:  # if nothing is good at least start at 20
            indG = 20
        if indR > 20:
            indR = 20
        for curveName in AVG_CURVES_MS:  # loop over all curves
            if curveName == 'tIntG' and self.checkGoodParam('aG'):  # require good green signal
                minInd = np.argmin(self.tInt[0][:10])
                self.startPoint[curveName] = \
                    np.where(self.tInt[0][minInd:] - self.params['bG'] + self.params['aG'] > 0)[0][0] + minInd
                # find point when the signal is above the baseline (b-a), but search only after minimum in the intensity
            elif curveName == 'tIntG':  # if fit is bad use average values
                minInd = np.argmin(self.tInt[0][:10])
                if (np.where(self.tInt[0][minInd:] - AVG_G_BASE_MS > 0)[0]).size > 0:
                    self.startPoint[curveName] = np.where(self.tInt[0][minInd:] - AVG_G_BASE_MS > 0)[0][0] + minInd
                else:
                    self.startPoint[curveName] = 0
            if curveName == 'tIntR' and self.checkGoodParam('aR'):  # same as in red
                minInd = np.argmin(self.tInt[1][:10])
                self.startPoint[curveName] = \
                    np.where(self.tInt[1][minInd:] - self.params['bR'] + self.params['aR'] > 0)[0][0] + minInd
            elif curveName == 'tIntR':
                minInd = np.argmin(self.tInt[1][:10])
                if (np.where(self.tInt[1][minInd:] - AVG_R_BASE_MS > 0)[0]).size > 0:
                    self.startPoint[curveName] = np.where(self.tInt[1][minInd:] - AVG_R_BASE_MS > 0)[0][0] + minInd
                else:
                    self.startPoint[curveName] = 0
            if curveName == 'MoI0G' or curveName == 'MoI1G' or curveName == 'CoM0G' or curveName == 'CoM1G':
                # reliable after good green signal
                self.startPoint[curveName] = max(indG, self.startPoint['tIntG']) + 2
            if curveName == 'MoI0R' or curveName == 'MoI1R' or curveName == 'CoM0R' or curveName == 'CoM1R' or curveName == 'SKR':
                # reliable after good red signal
                self.startPoint[curveName] = max(indR, self.startPoint['tIntR']) + 2
            if curveName == 'SKR' and np.sum(self.SK[1]) > 0:
                # reliable when there is something nonzero
                self.startPoint[curveName] = np.nonzero(self.SK[1])[0][0]
            if curveName == 'lengthR' and any(self.dims[1, :, 2] > 100):
                # checks that the length is at least 100 pixels
                self.startPoint[curveName] = np.where(self.dims[1, :, 2] > 100)[0][0]
            if curveName == 'headInt' and any(self.headInt > 0.05 * self.dz):
                # when head intensity is above empirically determined threshold
                self.startPoint[curveName] = np.where(self.headInt > 0.05 * self.dz)[0][0]
            elif curveName == 'headInt':
                self.startPoint[curveName] = np.where(np.ones_like(self.headInt) - np.isnan(self.headInt))[0][0]
        for key in self.startPoint:  # make all start points integer
            self.startPoint[key] = int(self.startPoint[key])

    def setSigmParams(self):
        """
        Sets parameters relevant to intensity parametrization
        :return:
        """
        if self.params['maxG'] > 3e8:  # check that there is some signal in green
            popt, perr = self.fitSigm2Int(0)
            self.params['aG'] = popt[0]
            self.params['bG'] = popt[1]
            self.params['mG'] = popt[2]
            self.params['sG'] = popt[3]
            self.params['rG'] = popt[4]
            self.params['aGErr'] = perr[0]
            self.params['bGErr'] = perr[1]
            self.params['mGErr'] = perr[2]
            self.params['sGErr'] = perr[3]
            self.params['rGErr'] = perr[4]
        else:
            self.params['aG'] = 0.
            self.params['bG'] = 0.
            self.params['mG'] = np.nan
            self.params['sG'] = 0.
            self.params['rG'] = np.nan
            self.params['aGErr'] = 0.
            self.params['bGErr'] = 0.
            self.params['mGErr'] = np.nan
            self.params['sGErr'] = 0.
            self.params['rGErr'] = 0.
            self.params['movement'] = 0
            self.movement = False
            self.tMove = np.nan
        if self.params['maxR'] > 1.5e8:  # check that there is some signal in red
            popt, perr = self.fitSigm2Int(1)
            self.params['aR'] = popt[0]
            self.params['bR'] = popt[1]
            self.params['mR'] = popt[2]
            self.params['sR'] = popt[3]
            self.params['rR'] = popt[4]
            self.params['aRErr'] = perr[0]
            self.params['bRErr'] = perr[1]
            self.params['mRErr'] = perr[2]
            self.params['sRErr'] = perr[3]
            self.params['rRErr'] = perr[4]
        else:
            self.params['aR'] = 0.
            self.params['bR'] = 0.
            self.params['mR'] = np.nan
            self.params['sR'] = 0.
            self.params['rR'] = np.nan
            self.params['aRErr'] = 0.
            self.params['bRErr'] = 0.
            self.params['mRErr'] = np.nan
            self.params['sRErr'] = 0.
            self.params['rRErr'] = 0.
            self.params['movement'] = 0
            self.movement = False
            self.tMove = np.nan
        if self.checkGoodParam('mR') and self.checkGoodParam('mG'):  # if both Ms are good calculate their difference
            self.params['mRmG'] = self.params['mR'] - self.params['mG']

    def setMIParams(self):
        """
        sets parameters relevant to moment of inertia
        :return:
        """
        curveNames = ['MoI{0}{1}'.format(ax, c) for c in ['G', 'R'] for ax in range(2)]
        for cn in curveNames:
            self.setCM_MoIParams(cn, 'MS')

    def setCMParams(self):
        """
        sets parameters relevant to center of mass
        :return:
        """
        for ch in ['R', 'G']:
            curveNames = ['CoM{0}{1}'.format(ax, ch) for ax in range(2)]
            for cn in curveNames:
                self.setCM_MoIParams(cn, 'MS')

    def setHeadParams(self):
        """
        sets parameters relevant to head intensity
        :return:
        """
        if self.headInt[-1] < 0.65:
            ind = np.where(np.roll(self.headInt, -1) - self.headInt < -0.8)[0]
            if ind.size > 0:
                ind = ind[-1] - 1
            else:
                ind = self.headInt.size
        else:
            ind = 31
        x, y = self.getCurve('headInt')
        if np.nanmin(y) < 0.2 and np.nanmax(y) > 0.5:
            popt, perr = fitSigmBump(x, y)
            #         popt, perr = fitSigmBump(self.time[:ind], self.headInt[:ind])
            self.params['aSigHead'] = popt[0]
            self.params['mSigHead'] = popt[1]
            self.params['sSigHead'] = popt[2]
            self.params['aGaussHead'] = popt[3]
            self.params['mGaussHead'] = popt[4]
            self.params['sGaussHead'] = popt[5]
            # self.params['aSigHeadErr'] = perr[0]
            # self.params['mSigHeadErr'] = perr[1]
            # self.params['sSigHeadErr'] = perr[2]
            # self.params['aGaussHeadErr'] = perr[3]
            # self.params['mGaussHeadErr'] = perr[4]
            # self.params['sGaussHeadErr'] = perr[5]
            self.params['maxHead'] = np.nanmax(self.headInt)
            x, y = self.getCurve('headInt')
            self.params['tailHead'] = np.nanmedian(y[x > 10])
            self.params['scaleHead'] = self.scale['headInt']
            self.params['devHead'] = x[self.cutoff['headInt'] - 1]
        else:
            self.params['aSigHead'] = np.nanmax(y)
            self.params['mSigHead'] = np.nan
            self.params['sSigHead'] = np.nan
            self.params['aGaussHead'] = 0
            self.params['mGaussHead'] = np.nan
            self.params['sGaussHead'] = np.nan
            # self.params['aSigHeadErr'] = 0
            # self.params['mSigHeadErr'] = np.nan
            # self.params['sSigHeadErr'] = np.nan
            # self.params['aGaussHeadErr'] = 0
            # self.params['mGaussHeadErr'] = np.nan
            # self.params['sGaussHeadErr'] = np.nan
            self.params['maxHead'] = np.nanmax(self.headInt)
            x, y = self.getCurve('headInt')
            # self.params['tailHead'] = 0  #these were enabled prior to 040722. Changed to negate the if/else to get vals for this param
            self.params['tailHead'] = np.nanmedian(y[x > 10])
            # self.params['scaleHead'] = np.nan #these were enabled prior to 040722. Changed to negate the if/else to get vals for this param
            self.params['scaleHead'] = self.scale['headInt']
            # self.params['devHead'] = 0  #these were enabled prior to 040722. Changed to negate the if/else to get vals for this param
            self.params['devHead'] = x[self.cutoff['headInt'] - 1]


    def setLengthParams(self):
        """
        sets parameters relevant to length in red
        :return:
        """
        self.params['scaleLength'] = self.scale['lengthR']
        x, y = self.getCurve('lengthR')
        self.params['devLength'] = x[self.cutoff['lengthR'] - 1]
        self.params['tailLength'] = np.nanmedian(y[x > 10])

    def setThresh(self):
        """
        set self.thresh, an automatic threshold for intensity, determined by the median of yen thresholds of 3
        consecutive images near onset of green signal
        :return:
        """
        for ch in range(2):
            t = []
            if self.checkGoodParam('mG') and self.checkGoodParam('aG'):
                t.append(self.params['mG'] - 2. * self.params['sG'])
            if self.checkGoodParam('mR') and self.checkGoodParam('aR'):
                t.append(self.params['mR'] - 1.5 * self.params['sR'])
            if len(t) == 0:
                t = 0
            elif ch == 1:
                t = min(t)
            else:
                t = t[0]
            if t < 0:
                t = 0
            if self.image is None:
                self.loadImages()
            if 25 > t > 0:
                ts = [t - 1, t, t + 1]
            elif t == 0:
                ts = range(3)
            else:
                ts = [t - 2, t - 1, t]
            th = []
            for t in ts:  # loop over 3 consecutive time points
                im = self.image.getMaxProject(np.round(t).astype(int), ch, 0)
                if ch == 0:  # if green signal, cut the edges to remove possible interference from adjacent embryos
                    l = im.shape[1]
                    im = im[10:-10, l / 4:3 * l / 4]
                thtmp = threshold_yen(im)
                if 1000 < thtmp < 10000:  # check that threshold is within normal range
                    th.append(thtmp)
            if len(th) > 0:
                th = np.median(th)
            else:
                th = 5000
            self.thresh[ch] = np.round(th).astype(int)

    def setIntTh(self):
        """
        sets self.tIntTh, the total intensity of every image above a threshold.
        Note that threshold is subtracted from each image before evaluation.
        :return:
        """
        for c in range(2):
            self.tIntTh[c, :] = self.image.getAllIntens(c, self.thresh[c])
            if c == 1:
                #                 self.SK[c,:] = np.sum(self.image.getSkeleton(c, self.thresh[c]), axis=(1,2,3))
                for slide in self.time.astype(np.int):
                    #                     rsk = self.image.getSkeleton(c, self.thresh[c], slide)
                    #                     rsk[:,:,:rsk.shape[2]/2]=0
                    #                     self.SK[c,slide] = SkeletonAnalysis.getLongestPath(rsk)
                    self.dims[c, slide] = np.round(self.image.getDims(c, self.thresh[c], slide)).astype(int)

    def setHead(self):
        """
        Calculates number of pixels above threshold projected on the end-on view from 1/3 of the embryo and nurmalized
        by the area of the embryo.
        :return:
        """
        from scipy.ndimage import zoom

        if self.image is None:
            self.loadImages()
        for t in self.time:
            ims = self.image.images[int(t), :, 0]
            length = ims[0].shape[1]
            ims = ims[:, :, :int(length / 3.)]
            imProj = np.max(ims, axis=2)
            imProj = (imProj > self.thresh[0])

            # if (t-self.t0)*self.tScale > 8:
            #     fig = myFigure()
            #     fig.makeBlack()
            #     fig.set_axes_equal()
            #     fig.noAxis()
            #     fig1 = myFigure()
            #     fig1.makeBlack()
            #     fig1.set_axes_equal()
            #     fig1.noAxis()
            #     imProj = zoom(imProj, (self.dz, 1), order=0)
            #     fig.imshow(imProj, bw=True, colorbar=False)
            #     fig1.imshow((self.image.getMaxProject(int(t), 0, self.thresh[0])>0), bw=True, colorbar=False)
            #     print('head Int = {}'.format(np.sum(imProj)/self.radius**2))
            #     fig.show()
            #     fig.save('Z:/EndOn_headInt_{d}_{l}_slide={t}.png'.format(d=self.date,l=self.label,t=t))
            # fig1.save('Z:/MIP_headInt_{d}_{l}_slide={t}.png'.format(d=self.date,l=self.label,t=t))

            self.headInt[int(t)] = np.sum(imProj)
        self.headInt /= self.radius ** 2 / self.dz  # normalize by embryo radius square and z sampling
        self.headInt = np.round(self.headInt, 3)

    def getDistFromGreen(self, ch, axis):
        """
        Calculates the euclidean distance from CoM in green to other channels. Only projection on the AP axis is
        considered for AP distance component.
        :param ch: channel (int)
        :param axis: end-o/ap (0/1)
        :return: list of distances for every time point
        """
        d = np.zeros(31)
        xe, ye, ze = self.getEmbCenter()
        for slide in range(31):
            if self.loaded or np.max(self.CMPos[ch, slide]) > 0:
                xc, yc, zc = self.CMPos[ch, slide]
            else:
                if self.image is None:
                    self.loadImages()
                xc, yc, zc = self.image.getIntCenter(slide, ch, self.thresh[ch])
                zc *= self.dz
                self.CMPos[ch, slide] = np.array([xc, yc, zc])
            if axis == 0:
                d[slide] = np.sqrt((yc - ye) ** 2 + (zc - ze) ** 2) / self.radius
            elif axis != 0:
                d[slide] = np.sqrt((xc - xe) ** 2) / self.length
        return np.round(d, 3)

    def getEmbCenter(self):
        if self.embCenter is None or np.isnan(self.embCenter).any():
            if self.image is None:
                self.loadImages()
            xc, yc, zc = [], [], []
            im = self.image.images[0, 9, 0]
            th = threshold_li(im)
            for slide in range(10):
                x, y, z = self.image.getIntCenter(slide, 0, th)
                xc.append(x)
                yc.append(y)
                zc.append(z * self.dz)
            self.embCenter = np.round([np.mean(xc), np.mean(yc),
                              np.mean(yc)]).astype(int)  # use y in place of z, because depth attenuation makes z center inacurate.
        return self.embCenter

    #     def getDistFromGreen(self, ch, axis):
    #         if self.params['mG']!=np.nan and self.params['mR']!=np.nan:
    #             d = np.zeros(31)
    #             for slide in range(31):
    #                 if self.loaded or np.max(self.CMPos[1, slide])>0: xc, yc, zc = self.CMPos[1,slide]
    #                 else:
    #                     xc, yc, zc = self.image.getIntCenter(slide,1,self.thresh[1])
    #                     zc*=self.dz
    #                     self.CMPos[1,slide]= np.array([xc, yc, zc])
    #                 if self.loaded or np.max(self.CMPos[0, slide])>0: xe, ye, ze = self.CMPos[0,slide]
    #                 else:
    #                     xe, ye, ze = self.image.getIntCenter(slide,0,self.thresh[0])
    #                     ze*=self.dz
    #                     self.CMPos[0,slide]= np.array([xe, ye, ze])
    #                 if axis==0: d[slide] = np.sqrt((yc-ye)**2+(zc-ze)**2)/self.radius
    #                 elif axis!=0: d[slide] = np.sqrt((xc-xe)**2)/self.length
    #             return d
    #         else: return 0

    def getAllMoments(self, ch, axis):
        """
        outputs moments of inertia (populates self.MoI)
        :param ch: channel (G-0, R-1, DIC-2)
        :param axis: end-on/ap (0/1)
        :return: moments of inertia for a given channel around the specified axis
        """
        if self.loaded or np.max(self.MoI[ch, axis]) > 0:
            m = self.MoI[ch, axis]
        else:
            if self.image is None:
                self.loadImages()
            m = self.image.getAllMoments(ch, axis, th=self.thresh[ch])
            self.MoI[ch, axis] = m
        m = m / self.tIntTh[ch]
        m[np.isnan(m)] = 0
        if axis == 0:
            m /= self.radius ** 2
        else:
            m /= self.length ** 2
        return np.round(m, 3)

    def getCurveNames(self):
        return np.array(AVG_CURVES_MS)

    def getCurve(self, curveName):
        """
        Return curve calculated for this embryo, all values before start point are nans
        :param curveName: name
        :return: x (normalized and aligned time), y (value)
        """
        if curveName not in AVG_CURVES_MS:
            raise NameError('error passing wrong name to fit2Avg: ' + curveName)
        if curveName == 'tIntG':
            if self.checkGoodParam('aG'):
                y = np.copy(self.tInt[0]) - self.params['bG'] + self.params['aG']
            else:
                y = np.copy(self.tInt[0]) - AVG_G_BASE_MS
            y[np.where(y < 0)] = 0
        elif curveName == 'tIntR':
            if self.checkGoodParam('aR'):
                y = np.copy(self.tInt[1]) - self.params['bR'] + self.params['aR']
            else:
                y = np.copy(self.tInt[1]) - AVG_R_BASE_MS
            y[np.where(y < 0)] = 0
        elif curveName == 'CoM0R':
            y = self.getDistFromGreen(1, 0)
        elif curveName == 'CoM1R':
            y = self.getDistFromGreen(1, 1)
        elif curveName == 'CoM0G':
            y = self.getDistFromGreen(0, 0)
        elif curveName == 'CoM1G':
            y = self.getDistFromGreen(0, 1)
        elif curveName == 'MoI0G':
            y = self.getAllMoments(0, 0)
        elif curveName == 'MoI1G':
            y = self.getAllMoments(0, 1)
        elif curveName == 'MoI0R':
            y = self.getAllMoments(1, 0)
        elif curveName == 'MoI1R':
            y = self.getAllMoments(1, 1)
        elif curveName == 'SKR':
            y = np.copy(self.SK[1])
        elif curveName == 'lengthR':
            y = np.copy(self.dims[1, :, 2])
        elif curveName == 'headInt':
            y = np.copy(self.headInt)
        y[:self.startPoint[curveName]] = np.nan
        return self.tScale * (self.time - self.t0), y * self.scale[curveName]

    def fitSigm2Int(self, ch):
        """
        Fitting total intensity values with a sigmoidal.
        :param ch: channel
        :return: fitted parameters
        """
        def sigmoidal(x, a, b, m, s, r):  # define sigmoidal equation. Note the baseline is b-a and the plateau is b.
            return b - a * (1. + np.exp((x - m) / s)) ** (-r)

        def getErr(y):
            """
            function to estimate error from second derivative (high for rapid up/down) of the data.
            :param y: data
            :return: errors
            """
            dy = (np.roll(y, -1) - y)[:-1]  # find difference between consecutive points
            ddy = (np.roll(dy, -1) - dy)[:-1]
            err = np.interp(np.arange(y.size), np.arange(ddy.size) + 1, np.abs(ddy))
            #     err= median_filter(err,3)
            errMin = 0.5e8
            err[np.where(err < errMin)] = errMin  # re-assign values in err that are below min to be 0.5e8
            # (empirically determined)
            return err

        show = False
        x = self.time
        y = self.tInt[ch]
        endInd = min(y.size, self.tMove + 2)
        dy = (np.roll(y, -1) - y)[:-1]
        ind0 = np.where(dy > 0)[0][0]  # find the first time point where the slope is positive to avoid initial drop
        x = x[max(0, ind0 - 3):]  # start fitting three points before the drop
        y = y[max(0, ind0 - 3):]
        dy = dy[max(0, ind0 - 3):]

        # This block defines stop point for sigmoidal fit in cases when intensity drops
        ind = np.where(dy < -3e8)[0]  # index positions (in dy) when embryo has hatched
        if len(ind) > 0 and ind[0] > 20 - ind0:
            x = x[:ind[0] + 1]  # + 1 due to dy being 1 short compared to y
            y = y[:ind[0] + 1]
            if ind[0] + 1 < self.stopInt:
                self.stopInt = ind[0] + 1
        elif y[-1] < np.max(y) / 2.:
            ind = np.where(dy * (y > np.max(y) / 2.)[:-1] > 0)[0]  # indices when curve is going up and above max/2
            if len(ind) > 0:
                x = x[:ind[-1] + 1]
                y = y[:ind[-1] + 1]
                if ind[0] + 1 < self.stopInt:
                    self.stopInt = ind[0] + 1

                #         endInd = y.size
        endInd = min(y.size, endInd)
        r0 = 1.  # exponent
        b0 = y[endInd - 1]  # np.max(y) #plateau at the end
        minInd = np.argmin(y[:y.size / 2])
        if endInd <= minInd:
            printLog('could not fit intensity')
            popt = np.ones(5) * np.nan  # expects 5 params: a, b, m, s, r
            perr = np.ones(5) * np.nan
            residsEarly = np.nan
            residsLate = np.nan
            moveT = np.nan
            fig = myFigure()
            return popt, perr

        a0 = b0 - np.min(y)  # b-a plateau at the beginning
        m0 = x[minInd:endInd][
            np.argmin(np.abs(y[minInd:endInd] - (b0 - a0 / 2.)))]  # x position of the middle intensity
        s0 = 3  # emperically determined width
        p0 = (a0, b0, m0, s0, r0)  # initially guessed parameters
        params = Parameters()  # define parameters to be fitted
        params.add('a', p0[0], min=0.5 * (max(y[:endInd]) - min(y[:endInd])),
                   max=1.2 * (max(y[:endInd]) - min(y[:endInd])))  # min/max defines boundaries for a parameter
        params.add('base', p0[1] - p0[0], min=0)
        params.add('b', p0[1], expr='base+a')
        params.add('m', p0[2], min=0, max=max(x[:endInd]))
        params.add('s', p0[3], min=0, max=max(x[:endInd]))
        params.add('r', p0[4], min=0, max=10, vary=False)
        mod = Model(sigmoidal)  # assign curve to be fitted
        err = getErr(y)
        #         err= np.ones_like(getErr(y))
        if show:
            color = ['g', 'r']
            fig = myFigure()
            fig.errorbar(self.time, self.tInt[ch], getErr(self.tInt[ch]), join=False, color=color[ch], alpha=0.5)
            fig.errorbar(x[:endInd], y[:endInd], err[:endInd], join=False, color=color[ch])
        try:
            # fits sigmoidal to the data
            res = mod.fit(y[:endInd], x=x[:endInd], params=params, weights=1. / err[:endInd] ** 2)
            popt = np.array([res.best_values['a'], res.best_values['b'], res.best_values['m'], res.best_values['s'],
                             res.best_values['r']])
            #             printLog(report_fit(res))
            if show: fig.plot(x, sigmoidal(x, *popt), color=color[ch])
        except Exception as e:
            printLog('could not fit intensity: {0}'.format(str(e)))
            if show: fig.plot(x, sigmoidal(x, *p0), 'k')
            popt = np.ones(5) * np.nan
            perr = np.ones(5) * np.nan
            residsEarly = np.nan
            residsLate = np.nan
            moveT = np.nan
        try:
            perr = np.array(
                [res.params['a'].stderr, res.params['b'].stderr, res.params['m'].stderr, res.params['s'].stderr,
                 res.params['r'].stderr])
            perr[np.where(perr == 0)] = np.nan
        except:
            printLog('bad parameter errors in intensity fit')
            perr = np.nan * np.ones_like(popt)
        if show:
            print(report_fit(res))
            plt.show()
        return popt, perr

    def show(self, curveName, fig=None, setColor=True):
        """
        plots curves on a figure
        :param curveName: name (str)
        :param fig: myFigure object (optional)
        :param setColor: if true green is used for the 1st channel and red is used for the 2nd. If False, colors are
        random.
        :return: myFigure object
        """
        if curveName == 'tIntG': return self.showSigmFit(fig, [0], setColor)
        if curveName == 'tIntR': return self.showSigmFit(fig, [1], setColor)
        if curveName == 'CoM0R': return self.showDistCentMass(0, fig, [1], setColor)
        if curveName == 'CoM1R': return self.showDistCentMass(1, fig, [1], setColor)
        if curveName == 'CoM0G': return self.showDistCentMass(0, fig, [0], setColor)
        if curveName == 'CoM1G': return self.showDistCentMass(1, fig, [0], setColor)
        if curveName == 'MoI0G': return self.showMOI(0, fig, [0], setColor)
        if curveName == 'MoI1G': return self.showMOI(1, fig, [0], setColor)
        if curveName == 'MoI0R': return self.showMOI(0, fig, [1], setColor)
        if curveName == 'MoI1R': return self.showMOI(1, fig, [1], setColor)
        if curveName == 'SKR': return self.showSkeleton(1, fig)
        if curveName == 'lengthR': return self.showDims(1, fig)
        if curveName == 'headInt': return self.showHead(fig)

    def showMOI(self, axis, fig=None, chs=range(2), setColor=True):
        """
        plots MOI for MS. If setColor is False, uses autocolor (would use this if plotting lots of embryos on one plot)
        :param axis: axis (0-end-on, 1-ap)
        :param fig: myFigure object (optional)
        :param chs: channel (if None, plots both)
        :param setColor: use channel color
        :return: myFigure object
        """
        if not fig: fig = myFigure()  # if no figure exists, creates figure
        color = ['g', 'r']  # specifies colors if setColor is True
        for ch in chs:  # iterates through channels
            curveName = 'MoI{0}{1}'.format(axis, ['G', 'R'][ch])
            t, MOI = self.getCurve(curveName)
            if setColor:
                fig.plot(t, MOI, label=self.label, color=color[
                    ch])  # plots in green, red... use for one embryo, or population if dont need to ID individual embs
                fig.plot(t[:self.cutoff[curveName]], MOI[:self.cutoff[curveName]], 'k--', label=self.label)
                if curveName in self.avgCurves:
                    x, y, yerr = self.avgCurves[curveName]
                    fig.errorbar(x[::10], y[::10], yerr[::10], color='gray')
            else:
                fig.plot(t, MOI,
                         label=self.label)  # plots in multiple colors... use for multiple embryos to trace back individual embs
            fig.title(curveName)
            fig.xlim([0, 30])
            fig.ylim([0, 1])
            fig.legend()
        return fig

    def showSigmFit(self, fig=None, chs=[0, 1], setColor=True):
        """
        plots sigmoidal fit to intensity for MS. If set color arg is False, uses autocolor (would use this if plotting
        lots of embryos on one plot)
        :param fig: myFigure object (optional)
        :param chs: channel (if None, plots both)
        :param setColor: use channel color
        :return: myFigure object
        """
        def sigmoidal(x, a, b, m, s, r):
            #             return b-a*(1.+np.exp((x-m)/s))**(-r)
            return a - a * (1. + np.exp((x - m) / s)) ** (-r)

        if fig is None: fig = myFigure()
        for ch in chs:
            curveName = 'tInt{col}'.format(col=['G', 'R'][ch])
            x, y = self.getCurve(curveName)
            if ch == 0:
                keys = ['aG', 'bG', 'mG', 'sG', 'rG']
                color = 'g'
            elif ch == 1:
                keys = ['aR', 'bR', 'mR', 'sR', 'rR']
                color = 'r'
            dy = (np.roll(self.tInt[0], -1) - self.tInt[0])[:-1]
            ind = np.where(dy < -3e8)[0]  # check that worm didn't hatch
            if len(ind) > 0:
                x = x[:ind[0] + 1]
                y = y[:ind[0] + 1]
            a, b, m, s, r = [self.params[k] for k in keys]
            if setColor:
                fig.plot(x, y, color=color)
                fig.plot(x[:self.cutoff[curveName]], y[:self.cutoff[curveName]], 'k--')
                fig.plot(x, sigmoidal(x / self.tScale, a, b, m / self.tScale, s, r) * self.scale[curveName], color='k')
                fig.plot(((self.tMove - self.t0) * self.tScale, (self.tMove - self.t0) * self.tScale), (0, a),
                         color='b')
                #                 xl = np.arange(0,5,0.1)
                #                 alpha = a/4./s
                #                 if ch==1: fig.plot(xl, xl*alpha/self.tScale, color='b')
                if curveName in self.avgCurves:
                    x, y, yerr = self.avgCurves[curveName]
                    fig.errorbar(x[::10], y[::10], yerr[::10], color=color)
                    fig.plot((T0_2_MOVE_AVG_MS, T0_2_MOVE_AVG_MS), (0, a), color='gray')
            else:
                fig.plot(x, y, label=self.label)
                fig.plot(x[:self.cutoff[curveName]], y[:self.cutoff[curveName]], 'k--')
                # fig.plot(((self.tMove - self.t0) * self.tScale, (self.tMove - self.t0) * self.tScale), (0, a),
                #          color='b')
            fig.title(curveName)
        return fig

    def showSkeleton(self, ch, fig=None):
        """
        plots longest shortest path of skeletonized signal.
        :param ch: channel (should be 1, not used)
        :param fig: myFigure object (optional)
        :return: myFigure object
        """
        if fig is None: fig = myFigure()
        x, y = self.getCurve('SKR')
        fig.title('SKR')
        fig.plot(x, y, label=self.label)
        fig.plot(x[:self.cutoff['SKR']], y[:self.cutoff['SKR']], 'k--')
        if self.curvesLoaded:
            xR, yR, yerrR = self.avgCurves['SKR']
            fig.errorbar(xR[::10], yR[::10], yerrR[::10], color='k')
        # yRt = np.interp(x, xR, yR)
        #             print('{d}/{l} res={r}'.format(d=self.date, l=self.label, r=np.nansum((y-yRt)**2)))
        return fig

    def showDims(self, ch, fig=None):
        """
        plots length of red signal.
        :param ch: channel (should be 1, not used)
        :param fig: myFigure object (optional)
        :return: myFigure object
        """
        if fig is None: fig = myFigure()
        x, y = self.getCurve('lengthR')
        fig.plot(x, y, label=self.label)
        fig.plot(x[:self.cutoff['lengthR']], y[:self.cutoff['lengthR']], 'k--')
        if 'lengthR' in self.avgCurves:
            x, y, yerr = self.avgCurves['lengthR']
            fig.errorbar(x[::10], y[::10], yerr[::10], color='k')
        fig.title('length of red')
        return fig

    def showHead(self, fig=None):
        """
        plots area of the end-on projection of the green head fluorescence.
        :param fig: myFigure object (optional)
        :return: myFigure object
        """
        def cost(x, a, m1, s1, b, m2, s2):  # sigmoidal with a bump (sigmoidal + gaussian)
            return a / (1. + np.exp(-(x - m1) / s1)) + b * np.exp(-(x - m2) ** 2 / s2 ** 2)

        if fig is None: fig = myFigure()
        x, y = self.getCurve('headInt')
        fig.plot(x, y, label=self.label)
        fig.plot(x[:self.cutoff['headInt']], y[:self.cutoff['headInt']], 'k--')
        xt = np.arange(x[0] / self.tScale + self.t0, x[-1] / self.tScale + self.t0, 0.1)
        # fig.plot((xt - self.t0) * self.tScale,
        #          self.scale['headInt'] * cost(xt, self.params['aSigHead'], self.params['mSigHead'],
        #                                       self.params['sSigHead'], \
        #                                       self.params['aGaussHead'], self.params['mGaussHead'],
        #                                       self.params['sGaussHead']))
        if 'headInt' in self.avgCurves:
            xR, yR, yerrR = self.avgCurves['headInt']
            fig.errorbar(xR[::10], yR[::10], yerrR[::10], color='k')
        #     yRt = np.interp(x, xR, yR)
            # print('{d}/{l} res={r}'.format(d=self.date, l=self.label, r=np.nansum((y - yRt) ** 2)))
        fig.title('head intensity')
        return fig

    def fix_head_params_error(self):
        '''
        This function updates headInt and tailHead parameters when an embryo object has been called.
         These required updating due to a bug in the code for pulling params from the headInt curve.
        :return:
        '''
        x, y = self.getCurve('headInt')
        # print(self.params['tailHead'])
        # print(self.params['devHead'])

        # self.params['tailHead'] = 0  #these were enabled prior to 040722. Changed to negate the if/else to get vals for this param
        tailHead_fix = np.nanmedian(y[x > 10])
        self.params['tailHead'] = tailHead_fix

        # self.params['devHead'] = 0  #these were original settings within the else statement. Enabled prior to 040722. Changed to negate the if/else to get vals for this param
        devHead_fix = x[self.cutoff['headInt'] - 1]
        self.params['devHead']=devHead_fix

        # print(self.params['tailHead'])
        # print(self.params['devHead'])
        self.save()

def fix_headInt_params_by_rnai_number():
    from db_utils_embryos import get_embryo_number_by_rna
    from varLookup import FOLDER_IN
    from os.path import exists


    rnaiList = [i for i in range(216, 255)] #RNAi numbers to update parameter values for (this has been run for 1-510)
    for rna in rnaiList:
        rnai = rna
        print(rnai)
        # r1 = RNAiClass(rnai)
        embs_list = get_embryo_number_by_rna(rnai, strain='MS') #looks at mySQL for specified RNAi condition and pulls a list of emb numbers to identify which objects to pull
        for emb in embs_list:  # emb is an integer
            folders = [FOLDER_IN + 'cropped/EMBD{r:04}/MS/Emb{e}/'.format(r=rnai, e=emb)]  # identify embryo object
            file_exists = exists(folders[0])
            if file_exists:
                embryo_obj = MSEmbryo(folders[0], check_version=False)  # access embryo object
                embryo_obj.fix_head_params_error() #updates parameters
            else:
                print('------------EMBRYO MISSING RNA={x}, emb={y}---------'.format(x=rnai,y=emb))


def fix_headInt_params_control_embs():
    from db_utils_embryos import get_embryo_number_by_rna, get_file_path_by_rna
    from varLookup import FOLDER_IN
    from os.path import exists

    rnaiList = [0]  # RNAi numbers to update parameter values for
    for rna in rnaiList:
        rnai = rna
        file_path_list = get_file_path_by_rna(rnai, strain='MS') # looks at mySQL for control condition and pulls a list of filepaths to identify which objects to pull
        # file_path_list = ['Z:/Emb2/'] #exclude
        exclude_list= ['140220, 140219', '140213', '140206', '140205', '140129', '140117', '140123', '140115', '140110','131219','131218','131213','131122', '150513', '140627', '150618', '131212', '131220', '140109', '140122']
        for emb in file_path_list:  # emb is an integer
            # folders = [FOLDER_IN + 'cropped/EMBD{r:04}/MS/Emb{e}/'.format(r=rnai, e=emb)]  # identify embryo object
            if emb[25:31] in exclude_list:
                pass
            else:
                file_exists = exists(emb)
                if file_exists:
                    embryo_obj = MSEmbryo(emb, check_version=False)  # access embryo object
                    embryo_obj.fix_head_params_error()  # updates parameters

            # else:
            #     print('------------EMBRYO MISSING RNA={x}, emb={y}---------'.format(x=rnai, y=emb))

class GSEmbryo(Embryo):
    """Class for GS strain"""

    def __init__(self, folder=None, verb=True, check_version=True):
        Embryo.__init__(self, folder, verb, check_version)
        self.iniGSVars()

    def iniGSVars(self):
        self.setSize()  # sets self.radius, self.length
        self.nLambda = nLambda_GLS
        self.tauScale = tauScale_GLS
        self.green = np.zeros_like(self.time)
        self.red = np.zeros_like(self.time)
        self.yellow = np.zeros_like(self.time)
        self.greenPos = None
        self.redPos = None
        self.yellowPos = None
        self.loadAllSpots()
        self.loadAvgCurves()  # sets self.avgCurves

    def refresh(self):
        Embryo.refresh(self)
        self.iniGSVars()
        self.tMove, self.movement = moveDetect.getMovementGS(self)
        self.updateParams()

    def updateParams(self):
        printLog('{d}/{r}/{l} updating parameters...'.format(d=self.date, l=self.label, r=self.RNAi))
        for cn in AVG_CURVES_GLS:
            self.cutoff[cn] = 0
            self.scale[cn] = 1.
            self.startPoint[cn] = 0
        if self.folder is not None:
            self.loadAllSpots()
        self.setSigmParams()
        self.t0 = 0.
        self.tScale = 1.
        if self.spotsLoaded:
            if np.sum(np.concatenate((self.green > 20, self.red > 20, self.yellow > 20))) > 10:
                if np.isnan(self.params['residsLateY']) or self.params['residsLateY'] < 30: self.movement = False
                if np.isnan(self.params['Ytail']) or np.abs(
                                self.params['Ytail'] - 110) < 10: self.movement = False  # 3 fold arrest
                self.setCurveStartPoints()
                timeAlign.alignTimeGLS(self)
                if self.curvesLoaded:  # curvesLoaded checks to see if avgCurves exist
                    self.setScaleCutoff()  # determines t0, td and x,y scale
            self.params['movement'] = self.movement
            self.setCMParams()
            self.setMIParams()
            self.params['tScale'] = self.tScale
            x, y = self.getSizeFromIm()
            # self.params['aspect'] = y / x
            self.params['tMove'] = (self.tMove - self.t0) * self.tScale
            self.params['mG'] = (self.params['mG'] - self.t0) * self.tScale
            self.params['mR'] = (self.params['mR'] - self.t0) * self.tScale
            self.params['mY'] = (self.params['mY'] - self.t0) * self.tScale
            self.params['mYmG'] *= self.tScale
            self.params['mRmG'] *= self.tScale
            self.params['mYmR'] *= self.tScale
            for key in self.params:
                if not isinstance(self.params[key], int):
                    self.params[key] = float(np.round(self.params[key], getSig(self.params[key]) + 5))
            self.save()
        self.image = None

    def loadAllSpots(self):
        """
        Loads spots number and positions from batch generated csv files. If any spots folders are missing,
        then self.spotsLoaded is set to False and none are loaded.
        :return:
        """

        def loadSpots(fileName):
            csvFile = csv.reader(open(fileName, 'rU'), delimiter=',')
            nSpots = np.zeros_like(self.time)
            flag = False
            for row in csvFile:
                if flag == True:
                    nSpots[np.where(self.time == int(row[2]) - 1)] = int(row[0])
                elif len(row) > 0 and row[0] == 'Value':
                    flag = True
            return nSpots

        def loadPositions(fileName):
            csvFile = csv.reader(open(fileName, 'rU'), delimiter=',')
            x, y, z, t = [], [], [], []
            flag = False
            for row in csvFile:
                if flag == True:
                    x.append(float(row[0]))
                    y.append(float(row[1]))
                    z.append(float(row[2]))
                    t.append(float(row[6]))
                elif len(row) > 0 and row[0] == 'Position X':
                    flag = True
            return np.array([x, y, z, t])

        ftmp = glob.glob('{0}*'.format(self.batchFolder))  # list of everything
        ftmp1 = glob.glob('{0}*.*'.format(self.batchFolder))  # list of files
        ftmp = [t for t in ftmp if t not in ftmp1]  # use only folders
        sort_nicely(ftmp)
        if len(ftmp) < 3:
            self.spotsLoaded = False
            printLog('Msg: Some spots folders are missing: {0}'.format(self.folder))

        # else:
            # self.spotsLoaded = True
            # tmp = glob.glob('{1}{0}*Number_of_Spots_per_Time_Point*'.format(self.splitSym, ftmp[0]))
            # if len(tmp) > 0:
            #     self.green = loadSpots(tmp[0])
            #     self.greenPos = loadPositions(glob.glob('{1}{0}*Position*'.format(self.splitSym, ftmp[0]))[0])
            # else:
            #     printLog('no green spots file')
            #
            # tmp = glob.glob('{1}{0}*Number_of_Spots_per_Time_Point*'.format(self.splitSym, ftmp[1]))
            # if len(tmp) > 0:
            #     self.red = loadSpots(tmp[0])
            #     self.redPos = loadPositions(glob.glob('{1}{0}*Position*'.format(self.splitSym, ftmp[1]))[0])
            # else:
            #     printLog('no red spots file')
            #
            # tmp = glob.glob('{1}{0}*Number_of_Spots_per_Time_Point*'.format(self.splitSym, ftmp[2]))
            # if len(tmp) > 0:
            #     self.yellow = loadSpots(tmp[0])
            #     self.yellowPos = loadPositions(glob.glob('{1}{0}*Position*'.format(self.splitSym, ftmp[2]))[0])
            # else:
            #     printLog('no yellow spots file')


        else:  # made some changes to recognizing file names because in cases where spots were run more than once, the
            # file assignment was incorrect (it was picked based on order, not name-- this affected hundreds of RNAi)
            for i in range(len(ftmp)):
                if "Ch1" in ftmp[i]:
                    green = i
                elif "Ch2" in ftmp[i]:
                    red = i
                elif "Ch3" in ftmp[i]:
                    yellow = i

            self.spotsLoaded = True
            tmp = glob.glob('{1}{0}*Number_of_Spots_per_Time_Point*'.format(self.splitSym, ftmp[green]))
            if len(tmp) > 0:
                self.green = loadSpots(tmp[0])
                self.greenPos = loadPositions(glob.glob('{1}{0}*Position*'.format(self.splitSym, ftmp[green]))[0])
            else:
                printLog('no green spots file')

            tmp = glob.glob('{1}{0}*Number_of_Spots_per_Time_Point*'.format(self.splitSym, ftmp[red]))
            if len(tmp) > 0:
                self.red = loadSpots(tmp[0])
                self.redPos = loadPositions(glob.glob('{1}{0}*Position*'.format(self.splitSym, ftmp[red]))[0])
            else:
                printLog('no red spots file')

            tmp = glob.glob('{1}{0}*Number_of_Spots_per_Time_Point*'.format(self.splitSym, ftmp[yellow]))
            if len(tmp) > 0:
                self.yellow = loadSpots(tmp[0])
                self.yellowPos = loadPositions(glob.glob('{1}{0}*Position*'.format(self.splitSym, ftmp[yellow]))[0])
            else:
                printLog('no yellow spots file')


    def loadAvgCurves(self):
        """
        creates a dictionary of average curves in self.avgCurves. Loads average curves from pickles
        :return:
        """
        self.avgCurves = {}
        if os.path.exists(FILENAME_AVG_CURVES_GLS):
            print('loadAvgCurves: curves loaded')
            with open(FILENAME_AVG_CURVES_GLS, 'rb') as input:
                self.avgCurves = pickle.load(input)
            self.curvesLoaded = True
        else:
            print('loadAvgCurves: NO curves loaded')
            self.curvesLoaded = False

    def setCurveStartPoints(self):
        """
        Assigns the beginning point when each curve is reliable (requires max is at least 20 spots, uses timepoints
        where spots > 10)
        :return:
        """
        indG, indR, indY = [0, 0, 0]
        if max(self.green) > 20: indG = np.where(self.green > 10)[0][0]
        if max(self.red) > 20: indR = np.where(self.red > 10)[0][0]
        if max(self.yellow) > 20: indY = np.where(self.yellow > 10)[0][0]
        for curveName in AVG_CURVES_GLS:
            if curveName in ['spotsG']: self.startPoint[curveName] = indG
            if curveName in ['spotsR']: self.startPoint[curveName] = indR
            if curveName in ['spotsY']: self.startPoint[curveName] = indY
            if curveName in ['MoI0G', 'MoI1G']: self.startPoint[curveName] = indG + 3
            if curveName in ['MoI0R', 'MoI1R']: self.startPoint[curveName] = indR + 3
            if curveName in ['MoI0Y', 'MoI1Y']: self.startPoint[curveName] = indY + 3
            if curveName in ['CoM0R', 'CoM1R']: self.startPoint[curveName] = max(indG, indR) + 3
            if curveName in ['CoM0Y', 'CoM1Y']: self.startPoint[curveName] = max(indG, indY) + 3
        for key in self.startPoint:
            self.startPoint[key] = int(self.startPoint[key])

    def setCMParams(self):
        curveNames = ['CoM{ax}{col}'.format(ax=ax, col=ch) for ax in range(2) for ch in ['R', 'Y']]
        for cn in curveNames:
            self.setCM_MoIParams(cn, 'GLS')

    def setMIParams(self):
        curveNames = ['MoI{0}{1}'.format(ax, c) for c in ['G', 'R', 'Y'] for ax in range(2)]
        for cn in curveNames:
            self.setCM_MoIParams(cn, 'GLS')

    def setSigmParams(self):
        if np.sum(self.green > 20) > 5:  # markers in green are on (5 timepoints with 20 or more spots)
            popt, perr, residsEarly, residsLate, tail, tailErr = self.fitSigmSpots(0)
            self.params['aG'] = popt[0]
            self.params['mG'] = popt[1]
            self.params['sG'] = popt[2]
            self.params['rG'] = popt[3]
            self.params['residsEarlyG'] = residsEarly
            self.params['residsLateG'] = residsLate
            self.params['aGErr'] = perr[0]
            self.params['mGErr'] = perr[1]
            self.params['sGErr'] = perr[2]
            self.params['rGErr'] = perr[3]
            self.params['Gtail'] = tail
            self.params['GtailErr'] = tailErr
        else:
            self.params['aG'] = 0.
            self.params['mG'] = np.nan
            self.params['sG'] = 0.
            self.params['rG'] = np.nan
            self.params['residsEarlyG'] = np.nan
            self.params['residsLateG'] = np.mean(self.avgCurves['spotsG'][1])
            self.params['aGErr'] = 0.
            self.params['mGErr'] = np.nan
            self.params['sGErr'] = 0.
            self.params['rGErr'] = np.nan
            self.params['Gtail'] = 0.
            self.params['GtailErr'] = 0.
        if np.sum(self.red > 20) > 5:
            popt, perr, residsEarly, residsLate, tail, tailErr = self.fitSigmSpots(1)
            self.params['aR'] = popt[0]
            self.params['mR'] = popt[1]
            self.params['sR'] = popt[2]
            self.params['rR'] = popt[3]
            self.params['residsEarlyR'] = residsEarly
            self.params['residsLateR'] = residsLate
            self.params['aRErr'] = perr[0]
            self.params['mRErr'] = perr[1]
            self.params['sRErr'] = perr[2]
            self.params['rRErr'] = perr[3]
            self.params['Rtail'] = tail
            self.params['RtailErr'] = tailErr
        else:
            self.params['aR'] = 0.
            self.params['mR'] = np.nan
            self.params['sR'] = 0.
            self.params['rR'] = np.nan
            self.params['residsEarlyR'] = np.nan
            self.params['residsLateR'] = np.mean(self.avgCurves['spotsR'][1])
            self.params['aRErr'] = 0.
            self.params['mRErr'] = np.nan
            self.params['sRErr'] = 0.
            self.params['rRErr'] = np.nan
            self.params['Rtail'] = 0.
            self.params['RtailErr'] = 0.
        if np.sum(self.yellow > 20) > 5:
            popt, perr, residsEarly, residsLate, tail, tailErr = self.fitSigmSpots(2)
            self.params['aY'] = popt[0]
            self.params['mY'] = popt[1]
            self.params['sY'] = popt[2]
            self.params['rY'] = popt[3]
            self.params['residsEarlyY'] = residsEarly
            self.params['residsLateY'] = residsLate
            self.params['aYErr'] = perr[0]
            self.params['mYErr'] = perr[1]
            self.params['sYErr'] = perr[2]
            self.params['rYErr'] = perr[3]
            self.params['Ytail'] = tail
            self.params['YtailErr'] = tailErr
        else:
            self.params['aY'] = 0.
            self.params['mY'] = np.nan
            self.params['sY'] = 0.
            self.params['rY'] = np.nan
            self.params['residsEarlyY'] = np.nan
            self.params['residsLateY'] = np.mean(self.avgCurves['spotsY'][1])
            self.params['aYErr'] = 0.
            self.params['mYErr'] = np.nan
            self.params['sYErr'] = 0.
            self.params['rYErr'] = np.nan
            self.params['Ytail'] = 0.
            self.params['YtailErr'] = 0.

        if self.checkGoodParam('mY') and self.checkGoodParam('mG'):
            self.params['mYmG'] = self.params['mY'] - self.params['mG']
        if self.checkGoodParam('mR') and self.checkGoodParam('mG'):
            self.params['mRmG'] = self.params['mR'] - self.params['mG']
        if self.checkGoodParam('mY') and self.checkGoodParam('mR'):
            self.params['mYmR'] = self.params['mY'] - self.params['mR']
        aT = 0  # total number of spots
        if not np.isnan(self.params['aG']):
            aT += self.params['aG']
        if not np.isnan(self.params['aR']):
            aT += self.params['aR']
        if not np.isnan(self.params['aY']):
            aT += self.params['aY']
        if aT != 0:
            self.params['fracG'] = 1. * self.params['aG'] / aT
            self.params['fracR'] = 1. * self.params['aR'] / aT
            self.params['fracY'] = 1. * self.params['aY'] / aT
        else:
            self.params['fracG'] = 0
            self.params['fracR'] = 0
            self.params['fracY'] = 0

    def getCurve(self, curveName):
        """
        returns curve data with time normalized by t0 and tScale. Data values are set to nan before the start point.
        :param curveName: str
        :return:
        """
        if curveName == 'tIntG':
            y = np.copy(self.tInt[0])
        elif curveName == 'tIntR':
            y = np.copy(self.tInt[1])
        elif curveName == 'spotsG':
            y = np.copy(self.green)
        elif curveName == 'spotsR':
            y = np.copy(self.red)
        elif curveName == 'spotsY':
            y = np.copy(self.yellow)
        elif curveName == 'CoM0R':
            y = self.getDistFromGreen(1, 0)
        elif curveName == 'CoM1R':
            y = self.getDistFromGreen(1, 1)
        elif curveName == 'CoM0Y':
            y = self.getDistFromGreen(2, 0)
        elif curveName == 'CoM1Y':
            y = self.getDistFromGreen(2, 1)
        elif curveName == 'MoI0G':
            y = self.getAllMoments(0, 0)
        elif curveName == 'MoI1G':
            y = (self.getAllMoments(0, 1) + self.getAllMoments(0, 1)) / 2.
        elif curveName == 'MoI0R':
            y = self.getAllMoments(1, 0)
        elif curveName == 'MoI1R':
            y = (self.getAllMoments(1, 1) + self.getAllMoments(1, 2)) / 2.
        elif curveName == 'MoI0Y':
            y = self.getAllMoments(2, 0)
        elif curveName == 'MoI1Y':
            y = (self.getAllMoments(2, 1) + self.getAllMoments(2, 2)) / 2.
        else:
            return np.nan, np.nan
        y[:self.startPoint[curveName]] = np.nan
        return self.tScale * (self.time - self.t0), y * self.scale[curveName]

    def getCurveNames(self):
        return np.array(AVG_CURVES_GLS)

    def getSpotsPositions(self, ch, slide):
        """
        calculates positions of spots at given slide
        :param ch: channel
        :param slide: slide number (int)
        :return: array (N, 3) where the first index is point number and the second is x, y, z poisition
        """
        if ch == 0:
            points = self.greenPos
        elif ch == 1:
            points = self.redPos
        elif ch == 2:
            points = self.yellowPos
        else:
            printLog('no specified channel')
            return None
        if points is not None:
            x, y, z, t = points
            if slide in t:
                return np.transpose([x[t == slide], y[t == slide], z[t == slide]])
            else:
                None
        else:
            None

    def getAllMoments(self, ch, axis):
        if self.loaded or np.max(self.MoI[ch, axis]) == 0:  # check if there are any values, i.e. it was calculated
            # previously. If yes then use the values, if no then calculate new.
            moments = []
            for slide in range(31):
                points = self.getSpotsPositions(ch, slide)
                if points is not None and points.size > 0:
                    center = getCenterOfMass(points)
                    moment = calcMoment(points, center, axis)
                    moment /= points.shape[0]
                    if axis == 0:
                        moment /= self.radius ** 2
                    else:
                        moment /= self.length ** 2
                    moments.append(moment)
                else:
                    moments.append(np.nan)
            self.MoI[ch, axis] = np.array(moments)
        return self.MoI[ch, axis]

    def getDistFromGreen(self, ch, axis):
        d = np.ones(31) * np.nan
        for slide in range(31):
            points = self.getSpotsPositions(ch, slide)
            pointsG = self.getSpotsPositions(0, slide)
            if points is not None and pointsG is not None and points.shape[0] > 10:
                if np.max(self.CMPos[ch, slide]) > 0:
                    xc, yc, zc = self.CMPos[ch, slide]
                else:
                    xc, yc, zc = getCenterOfMass(points)
                    self.CMPos[ch, slide] = np.array([xc, yc, zc])
                if np.max(self.CMPos[0, slide]) > 0:
                    xe, ye, ze = self.CMPos[0, slide]
                else:
                    xe, ye, ze = getCenterOfMass(pointsG)
                    self.CMPos[0, slide] = np.array([xe, ye, ze])
                    #             xe, ye, ze = self.getEmbCenter()
                if xe is not None:
                    if axis == 0:
                        d[slide] = np.sqrt((yc - ye) ** 2 + (zc - ze) ** 2) / self.radius
                    else:
                        d[slide] = np.sqrt((xc - xe) ** 2) / self.length  # Use rotation independent projection
                        # (np.sqrt((xc-xe)**2+(zc-ze)**2)+np.sqrt((yc-ye)**2+(xc-xe)**2))/2./self.length
                else:
                    printLog('can not get embryo center for slide={0}'.format(slide))
        return d

    def fitSigmSpots(self, ch):
        """
        Fit sigmoidal curve to the spots data
        :param ch: channel
        :return:
        """
        def sigmoidal(x, a, m, s, r):
            """
            sigmoidal curve without the baseline
            :param x: time coordinate
            :param a: plateau
            :param m: mid point
            :param s: slope
            :param r: 1
            :return:
            """
            return a * (1. - (1. + np.exp((x - m) / s)) ** (-r))

        def getErr(y):
            dy = (np.roll(y, -1) - y)[:-1]
            ddy = (np.roll(dy, -1) - dy)[:-1]
            err = np.interp(np.arange(y.size), np.arange(ddy.size) + 1, np.abs(ddy))
            #     err= median_filter(err,3)
            errMin = 10
            err[np.where(err < errMin)] = errMin
            return err

        show = False
        x = self.time
        if ch == 0:
            y = self.green
        elif ch == 1:
            y = self.red
        elif ch == 2:
            y = self.yellow
        dy = (np.roll(self.yellow, -1) - self.yellow)[:-1]
        ind = np.where(dy < -60)[0]  # check that worm didn't hatch. if drop is more then 50 spots, the worm has
        # hatched, then use point before that happens.
        if len(ind) > 0 and ind[0] > 20:
            x = x[:ind[0] + 1]
            y = y[:ind[0] + 1]
        endInd = y.size
        if self.tMove is not None:
            endInd = self.tMove - 1 - GLS_MOVE_PRECISION  # GLS_MOVE_PRECISION is the precision that we add in
            # moveDetect.getMovementGLS(emb); -1 because movement is the next slide after cutoff
        else:  # find a peak in the counts (maximum value)
            maxV = 0
            peak = False
            for val in y:
                if val > maxV:
                    maxV = val
                if val < 0.8 * maxV and maxV > 0.8 * np.max(y):
                    peak = True
                    break
            if peak:
                endInd = np.argmin(np.abs(y - maxV)) + 1

        r0 = 1.  # exponent
        startInd = 0  # max(0,np.where(y>0)[0][0]-1)
        endInd = max(startInd + 1, int(endInd))
        b0 = np.mean(y[endInd - 5:endInd])
        m0 = x[startInd:endInd][np.argmin(np.abs(y[startInd:endInd] - 0.5 * b0))]  # x position of the middle intensity
        s0 = 2  # width
        p0 = (b0, m0, s0, r0)
        params = Parameters()
        params.add('a', p0[0], min=0.5 * b0, max=np.max(y[:endInd]))
        params.add('m', p0[1], min=-5, max=max(x))
        params.add('s', p0[2], min=0, max=max(x))
        params.add('r', p0[3], min=0, max=10, vary=False)
        mod = Model(sigmoidal)
        err = getErr(y)
        #         err=np.ones_like(err)
        if np.max(y) < 20:
            popt = np.ones(4) * np.nan
            perr = np.ones(4) * np.nan
            residsEarly = np.nan
            residsLate = np.nan
            perr = np.nan * np.ones_like(popt)
        else:
            try:
                res = mod.fit(y[:endInd], x=x[:endInd], params=params, weights=1. / err[:endInd] ** 2)
                popt = np.array(
                    [res.best_values['a'], res.best_values['m'], res.best_values['s'], res.best_values['r']])
                if show:
                    fig = myFigure()
                    fig.errorbar(x[:endInd], y[:endInd], err[:endInd])
                    fig.plot(x, sigmoidal(x, *popt))
                    plt.show()
                th = popt[0]
                ind0 = np.ceil(popt[1] + 3)
                residsEarly = np.mean(1. * np.abs(y - sigmoidal(x, *popt))[:endInd])
                if endInd < x.size:
                    residsLate = np.sign(np.mean(1. * (y - sigmoidal(x, *popt))[endInd:])) * np.mean(
                        1. * np.abs(y - sigmoidal(x, *popt))[endInd:])
                else:
                    residsLate = np.nan
            except:
                printLog('{1}: could not fit sigmoidal to channel {0}'.format(ch, self.label))
                popt = np.ones(4) * np.nan
                perr = np.ones(4) * np.nan
                residsEarly = np.nan
                residsLate = np.nan
            try:
                perr = np.sqrt(np.diag(res.covar))
                perr = np.concatenate((perr, np.zeros(popt.size - perr.size)))
            except:
                printLog('{1}: bad parameter errors for fit sigmoidal to channel {0}'.format(ch, self.label))
                perr = np.nan * np.ones_like(popt)
        tail, tailErr = np.nan, np.nan
        if ch == 0:
            if y.size > 5:
                tail = np.array([np.mean(y[-5:])])[0]
                tailErr = np.array([np.std(y[-5:])])[0]
        elif ch == 1:
            if y.size > 5:
                tail = np.array([np.mean(y[-5:])])[0]
                tailErr = np.array([np.std(y[-5:])])[0]
        if ch == 2:
            if y.size > 5:
                tail = np.array([np.mean(y[-7:])])[0]
                tailErr = np.array([np.std(y[-7:])])[0]
        return popt, perr, residsEarly, residsLate, tail, tailErr

    def show(self, curveName, fig=None, setColor=True):
        """
        Shows various curves
        :param curveName: name of the curve to show
        :param fig: myFigure object to plot the curve on
        :param setColor: use channel color (0 green(GFP), 1 red(RFP)) if True
        :return: myFigure object
        """


        if curveName == 'spotsG': return self.showSigmFit(fig, [0], setColor)
        if curveName == 'spotsR': return self.showSigmFit(fig, [1], setColor)
        if curveName == 'spotsY': return self.showSigmFit(fig, [2], setColor)
        if curveName == 'CoM0R': return self.showDistCentMass(0, fig, [1], setColor)
        if curveName == 'CoM1R': return self.showDistCentMass(1, fig, [1], setColor)
        if curveName == 'CoM0Y': return self.showDistCentMass(0, fig, [2], setColor)
        if curveName == 'CoM1Y': return self.showDistCentMass(1, fig, [2], setColor)
        if curveName == 'MoI0G': return self.showMOI(0, fig, [0], setColor)
        if curveName == 'MoI1G': return self.showMOI(1, fig, [0], setColor)
        if curveName == 'MoI0R': return self.showMOI(0, fig, [1], setColor)
        if curveName == 'MoI1R': return self.showMOI(1, fig, [1], setColor)
        if curveName == 'MoI0Y': return self.showMOI(0, fig, [2], setColor)
        if curveName == 'MoI1Y': return self.showMOI(1, fig, [2], setColor)

    def showSigmFit(self, fig=None, chs=[0, 1, 2], setColor=True):
        def sigmoidal(x, a, m, s, r):
            return a * (1. - (1. + np.exp((x - m) / s)) ** (-r))
        showControls = True
        

        if fig is None: fig = myFigure()
        for ch in chs:
            curveName = 'spots{col}'.format(col=['G', 'R', 'Y'][ch])
            x, y = self.getCurve(curveName)
            if ch == 0:
                keys = ['aG', 'mG', 'sG', 'rG']
                color = 'g'
            elif ch == 1:
                keys = ['aR', 'mR', 'sR', 'rR']
                color = 'r'
            elif ch == 2:
                keys = ['aY', 'mY', 'sY', 'rY']
                color = 'y'
            a, m, s, r = [self.params[k] for k in keys]
            if setColor:
                fig.plot(x, y, color=color)  # data
                fig.plot(x[:self.cutoff[curveName]], y[:self.cutoff[curveName]], 'k--')  # data until td (dashed)
                # fig.plot(x[:self.cutoff[curveName]], y[:self.cutoff[curveName]], 'lightgray')  # data until td (dashed)
                fig.plot(x, sigmoidal(x / self.tScale, a, m / self.tScale, s, r) * self.scale[curveName], color=color)
                # sigmoidal fit
            if showControls:
                if curveName in self.avgCurves:
                    x, y, yerr = self.avgCurves[curveName]
                    fig.errorbar(x[::10], y[::10], yerr[::10], color=color)  # average curve
            else:  # else plot only data
                fig.plot(x, y, label=self.label)
        return fig

    def showMaxProj(self, t):
        """
        returns figure with maximum projection with centers of mass indicated with dots
        :param t: time point to show (normalized developmental time)
        :return: myFigure object
        """
        self.loadAllSpots()
        slide = int(np.round(t / self.tScale + self.t0))
        if slide <= 31:
            fig = myFigure()
            fig.markerSize = 150
            if self.image is None:
                self.loadImages()
            g = np.array(self.image.getMaxProject(slide, 0, 0), dtype=np.float32) / 20000.  # dividing by 20000 to
            # normalize to a 0-1 range. Empirically determined 20000 to optimize display.
            r = np.array(self.image.getMaxProject(slide, 1, 0), dtype=np.float32) / 30000.
            xG, yG, zG = getCenterOfMass(self.getSpotsPositions(0, slide))
            xR, yR, zR = getCenterOfMass(self.getSpotsPositions(1, slide))
            xY, yY, zY = getCenterOfMass(self.getSpotsPositions(2, slide))
            r[np.where(r > 1)] = 1.  # sets the values that were above 20000, before normalization, to a max value of 1
            g[np.where(g > 1)] = 1.
            b = np.zeros_like(g)
            fig.imshow(np.dstack([r, g, b]), colorbar=False)
            fig.scatter([xG], [yG], color='darkgray')  # sets a dark bkgd for CoM spot (next step) to make more visible
            fig.scatter([xR], [yR], color='darkgray')
            fig.scatter([xY], [yY], color='darkgray')
            fig.markerSize = 100
            fig.scatter([xG], [yG], color='springgreen')
            fig.scatter([xR], [yR], color='red')
            fig.scatter([xY], [yY], color='gold')
            fig.noAxis()
            fig.setBGColor('k')
            fig.title('{2} t={0}, d={1}'.format((slide - self.t0)*self.tScale,
                                                np.sqrt((xR - xG) ** 2 + (yR - yG) ** 2 + (zR - zG) ** 2) / self.length,
                                                self.label))
            return fig

    def showMOI(self, axis, fig=None, chs=range(3),
                setColor=True):
        """
        plots MOI for GLS. If set color arg is False, uses autocolor
        (would use this if plotting lots of embryos on one plot)
        :param axis: axis (0-end-on, 1-ap)
        :param fig: myFigure object
        :param chs: channels
        :param setColor: use channel color (0 green(GFP), 1 red(RFP)) if True
        :return: myFigure object
        """
        if not fig: fig = myFigure()  # if no figure exists, creates figure
        color = ['g', 'r', 'y']  # specifies colors if setColor is True
        for ch in chs:  # iterates through channels
            curveName = 'MoI{0}{1}'.format(axis, ['G', 'R', 'Y'][ch])
            fig.title(curveName)
            t, MOI = self.getCurve(curveName)
            if setColor:
                fig.plot(t, MOI, label=self.label, color=color[
                    ch])  # plots in green, red... use for one embryo, or population if dont need to ID individual embs
                fig.plot(t[:self.cutoff[curveName]], MOI[:self.cutoff[curveName]], 'k--', label=self.label)
                if curveName in self.avgCurves:
                    x, y, yerr = self.avgCurves[curveName]
                    fig.errorbar(x[::10], y[::10], yerr[::10], color=color[ch])
            else:
                fig.plot(t, MOI,
                         label=self.label)  # plots in multiple colors...
                # use for multiple embryos to trace back individual embs
        return fig

def get_avg_counts_pilot(folders, show=False):
    '''
    gets the spots data and calculates the average spots (3 consecutive timepoints) either just prior to movement
    (if moves) or the last 3 timepoints (for non-movement). This is used to get the data needed for figures in 50-gene
    pilot paper.
    :param folders: examples :    # folders = [FOLDER_IN + 'cropped/EMBD0000/GLS/20140515T131046/Emb1/']
    # folders = [FOLDER_IN + 'cropped/EMBD0141/GLS/Emb2/']
    :return:
    '''
    for folder in folders[:]:
        # emb = MSEmbryo(folder, check_version=False)
        emb = GSEmbryo(folder)
        red = emb.red
        yellow = emb.yellow
        green = emb.green
        print(emb.params['tMove'], emb.tMove, emb.movement)
        if show:
            print('green = {0}').format(green)
            print('red = {0}').format(red)
            print('yellow = {0}').format(yellow)
        if emb.movement:
            greenAv = np.mean(green[emb.tMove-7:emb.tMove-4])
            redAv = np.mean(red[emb.tMove -7:emb.tMove - 4])
            yellAv = np.mean(yellow[emb.tMove - 7:emb.tMove - 4])
            print(greenAv, redAv, yellAv)
        elif not emb.movement:
            greenAv = np.mean(green[-3:])
            redAv = np.mean(red[-3:])
            yellAv = np.mean(yellow[-3:])
            print('Green = {0}, Red = {1}, Yellow = {2}'.format(greenAv, redAv, yellAv))

def show_frames_control_ims(self, strain):
    if strain =='GLS':
        folders = [FOLDER_IN + 'cropped/EMBD0088/GLS/Emb5']
    for folder in folders[:]:
        emb = GSEmbryo(folder, check_version=False)
        for i in range(18):
            fig = emb.showMaxProj(emb.tMove-i)
            # fig.show()
            fig.save(
                'Z:/control_timepoints/GLS/{e}_{t}.tif'.format(e=emb.label,
                                                                                                            t=i))
    if strain =='MS':
        folders = [FOLDER_IN + 'cropped/EMBD0088/MS/Emb5']
        for folder in folders[:]:
            emb = MSEmbryo(folder, check_version=False)
            for i in range(18):
                slide = int(np.round(self.time / self.tScale + self.t0))
                if slide <= 31:
                    fig = myFigure()
                    fig.markerSize = 150
                    if self.image is None:
                        self.loadImages()
                    g = np.array(self.image.getMaxProject(slide, 0, 0), dtype=np.float32) / 20000.  # dividing by 20000 to
                    # normalize to a 0-1 range. Empirically determined 20000 to optimize display.
                    r = np.array(self.image.getMaxProject(slide, 1, 0), dtype=np.float32) / 30000.
                    # xG, yG, zG = getCenterOfMass(self.getSpotsPositions(0, slide))
                    r[np.where(r > 1)] = 1.  # sets the values that were above 20000, before normalization, to a max value of 1
                    g[np.where(g > 1)] = 1.
                    b = np.zeros_like(g)
                    fig.imshow(np.dstack([r, g, b]), colorbar=False)
                    fig.noAxis()
                    fig.setBGColor('k')
                    fig.title('{1} t={0}'.format((slide - self.t0) * self.tScale,
                                                        self.label))
                    fig.show()





if __name__ == '__main__':
    # folders = [FOLDER_IN + 'cropped/EMBD0000/MS/20140516T130559/Emb4']
    # folders = [FOLDER_IN + 'cropped/EMBD0000/MS/20140515T131046/Emb4/']
    # folders = [FOLDER_IN + 'cropped/EMBD0000/MS/20140430T140422/Emb2/']
    # folders = [FOLDER_IN + 'cropped/EMBD0000/MS/20140514T130911/Emb5/']
    # folders = [FOLDER_IN + 'cropped/EMBD0000/MS/20140402T140154/Emb2/']
    # folders = [FOLDER_IN + 'cropped/EMBD0000/MS/20140402T140154/Emb2/']
    # folders = [FOLDER_IN + 'cropped/EMBD0000/GLS/20140501T135409/Emb5/']
    # folders = [FOLDER_IN + 'cropped/EMBD0000/GLS/20140416T140401/Emb4/']
    # folders = [FOLDER_IN + 'cropped/EMBD0000/GLS/20140515T131046/Emb4/']
    # folders = [FOLDER_IN + 'cropped/EMBD0141/MS/Emb2/']
    # folders = [FOLDER_IN + 'cropped/EMBD0005/MS/Emb5/']
    folders = [FOLDER_IN + 'cropped/EMBD0002/MS/Emb7/']
    # folders = [FOLDER_IN + 'cropped/EMBD1360/GLS/Emb3']
    # folders = [FOLDER_IN + 'cropped/EMBD0088/GLS/Emb5']


    # get_avg_counts_pilot(folders)  #, show=True

    # for folder in folders[:]:
    #     # emb = GSEmbryo(folder, check_version=False)
    #     emb = MSEmbryo(folder, check_version=False)
    #     x,y = emb.getCurve('headInt')
    #     # print x,y
    #     # emb.updateParams()
    #     # curveName = 'spotsR'
    #     # x, y = emb.getCurve(curveName)
    #     # print(x, y)
    #     # emb.showSigmFit().show()
    #
    #     # for i in range(18):
    #     #     fig = emb.showMaxProj(emb.tMove-i)
    #     #     # fig.show()
    #     #     fig.save(
    #     #         'Z:/GLS/{e}_{t}.tif'.format(e=emb.label,
    #     #                                                                                                  t=i))
    #     emb.fix_head_params_error()

    fix_headInt_params_by_rnai_number()
    # fix_headInt_params_control_embs()








    #
    #     emb.showHead().show()

    #     print(emb.label, emb.RNAi)
    #     emb.showSigmFit().show()

        # params0 = copy(emb.params)
        # emb.load_from_pickle()
        # emb.refresh()
        # for key in params0:
        #     if params0[key] != emb.params[key] and not np.isnan(params0[key])*np.isnan(emb.params[key]):
        #         print('key =', key, params0[key], emb.params[key])
        #         pass
        # emb.updateParams()
        # emb.refresh()

        # emb.scale['tIntR'] = 1
        # emb.scale['tIntG'] = 1
        # emb.scale['lengthR'] = 1
        # emb.scale['MoI0R'] = 1
        # # emb.scale['headInt'] = 1
        # print('t0=', emb.t0)
        # print(emb.thresh[1])
        # print('G=', emb.params['aG'])
        # print('R=', emb.params['aR'])
        # print('Y=', emb.params['aY'])
        # curveName = 'spotsR'
        # x, y = emb.getCurve(curveName)
        # print(x, y)
        # emb.showSigmFit().show()
        # emb.show('tIntR')
        # emb.show('tIntG')
        # emb.show('headInt')
        # emb.show('MoI1R').show()
        # emb.show('MoI0R').show()
        # emb.show('MoI0G').show()
        # emb.show('MoI1G').show()
        # emb.show('lengthR').show()


        # folders = [FOLDER_IN + 'cropped/EMBD0052/MS/Emb5/']
        # for folder in folders[:]:
        #     emb = MSEmbryo(folder, check_version=False)
        # for i in range(5,15):
        #     fig = emb.showMaxProj(i)
        #     fig.show()
