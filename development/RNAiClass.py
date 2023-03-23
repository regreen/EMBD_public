"""
Created on Jul 8, 2016

@author: Becky
"""

from varLookup import FOLDER_IN, PARAM_NAMES, printLog, DT0_GLS_2_MS
from fileLookup import FOLDER_RNAI_DATA, FOLDER_TIME_ALIGN_MOVIES
from emb_handler import loadEmbs
from myFigure import myFigure
import pickle
from myMath import AvgCurve
from myFunc import saveImagesMulti
import cv2
import os
import numpy as np
from cv2 import FONT_ITALIC
from PIL import ImageFont, ImageDraw, Image
from copy import copy
import db_utils_embryos
import db_utils_genes
import emb_handler
from embdFunc import getGeneName, get_param_names, split_params, getNDDistance, getNDAvgVector, getGeneCommonName, \
    get_rnai_revision_version, make_dict_float
from params_use_assign import get_params_use

PARAM_NAMES_USE_GLS, PARAM_NAMES_USE_MS = get_params_use()


class RNAiClass(object):
    """
    Class that takes all embryos from a given RNAi in both GS and MS (rules and attributes that apply to RNAi condition)
    """

    def __init__(self, RNAi, verb=True, check_version=False):
        '''
        Constructor
        '''
        self.RNAi = RNAi  # EMBD number... i.e. 132
        self.verb = verb
        conn, curs = db_utils_embryos.initialize_db()
        self.id = db_utils_genes.get_id_from_embd(self.RNAi, curs)
        self.check_version = check_version
        if check_version:
            self.__version__ = get_rnai_revision_version()
        else:
            self.__version__ = None
        self.GLSfolder = FOLDER_IN + 'cropped/EMBD{0:04}/GLS/'.format(self.RNAi)
        self.MSfolder = FOLDER_IN + 'cropped/EMBD{0:04}/MS/'.format(self.RNAi)
        self.fileName = FOLDER_RNAI_DATA + 'EMBD{0:04}'.format(self.RNAi)
        self.paramsGLS = {}  # average, good parameters GLS
        self.paramsGLSstd = {}  # std, good parameters GLS
        self.paramsMS = {}  # average, good parameters MS
        self.paramsMSstd = {}  # std, good parameters MS
        self.paramsEmbsMS = {}
        self.paramsEmbsGLS = {}
        for pN in PARAM_NAMES:  # adds parameter and value, for each pN, to params dictionary
            self.paramsGLS[pN] = np.nan
            self.paramsGLSstd[pN] = np.nan
            self.paramsMS[pN] = np.nan
            self.paramsMSstd[pN] = np.nan
            self.paramsEmbsMS[pN] = np.nan
            self.paramsEmbsGLS[pN] = np.nan
        self.paramsUseGLS = []
        self.paramsUseMS = []
        self.GLSUseEmbs = []
        self.MSUseEmbs = []
        self.GLSMoveEmbs, self.GLSNoMoveEmbs = [], []
        self.MSMoveEmbs, self.MSNoMoveEmbs = [], []
        self.prev_GLS, self.prev_MS = 0, 0
        #         self.posND = np.zeros(len(PARAM_NAMES_USE_MS)+len(PARAM_NAMES_USE_GLS))
        self.posND = np.ones(2 * len(PARAM_NAMES)) * np.nan
        self.stdV = np.zeros(2 * len(PARAM_NAMES))
        self.pNamesUse = np.array([PARAM_NAMES_USE_GLS, PARAM_NAMES_USE_MS])
        self.origin = None
        self.paramNorms = None
        self.weights = None
        self.NDFlag = False
        self.label = 'EMBD{0:04}'.format(self.RNAi)
        self.refreshed = False
        self.load_params_data()
        self.setOrigin()

    def save_params_data(self):  # makes PICKLE
        conn, curs = db_utils_embryos.initialize_db()
        columns = {}
        for key in self.paramsGLS:
            columns[key + '_GLS'] = self.paramsGLS[key]
        for key in self.paramsMS:
            columns[key + '_MS'] = self.paramsMS[key]
        columns['prev_GLS'] = self.prev_GLS
        columns['prev_MS'] = self.prev_MS
        columns['version'] = self.__version__
        db_utils_embryos.update_row(self.id, columns, curs, 'genes')
        conn.commit()
        conn.close()

        # with open(self.fileName,
        #           'wb') as output:  # this means open file (filename) as the file to save to (output) and do XYZ with it and then close it
        #     pickle.dump(self.__version__, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsGLS, output,
        #                 pickle.HIGHEST_PROTOCOL)  # saves ave GLS params in first arg to fileName (output) leave 3rd arg alone.
        #     pickle.dump(self.paramsGLSstd, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsMS, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsMSstd, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsUseMS, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsUseGLS, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsEmbsMS, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsEmbsGLS, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.prev, output, pickle.HIGHEST_PROTOCOL)

    def load_params_data(self):
        # try:
        if self.verb:
            printLog('loading from MySQL...')
        conn, curs = db_utils_embryos.initialize_db()
        p_names = get_param_names() + ['prev_GLS', 'prev_MS', 'version']
        params = db_utils_embryos.get_row(self.id, curs, 'genes', p_names)
        self.paramsGLS, self.paramsMS, other = split_params(params)
        self.paramsUseGLS = [key for key in self.paramsGLS if not np.isnan(self.paramsGLS[key])]
        self.paramsUseMS = [key for key in self.paramsMS if not np.isnan(self.paramsMS[key])]
        self.prev_GLS = other['prev_GLS']
        self.prev_MS = other['prev_MS']
        ver = other['version']
        if ver is None:
            self.load_from_pickle()
            self.save_params_data()
        if self.check_version and self.__version__ != ver:
            printLog('old version, updating RNA')
            self.refresh_params_data()
        else:
            self.__version__ = ver
            # except Exception as e:
            #     printLog('Failed to load from mysql: {}'.format(e))
            #     self.refresh_params_data()

    def load_from_pickle(self):
        if not os.path.exists(self.fileName):
            self.refresh_params_data()
        else:
            try:
                if self.verb:
                    printLog('loading from {0}...'.format(self.fileName))
                with open(self.fileName, 'rb') as input:
                    self.__version__ = pickle.load(input)
                    self.paramsGLS = pickle.load(input)
                    self.paramsGLSstd = pickle.load(input)
                    self.paramsMS = pickle.load(input)
                    self.paramsMSstd = pickle.load(input)
                    self.paramsUseMS = pickle.load(input)
                    self.paramsUseGLS = pickle.load(input)
                    self.paramsEmbsMS = pickle.load(input)
                    self.paramsEmbsGLS = pickle.load(input)
                    self.prev_GLS = pickle.load(input)
                    if self.check_version and self.__version__ != get_rnai_revision_version():
                        printLog('old pickle vesion {0}...'.format(self.fileName))
                        self.refresh_params_data()
            except:
                printLog('Failed to load from pickle')
                self.refresh_params_data()

    def refresh_params_data(self):
        self.__version__ = get_rnai_revision_version()
        self.paramsGLS = {}  # average, good parameters GLS
        self.paramsGLSstd = {}  # std, good parameters GLS
        self.paramsMS = {}  # average, good parameters MS
        self.paramsMSstd = {}  # std, good parameters MS
        self.paramsEmbsMS = {}
        self.paramsEmbsGLS = {}
        for pN in PARAM_NAMES:
            self.paramsGLS[pN] = np.nan
            self.paramsGLSstd[pN] = np.nan
            self.paramsMS[pN] = np.nan
            self.paramsMSstd[pN] = np.nan
            self.paramsEmbsMS[pN] = np.nan
            self.paramsEmbsGLS[pN] = np.nan
        self.paramsUseGLS = []
        self.paramsUseMS = []
        self.getEmbryos()
        self.prev_GLS, self.prev_MS = self.getPrev()
        self.asignUse()
        self.calcAverageParams()
        self.paramsGLS = make_dict_float(self.paramsGLS)
        self.paramsGLSstd = make_dict_float(self.paramsGLSstd)
        self.paramsMS = make_dict_float(self.paramsMS)
        self.paramsMSstd = make_dict_float(self.paramsMSstd)
        self.save_params_data()
        self.refreshed = True

    def refresh_embs(self):
        self.getEmbryos()
        emb_handler.refreshEmbsMulti(self.GLSMoveEmbs + self.GLSNoMoveEmbs + self.MSMoveEmbs + self.MSNoMoveEmbs)

    def asignUse(self):
        if len(self.GLSNoMoveEmbs) > 1:
            self.GLSUseEmbs = self.GLSNoMoveEmbs
        else:
            self.GLSUseEmbs = self.GLSMoveEmbs + self.GLSNoMoveEmbs
        if len(self.MSNoMoveEmbs) > 1:
            self.MSUseEmbs = self.MSNoMoveEmbs
        else:
            self.MSUseEmbs = self.MSMoveEmbs + self.MSNoMoveEmbs

    def getEmbryos(self):
        embs = loadEmbs(self.GLSfolder, inds=None)
        # for emb in embs:
        #     emb.refresh()
        # embs = embdFunc.refreshEmbsMulti(embs)
        self.GLSMoveEmbs, self.GLSNoMoveEmbs = self.splitByMovement(
            embs)  # generates two lists- moving and not moving GLS embs
        embs = loadEmbs(self.MSfolder, inds=None)
        for emb in embs:
            # emb.refresh()
            emb.image = None
        # embs = embdFunc.refreshEmbsMulti(embs)
        self.MSMoveEmbs, self.MSNoMoveEmbs = self.splitByMovement(
            embs)  # generates two lists- moving and not moving MS embs

    def splitByMovement(self,
                        embs):  # checks if moving- if so, checks to see if it continues to move (maintains high late resids)... embs that move and maintain high resides are called "wt"
        mEmbs = []  # empty list of moving embryo objects
        nomEmbs = []  # empty list of non-moving (WT-like) embryo objects
        for emb in embs:  # evaluates for each embryo whether it moves and continues to move, or doesnt move/stops moving and appends embryo to appropriate list
            if emb.params['movement'] and not emb.dead: mEmbs.append(emb)
            if not emb.params['movement'] and not emb.dead:
                nomEmbs.append(emb)
            elif emb.dead:
                printLog('DEAD {0} embryo {1}'.format(self.label, emb.folder))
        return mEmbs, nomEmbs  # returns a list of moving and non-moving embs

    def calcAverageParams(self):
        """stores ave value for params which have small errors (less than 50% of parameter's value) in self.params"""

        for pN in PARAM_NAMES:  # check whether error value is present for individual param, if not assigns pNErr as pN
            goodParamsGLS = []  # empty list to add good GLS parameters to
            for emb in self.GLSUseEmbs:  # iterates through non-moving embryos objects in GLSUseEmbs list
                if emb.checkGoodParam(pN): goodParamsGLS.append(emb.params[pN])  # Appends value for each good parameter
            if len(goodParamsGLS) > 0:  # checks to be sure there are values in the list
                self.paramsUseGLS.append(pN)  # adds parameter names to paramsUseGLS list
                self.paramsGLS[pN] = np.median(goodParamsGLS)  # adds median good parameters to dictionary self.params
                self.paramsGLSstd[pN] = np.std(goodParamsGLS)  # adds std for RNAi to dictionary
                self.paramsEmbsGLS[pN] = goodParamsGLS
            # if pN == 'mG':
            #                     printLog('paramsGLS',self.paramsGLS[pN])

            goodParamsMS = []  # empty list to add good MS params to
            for emb in self.MSUseEmbs:
                if emb.checkGoodParam(pN): goodParamsMS.append(emb.params[pN])  # Appends value for each good parameter
            if len(goodParamsMS) > 0:  # checks to be sure there are values in the list
                self.paramsUseMS.append(pN)  # adds parameter names to paramsUseMS list
                self.paramsMS[pN] = np.median(goodParamsMS)  # adds median good parameters to dictionary self.params
                self.paramsMSstd[pN] = np.std(goodParamsMS)  # adds std for RNAi to dictionary
                self.paramsEmbsMS[pN] = goodParamsMS

    def setOrigin(self, origin=None, norms=None, weights=None):
        self.controlOrigin = {'bYErr_GLS': np.nan, 'YtailErr_GLS': 17.490743823529414,
                              'CoM1Rscale_MS_MS': 0.96693014683544309, 'maxG_MS': 758262979.74683547,
                              'tScale_MS': 1.0417721518987342, 'MoI0Rtd_MS_GLS': np.nan, 'RtailErr_MS': np.nan,
                              'MoI0Rscale_GLS_MS': np.nan, 'sG_MS': 2.2476907088607598, 'mG_GLS': 1.9215470171568629,
                              'bY_GLS': np.nan, 'bGErr_MS': 17081041.645569619, 'sGErr_MS': 0.27105817139240507,
                              'MoI0Ytd_GLS_GLS': 16.974855392156861, 'MoI1GstdTail_GLS_GLS': 0.06371408799019608,
                              'MoI0YstdTail_GLS_MS': np.nan, 'CoM1Rscale_GLS_GLS': 0.84664579166666665,
                              'sR_GLS': 1.4972319019607845, 'mG_MS': 4.0411543147208118, 'MoI1GstdTail_GLS_MS': np.nan,
                              'MoI1Gtd_GLS_MS': np.nan, 'CoM0Yscale_GLS_GLS': 0.93003753186274518,
                              'MoI0RavgRes_MS_MS': -0.019298576536708863, 'MoI0GstdTail_GLS_GLS': 0.060636927696078438,
                              'rYErr_GLS': 0.0, 'movement_MS': 0.96962025316455691, 'MoI0Rscale_MS_GLS': np.nan,
                              'MoI1Gscale_GLS_GLS': 0.97759795343137257, 'MoI1Rtd_MS_MS': 11.939508860759494,
                              'CoM1Rtd_GLS_GLS': 13.086759803921568, 'CoM0Rscale_GLS_MS': np.nan,
                              'MoI1Ytd_GLS_MS': np.nan, 'CoM1GstdTail_MS_MS': 0.059301411215189864,
                              'aGErr_GLS': 1.9975246648241209, 'mR_GLS': 5.175359558823529,
                              'CoM0Rscale_GLS_GLS': 0.90412982352941185, 'MoI1Ytd_GLS_GLS': 16.947639705882356,
                              'aRErr_GLS': 3.6984995359828012, 'mRErr_MS': 0.41311080888324875,
                              'MoI1Gtd_MS_MS': 11.972296202531645, 'aSigHead_GLS': np.nan,
                              'CoM0Gscale_MS_MS': 0.96263296962025302, 'MoI1Rscale_MS_GLS': np.nan,
                              'CoM1RstdTail_MS_GLS': np.nan, 'MoI1Rscale_GLS_GLS': 0.98146218137254893,
                              'maxR_GLS': 1679324583.3333333, 'CoM1Yscale_GLS_MS': np.nan,
                              'tailHead_MS': 1.5565271980529, 'tailLength_GLS': np.nan, #'tailHead_MS': 0.16143561620253163
                              'mYErr_GLS': 0.18583931348039215, 'bR_MS': 698311174.68354428, 'MoI1Rtd_MS_GLS': np.nan,
                              'CoM0Gscale_MS_GLS': np.nan, 'MoI0Ytd_GLS_MS': np.nan, 'scaleHead_MS': 1.0186421739130436,
                              'scaleLength_GLS': np.nan, 'Gtail_GLS': 77.895588235294127, 'bR_GLS': np.nan,
                              'rY_GLS': 1.0, 'CoM0Ytd_GLS_MS': np.nan, 'MoI1YstdTail_GLS_MS': np.nan,
                              'sYErr_MS': np.nan, 'tailHead_GLS': np.nan, 'Ytail_MS': np.nan,
                              'tMove_GLS': 18.210443627450982, 'MoI0Gtd_MS_MS': 11.973878481012658,
                              'aR_GLS': 109.07239362745098, 'CoM0Rtd_MS_GLS': np.nan, 'mYmG_GLS': 3.7506534528678306,
                              'MoI0RstdTail_GLS_MS': np.nan, 'MoI0GavgRes_GLS_GLS': -0.01533301942132353,
                              'CoM1YstdTail_GLS_MS': np.nan, 'MoI1GstdTail_MS_MS': 0.038012747468354431,
                              'fracY_MS': np.nan, 'fracG_GLS': 0.31976287745098042, 'mY_MS': np.nan,
                              'bGErr_GLS': np.nan, 'sY_GLS': 1.3668969950980392, 'MoI0GavgRes_GLS_MS': np.nan,
                              'MoI0Gtd_MS_GLS': np.nan, 'sSigHead_GLS': np.nan, 'CoM0GstdTail_MS_GLS': np.nan,
                              'bY_MS': np.nan, 'aYErr_MS': np.nan, 'CoM0RstdTail_MS_GLS': np.nan,
                              'MoI1RavgRes_GLS_MS': np.nan, 'sYErr_GLS': 0.16300924240196082,
                              'MoI0YavgRes_GLS_MS': np.nan, 'CoM0GstdTail_MS_MS': 0.077416272632911384,
                              'sGaussHead_GLS': np.nan, 'MoI0GavgRes_MS_GLS': np.nan,
                              'CoM1RstdTail_GLS_GLS': 0.043678223039215688, 'CoM0Gtd_MS_GLS': np.nan,
                              'fracY_GLS': 0.29092612254901962, 'mRmG_MS': -2.2533147079795395,
                              'CoM0YavgRes_GLS_GLS': -0.0091808152573529395, 'CoM1Rtd_GLS_MS': np.nan,
                              'aGaussHead_GLS': np.nan, 'CoM1GavgRes_MS_MS': -0.010589138903797471,
                              'CoM0RavgRes_GLS_MS': np.nan, 'residsLateG_MS': np.nan,
                              'movement_GLS': 0.93627450980392157, 'RtailErr_GLS': 19.72861911764706,
                              'CoM0YavgRes_GLS_MS': np.nan, 'YtailErr_MS': np.nan, 'MoI0GstdTail_GLS_MS': np.nan,
                              'CoM0Yscale_GLS_MS': np.nan, 'CoM1Gscale_MS_MS': 0.81602606835443048,
                              'MoI0RstdTail_MS_GLS': np.nan, 'MoI0RavgRes_GLS_MS': np.nan, 'aY_MS': np.nan,
                              'residsEarlyG_GLS': 4.4176435294117642, 'CoM1Rscale_GLS_MS': np.nan,
                              'mR_MS': 1.8125412450632914, 'mSigHead_GLS': np.nan, 'sRErr_GLS': 0.22093295724815723,
                              'MoI0Gscale_GLS_GLS': 0.98544978431372554, 'mGaussHead_GLS': np.nan, 'rG_GLS': 1.0,
                              'scaleLength_MS': 0.99433048354430376, 'MoI1GavgRes_MS_MS': -0.0068756340848101271,
                              'rRErr_GLS': 0.0, 'residsLateR_GLS': 53.464367156862735, 'rGErr_GLS': 0.0,
                              'devHead_GLS': np.nan, 'CoM1Gtd_MS_MS': 11.896779746835444,
                              'aRErr_MS': 41290218.581218272, 'aSigHead_MS': 1.9631509081218275,
                              'MoI1GstdTail_MS_GLS': np.nan, 'CoM1RavgRes_MS_MS': 3.0156893670888658e-05,
                              'CoM1Rtd_MS_MS': 11.819582278481015, 'sRErr_MS': 0.32665483502538067,
                              'residsEarlyG_MS': np.nan, 'rR_GLS': 1.0, 'aGErr_MS': 22247523.06329114,
                              'fracR_MS': np.nan, 'MoI1YavgRes_GLS_GLS': -0.0046244549024509803,
                              'MoI0Gscale_MS_GLS': np.nan, 'bG_MS': 646225751.89873421, 'mYmR_GLS': 0.50657209013350124,
                              'maxHead_MS': 2.1002177215189874, 'MoI0Rtd_GLS_GLS': 16.9829681372549,
                              'CoM1YavgRes_GLS_MS': np.nan, 'MoI0YstdTail_GLS_GLS': 0.035586934313725491,
                              'devHead_MS': 20.2912972305131, 'rR_MS': 1.0, 'CoM0RavgRes_GLS_GLS': -0.011830159460784313, #'devHead_MS': 2.3868
                              'MoI1Gscale_GLS_MS': np.nan, 'residsLateG_GLS': 2.2989975735294124,
                              'MoI1Rscale_MS_MS': 0.92609765316455683, 'MoI1Gtd_MS_GLS': np.nan,
                              'MoI1GavgRes_MS_GLS': np.nan, 'GtailErr_MS': np.nan, 'residsLateR_MS': np.nan,
                              'scaleHead_GLS': np.nan, 'CoM0RstdTail_GLS_MS': np.nan, 'CoM1Rtd_MS_GLS': np.nan,
                              'mGErr_GLS': 178.88216248492461, 'MoI1Rtd_GLS_GLS': 16.924737745098039,
                              'CoM1RavgRes_GLS_MS': np.nan, 'CoM1RstdTail_GLS_MS': np.nan, 'MoI1Rscale_GLS_MS': np.nan,
                              'maxG_GLS': 1349349889.7058823, 'aR_MS': 610724681.01265824, 'CoM0GavgRes_MS_GLS': np.nan,
                              'devLength_GLS': np.nan, 'MoI0YavgRes_GLS_GLS': -0.0059381861264705882,
                              'fracR_GLS': 0.38931102205882345, 'MoI1Rtd_GLS_MS': np.nan, 'tMove_MS': 12.16695918367347,
                              'CoM0Gtd_MS_MS': 11.81680506329114, 'CoM1GstdTail_MS_GLS': np.nan,
                              'Ytail_GLS': 147.92437181372549, 'bYErr_MS': np.nan,
                              'MoI0RstdTail_GLS_GLS': 0.047997642647058823, 'aY_GLS': 81.190458578431375,
                              'maxHead_GLS': np.nan, 'bRErr_MS': 33511061.649746194, 'rGErr_MS': 0.0,
                              'aYErr_GLS': 2.353977497782108, 'MoI1RstdTail_GLS_MS': np.nan, 'rG_MS': 1.0,
                              'CoM1Gtd_MS_GLS': np.nan, 'CoM1RavgRes_GLS_GLS': -0.012550435264215686,
                              'MoI1Gscale_MS_GLS': np.nan, 'CoM1YstdTail_GLS_GLS': 0.067946717156862735,
                              'residsEarlyY_GLS': 3.0268587867647057, 'residsLateY_MS': np.nan,
                              'MoI0GstdTail_MS_GLS': np.nan, 'CoM0Rscale_MS_MS': 1.0028619924050632,
                              'CoM0Rscale_MS_GLS': np.nan, 'MoI0GavgRes_MS_MS': -0.0073416472784810132,
                              'MoI0Yscale_GLS_MS': np.nan, 'mYmG_MS': np.nan, 'CoM1Yscale_GLS_GLS': 0.99358843627450988,
                              'mY_GLS': 5.6709817941176475, 'bG_GLS': np.nan, 'bRErr_GLS': np.nan,
                              'mSigHead_MS': 1.8601007999999999, 'CoM0RavgRes_MS_MS': 0.0093940502784810136,
                              'MoI1Yscale_GLS_GLS': 0.98771595098039222, 'MoI1GavgRes_GLS_GLS': -0.0062532987303921574,
                              'CoM1YavgRes_GLS_GLS': 0.00061024391495098098, 'mRErr_GLS': 0.22262055356265356,
                              'MoI0Rscale_MS_MS': 0.95406086835443038, 'MoI1RavgRes_MS_GLS': np.nan,
                              'CoM0YstdTail_GLS_GLS': 0.067493157598039227, 'aGaussHead_MS': 0.079761520304568526,
                              'mYmR_MS': np.nan, 'MoI0RstdTail_MS_MS': 0.050182506050632916,
                              'MoI0Rtd_MS_MS': 11.971068354430381, 'CoM0YstdTail_GLS_MS': np.nan, 'rY_MS': np.nan,
                              'CoM1Ytd_GLS_GLS': 16.249458333333333, 'aG_GLS': 89.344855392156859,
                              'CoM0Rtd_GLS_MS': np.nan, 'MoI1RstdTail_MS_GLS': np.nan,
                              'GtailErr_GLS': 13.211739093137254, 'mYErr_MS': np.nan, 'CoM1GavgRes_MS_GLS': np.nan,
                              'sGaussHead_MS': 3.9504039999999998, 'maxR_MS': 782318506.32911396,
                              'MoI0Gtd_GLS_GLS': 16.963723039215687, 'MoI1Gscale_MS_MS': 0.98882014177215205,
                              'sY_MS': np.nan, 'MoI1RstdTail_GLS_GLS': 0.028166627794117646,
                              'sSigHead_MS': 6.0331294444444454, 'CoM1Gscale_MS_GLS': np.nan,
                              'MoI1RavgRes_GLS_GLS': -0.0062760544156862741, 'MoI0Rtd_GLS_MS': np.nan,
                              'devLength_MS': 21.528567088607595, 'CoM0Ytd_GLS_GLS': 15.827249999999999,
                              'CoM0Rtd_GLS_GLS': 16.136029411764707, 'sR_MS': 2.5961415791189872,
                              'residsLateY_GLS': 66.775018651960778, 'MoI1Yscale_GLS_MS': np.nan,
                              'mRmG_GLS': 3.2407050253164558, 'mGaussHead_MS': 4.6074755555555553,
                              'CoM0Rtd_MS_MS': 11.718741772151899, 'residsEarlyR_MS': np.nan,
                              'MoI1YavgRes_GLS_MS': np.nan, 'tScale_GLS': 0.9948039215686274,
                              'CoM1RavgRes_MS_GLS': np.nan, 'MoI1GavgRes_GLS_MS': np.nan, 'MoI0Gtd_GLS_MS': np.nan,
                              'rYErr_MS': np.nan, 'residsEarlyR_GLS': 4.6629248529411766,
                              'MoI0Rscale_GLS_GLS': 0.99926214460784313, 'CoM0RstdTail_MS_MS': 0.15934806481012659,
                              'CoM1Rscale_MS_GLS': np.nan, 'MoI0Yscale_GLS_GLS': 1.0017748774509805,
                              'MoI1RavgRes_MS_MS': -0.011872208444556963, 'MoI1Gtd_GLS_GLS': 16.940968137254902,
                              'sGErr_GLS': 126.741469159799, 'rRErr_MS': np.nan, 'CoM0RavgRes_MS_GLS': np.nan,
                              'CoM0GavgRes_MS_MS': 0.013371890751898733, 'tailLength_MS': 192.49283921518989,
                              'CoM0RstdTail_GLS_GLS': 0.059845477450980393, 'MoI1RstdTail_MS_MS': 0.04254413318987342,
                              'fracG_MS': np.nan, 'Rtail_MS': np.nan, 'Rtail_GLS': 156.13186274509803,
                              'sG_GLS': 1.3158710054232845, 'aG_MS': 435264169.62025315,
                              'mGErr_MS': 0.29364657868020305, 'Gtail_MS': np.nan, 'residsEarlyY_MS': np.nan,
                              'MoI0RavgRes_MS_GLS': np.nan, 'MoI0RavgRes_GLS_GLS': -0.0070226766053921568,
                              'CoM1RstdTail_MS_MS': 0.17862180556962029, 'MoI1YstdTail_GLS_GLS': 0.025877258455882353,
                              'MoI0GstdTail_MS_MS': 0.039432833291139245, 'CoM1Ytd_GLS_MS': np.nan,
                              'MoI0Gscale_MS_MS': 0.98549217468354433, 'MoI0Gscale_GLS_MS': np.nan}
        self.controlStd = {'bYErr_GLS': np.nan, 'YtailErr_GLS': 5.7480885693648913,
                           'CoM1Rscale_MS_MS': 0.21534782912941827, 'maxG_MS': 156920541.55648118,
                           'tScale_MS': 0.13736728456969979, 'MoI0Rtd_MS_GLS': np.nan, 'RtailErr_MS': np.nan,
                           'MoI0Rscale_GLS_MS': np.nan, 'sG_MS': 0.4881897522070735, 'mG_GLS': 0.9645900583460707,
                           'bY_GLS': np.nan, 'bGErr_MS': 10422077.384296101, 'sGErr_MS': 0.11681285921733457,
                           'MoI0Ytd_GLS_GLS': 0.40598924955623944, 'MoI1GstdTail_GLS_GLS': 0.017964034354585372,
                           'MoI0YstdTail_GLS_MS': np.nan, 'CoM1Rscale_GLS_GLS': 0.29793624144386877,
                           'sR_GLS': 0.38488216432648248, 'mG_MS': 0.88170158178403513, 'MoI1GstdTail_GLS_MS': np.nan,
                           'MoI1Gtd_GLS_MS': np.nan, 'CoM0Yscale_GLS_GLS': 0.26656620287025445,
                           'MoI0RavgRes_MS_MS': 0.037455726796340787, 'MoI0GstdTail_GLS_GLS': 0.020698343674253085,
                           'rYErr_GLS': 0.0, 'movement_MS': 0.17162988614357763, 'MoI0Rscale_MS_GLS': np.nan,
                           'MoI1Gscale_GLS_GLS': 0.10732398055857115, 'MoI1Rtd_MS_MS': 0.39103142649122902,
                           'CoM1Rtd_GLS_GLS': 2.9488949660401973, 'CoM0Rscale_GLS_MS': np.nan, 'MoI1Ytd_GLS_MS': np.nan,
                           'CoM1GstdTail_MS_MS': 0.060500904708226444, 'aGErr_GLS': 1.4444028359567229,
                           'mR_GLS': 0.73940100742107162, 'CoM0Rscale_GLS_GLS': 0.17129375971537797,
                           'MoI1Ytd_GLS_GLS': 0.52627626392652604, 'aRErr_GLS': 3.7912058501909227,
                           'mRErr_MS': 1.9477542407066499, 'MoI1Gtd_MS_MS': 0.27561139558803499, 'aSigHead_GLS': np.nan,
                           'CoM0Gscale_MS_MS': 0.3973660867771257, 'MoI1Rscale_MS_GLS': np.nan,
                           'CoM1RstdTail_MS_GLS': np.nan, 'MoI1Rscale_GLS_GLS': 0.088094409973587737,
                           'maxR_GLS': 181214765.34283227, 'CoM1Yscale_GLS_MS': np.nan,
                           'tailHead_MS': 0.333737195341383, 'tailLength_GLS': np.nan, #'tailHead_MS': 0.49501328410853707
                           'mYErr_GLS': 0.088343345171148946, 'bR_MS': 110871068.01922636, 'MoI1Rtd_MS_GLS': np.nan,
                           'CoM0Gscale_MS_GLS': np.nan, 'MoI0Ytd_GLS_MS': np.nan, 'scaleHead_MS': 0.26085270289040635,
                           'scaleLength_GLS': np.nan, 'Gtail_GLS': 16.481474452978269, 'bR_GLS': np.nan, 'rY_GLS': 0.0,
                           'CoM0Ytd_GLS_MS': np.nan, 'MoI1YstdTail_GLS_MS': np.nan, 'sYErr_MS': np.nan,
                           'tailHead_GLS': np.nan, 'Ytail_MS': np.nan, 'tMove_GLS': 1.1167242149753591,
                           'MoI0Gtd_MS_MS': 0.27396693955167822, 'aR_GLS': 12.55081996986997, 'CoM0Rtd_MS_GLS': np.nan,
                           'mYmG_GLS': 0.81532531270393749, 'MoI0RstdTail_GLS_MS': np.nan,
                           'MoI0GavgRes_GLS_GLS': 0.046622891737998634, 'CoM1YstdTail_GLS_MS': np.nan,
                           'MoI1GstdTail_MS_MS': 0.025582317458318304, 'fracY_MS': np.nan,
                           'fracG_GLS': 0.02434005253097191, 'mY_MS': np.nan, 'bGErr_GLS': np.nan,
                           'sY_GLS': 0.31177207692373321, 'MoI0GavgRes_GLS_MS': np.nan, 'MoI0Gtd_MS_GLS': np.nan,
                           'sSigHead_GLS': np.nan, 'CoM0GstdTail_MS_GLS': np.nan, 'bY_MS': np.nan, 'aYErr_MS': np.nan,
                           'CoM0RstdTail_MS_GLS': np.nan, 'MoI1RavgRes_GLS_MS': np.nan, 'sYErr_GLS': 0.1264221689166215,
                           'MoI0YavgRes_GLS_MS': np.nan, 'CoM0GstdTail_MS_MS': 0.12968685659481313,
                           'sGaussHead_GLS': np.nan, 'MoI0GavgRes_MS_GLS': np.nan,
                           'CoM1RstdTail_GLS_GLS': 0.021576564177957031, 'CoM0Gtd_MS_GLS': np.nan,
                           'fracY_GLS': 0.021980258947110712, 'mRmG_MS': 1.2079530156208624,
                           'CoM0YavgRes_GLS_GLS': 0.046291354185156258, 'CoM1Rtd_GLS_MS': np.nan,
                           'aGaussHead_GLS': np.nan, 'CoM1GavgRes_MS_MS': 0.08920657643141322,
                           'CoM0RavgRes_GLS_MS': np.nan, 'residsLateG_MS': np.nan, 'movement_GLS': 0.24426328437845085,
                           'RtailErr_GLS': 8.6790378659274019, 'CoM0YavgRes_GLS_MS': np.nan, 'YtailErr_MS': np.nan,
                           'MoI0GstdTail_GLS_MS': np.nan, 'CoM0Yscale_GLS_MS': np.nan,
                           'CoM1Gscale_MS_MS': 0.32466174379201784, 'MoI0RstdTail_MS_GLS': np.nan,
                           'MoI0RavgRes_GLS_MS': np.nan, 'aY_MS': np.nan, 'residsEarlyG_GLS': 1.3441842429774706,
                           'CoM1Rscale_GLS_MS': np.nan, 'mR_MS': 1.2035508050538217, 'mSigHead_GLS': np.nan,
                           'sRErr_GLS': 0.43975814713912981, 'MoI0Gscale_GLS_GLS': 0.1220883147223865,
                           'mGaussHead_GLS': np.nan, 'rG_GLS': 0.0, 'scaleLength_MS': 0.079787132428234417,
                           'MoI1GavgRes_MS_MS': 0.043455430415501123, 'rRErr_GLS': 0.0,
                           'residsLateR_GLS': 16.769497953338067, 'rGErr_GLS': 0.0, 'devHead_GLS': np.nan,
                           'CoM1Gtd_MS_MS': 0.52734104891264666, 'aRErr_MS': 130734786.61666673,
                           'aSigHead_MS': 0.2371551437728916, 'MoI1GstdTail_MS_GLS': np.nan,
                           'CoM1RavgRes_MS_MS': 0.14567759991378806, 'CoM1Rtd_MS_MS': 0.8547360005124528,
                           'sRErr_MS': 0.95864935432436593, 'residsEarlyG_MS': np.nan, 'rR_GLS': 0.0,
                           'aGErr_MS': 10148490.936777456, 'fracR_MS': np.nan,
                           'MoI1YavgRes_GLS_GLS': 0.016008124328654578, 'MoI0Gscale_MS_GLS': np.nan,
                           'bG_MS': 141137900.6802552, 'mYmR_GLS': 0.63598463546225792,
                           'maxHead_MS': 0.35534712312581884, 'MoI0Rtd_GLS_GLS': 0.24166016177410587,
                           'CoM1YavgRes_GLS_MS': np.nan, 'MoI0YstdTail_GLS_GLS': 0.022988165929896453,
                           'devHead_MS': 4.59321836832564, 'rR_MS': 0.0, 'CoM0RavgRes_GLS_GLS': 0.030407593282412819, #'devHead_MS': 6.8673793577713687
                           'MoI1Gscale_GLS_MS': np.nan, 'residsLateG_GLS': 20.249268826288294,
                           'MoI1Rscale_MS_MS': 0.14434198110003285, 'MoI1Gtd_MS_GLS': np.nan,
                           'MoI1GavgRes_MS_GLS': np.nan, 'GtailErr_MS': np.nan, 'residsLateR_MS': np.nan,
                           'scaleHead_GLS': np.nan, 'CoM0RstdTail_GLS_MS': np.nan, 'CoM1Rtd_MS_GLS': np.nan,
                           'mGErr_GLS': 2661.5604041133693, 'MoI1Rtd_GLS_GLS': 0.69899351826415934,
                           'CoM1RavgRes_GLS_MS': np.nan, 'CoM1RstdTail_GLS_MS': np.nan, 'MoI1Rscale_GLS_MS': np.nan,
                           'maxG_GLS': 184546982.59325346, 'aR_MS': 114123901.69850771, 'CoM0GavgRes_MS_GLS': np.nan,
                           'devLength_GLS': np.nan, 'MoI0YavgRes_GLS_GLS': 0.02353256417872613,
                           'fracR_GLS': 0.028608907050934735, 'MoI1Rtd_GLS_MS': np.nan, 'tMove_MS': 2.2324142894092303,
                           'CoM0Gtd_MS_MS': 0.98406270510028981, 'CoM1GstdTail_MS_GLS': np.nan,
                           'Ytail_GLS': 20.6230779225703, 'bYErr_MS': np.nan,
                           'MoI0RstdTail_GLS_GLS': 0.024685408058469121, 'aY_GLS': 6.7519340010245301,
                           'maxHead_GLS': np.nan, 'bRErr_MS': 149378290.78571212, 'rGErr_MS': 0.0,
                           'aYErr_GLS': 1.5538702741432096, 'MoI1RstdTail_GLS_MS': np.nan, 'rG_MS': 0.0,
                           'CoM1Gtd_MS_GLS': np.nan, 'CoM1RavgRes_GLS_GLS': 0.026989896893836351,
                           'MoI1Gscale_MS_GLS': np.nan, 'CoM1YstdTail_GLS_GLS': 0.026503766184760805,
                           'residsEarlyY_GLS': 0.92491148648232402, 'residsLateY_MS': np.nan,
                           'MoI0GstdTail_MS_GLS': np.nan, 'CoM0Rscale_MS_MS': 0.30616512922506367,
                           'CoM0Rscale_MS_GLS': np.nan, 'MoI0GavgRes_MS_MS': 0.048014605500989474,
                           'MoI0Yscale_GLS_MS': np.nan, 'mYmG_MS': np.nan, 'CoM1Yscale_GLS_GLS': 0.17382459869419206,
                           'mY_GLS': 0.81367193840034402, 'bG_GLS': np.nan, 'bRErr_GLS': np.nan,
                           'mSigHead_MS': 4.9024230132443627, 'CoM0RavgRes_MS_MS': 0.18027275065692661,
                           'MoI1Yscale_GLS_GLS': 0.077812141564227125, 'MoI1GavgRes_GLS_GLS': 0.034494641868752084,
                           'CoM1YavgRes_GLS_GLS': 0.031809157156158237, 'mRErr_GLS': 0.17754624905607153,
                           'MoI0Rscale_MS_MS': 0.12810844540965624, 'MoI1RavgRes_MS_GLS': np.nan,
                           'CoM0YstdTail_GLS_GLS': 0.024890512558360132, 'aGaussHead_MS': 0.28071666225918057,
                           'mYmR_MS': np.nan, 'MoI0RstdTail_MS_MS': 0.025794880817173814,
                           'MoI0Rtd_MS_MS': 0.27932493285147669, 'CoM0YstdTail_GLS_MS': np.nan, 'rY_MS': np.nan,
                           'CoM1Ytd_GLS_GLS': 1.6456164180399733, 'aG_GLS': 8.4202533348969659,
                           'CoM0Rtd_GLS_MS': np.nan, 'MoI1RstdTail_MS_GLS': np.nan, 'GtailErr_GLS': 6.4594125561361633,
                           'mYErr_MS': np.nan, 'CoM1GavgRes_MS_GLS': np.nan, 'sGaussHead_MS': 3.2261624020976445,
                           'maxR_MS': 108825576.14431038, 'MoI0Gtd_GLS_GLS': 0.40934400444980074,
                           'MoI1Gscale_MS_MS': 0.09654380260667754, 'sY_MS': np.nan,
                           'MoI1RstdTail_GLS_GLS': 0.015791620569142475, 'sSigHead_MS': 6.2626621189444975,
                           'CoM1Gscale_MS_GLS': np.nan, 'MoI1RavgRes_GLS_GLS': 0.017892133524464823,
                           'MoI0Rtd_GLS_MS': np.nan, 'devLength_MS': 3.8636836935072472,
                           'CoM0Ytd_GLS_GLS': 1.7953732697397715, 'CoM0Rtd_GLS_GLS': 1.9244563706835134,
                           'sR_MS': 0.62628520427723078, 'residsLateY_GLS': 12.101716397482607,
                           'MoI1Yscale_GLS_MS': np.nan, 'mRmG_GLS': 0.82061308949762357,
                           'mGaussHead_MS': 5.5375609403733996, 'CoM0Rtd_MS_MS': 1.2468343242197639,
                           'residsEarlyR_MS': np.nan, 'MoI1YavgRes_GLS_MS': np.nan, 'tScale_GLS': 0.098930361644995191,
                           'CoM1RavgRes_MS_GLS': np.nan, 'MoI1GavgRes_GLS_MS': np.nan, 'MoI0Gtd_GLS_MS': np.nan,
                           'rYErr_MS': np.nan, 'residsEarlyR_GLS': 1.3170938499672478,
                           'MoI0Rscale_GLS_GLS': 0.11738755268794418, 'CoM0RstdTail_MS_MS': 0.09597025690116151,
                           'CoM1Rscale_MS_GLS': np.nan, 'MoI0Yscale_GLS_GLS': 0.1090636163823291,
                           'MoI1RavgRes_MS_MS': 0.02957702124041043, 'MoI1Gtd_GLS_GLS': 0.54494329302434441,
                           'sGErr_GLS': 1783.7292779284614, 'rRErr_MS': np.nan, 'CoM0RavgRes_MS_GLS': np.nan,
                           'CoM0GavgRes_MS_MS': 0.16860942605741391, 'tailLength_MS': 27.641686566035059,
                           'CoM0RstdTail_GLS_GLS': 0.02122608469904912, 'MoI1RstdTail_MS_MS': 0.018958129358575024,
                           'fracG_MS': np.nan, 'Rtail_MS': np.nan, 'Rtail_GLS': 29.531053188339449,
                           'sG_GLS': 0.49605203297762485, 'aG_MS': 112857175.45921478, 'mGErr_MS': 0.13749815347876937,
                           'Gtail_MS': np.nan, 'residsEarlyY_MS': np.nan, 'MoI0RavgRes_MS_GLS': np.nan,
                           'MoI0RavgRes_GLS_GLS': 0.029993843670286676, 'CoM1RstdTail_MS_MS': 0.064883010838061611,
                           'MoI1YstdTail_GLS_GLS': 0.013718286031999007, 'MoI0GstdTail_MS_MS': 0.030305489849142025,
                           'CoM1Ytd_GLS_MS': np.nan, 'MoI0Gscale_MS_MS': 0.12938743031381575, 'MoI0Gscale_GLS_MS': np.nan}
        self.rnaiMean = {'bYErr_GLS': np.nan, 'MoI1Geval_7_GLS_GLS': 0.11848254756035149,
                         'residsEarlyY_GLS': 4.2341044107221109, 'bY_GLS': np.nan,
                         'MoI0RstdTail_MS_MS': 0.040663594822456343, 'MoI1Gtd_GLS_MS': np.nan,
                         'MoI0Rscale_MS_GLS': np.nan, 'mRErr_GLS': 1.3322684516071213,
                         'CoM0YstdTail_GLS_GLS': 0.057288361332520736, 'MoI1Reval_15_MS_GLS': np.nan,
                         'MoI1Rscale_GLS_GLS': 0.98791091586665902, 'MoI1Yeval_22_GLS_MS': np.nan,
                         'scaleHead_MS': 0.98310319677491254, 'bR_GLS': np.nan, 'rY_GLS': 1.0, 'tailHead_GLS': np.nan,
                         'MoI0RavgRes_MS_MS': 0.003798052376218618, 'tMove_GLS': 17.519716456153802,
                         'MoI0Gtd_MS_MS': 11.541416691012405, 'MoI1Gtd_MS_GLS': np.nan, 'MoI0Gtd_MS_GLS': np.nan,
                         'sYErr_GLS': 0.16486663263426368, 'MoI1Yeval_10_GLS_MS': np.nan, 'MoI1Reval_3_MS_GLS': np.nan,
                         'MoI1Geval_1_MS_GLS': np.nan, 'CoM0RavgRes_GLS_MS': np.nan, 'MoI1Rtd_MS_GLS': np.nan,
                         'CoM0Ytd_GLS_MS': np.nan, 'MoI1RstdTail_GLS_GLS': 0.019960459955118057,
                         'CoM1Reval_22_GLS_GLS': 0.062016536636750573, 'MoI1Yeval_7_GLS_GLS': 0.14747770343531885,
                         'MoI0Gscale_GLS_GLS': 0.97257222545197486, 'rRErr_GLS': 0.0, 'MoI1Reval_7_GLS_MS': np.nan,
                         'CoM1Geval_15_MS_MS': 0.090247493420280489, 'CoM1RavgRes_MS_MS': -0.033831516368099904,
                         'sRErr_MS': 0.27255329992478317, 'MoI0Gscale_MS_GLS': np.nan, 'MoI1Gscale_GLS_MS': np.nan,
                         'MoI1Geval_3_MS_GLS': np.nan, 'GtailErr_MS': np.nan, 'CoM1Rtd_MS_GLS': np.nan,
                         'CoM0GavgRes_MS_GLS': np.nan, 'CoM1GstdTail_MS_GLS': np.nan, 'MoI0Yscale_GLS_MS': np.nan,
                         'aY_GLS': 80.997275646151962, 'aYErr_GLS': 2.0730166996064017, 'MoI1RstdTail_GLS_MS': np.nan,
                         'rG_MS': np.nan, 'MoI0RavgRes_MS_GLS': np.nan, 'CoM1Geval_1_MS_MS': 0.0, 'bY_MS': np.nan,
                         'MoI1GstdTail_GLS_GLS': 0.048171308495962024, 'MoI0GavgRes_MS_MS': 0.017456303632165559,
                         'MoI1GavgRes_MS_GLS': np.nan, 'mSigHead_MS': 1.5719212191318181,
                         'CoM1YavgRes_GLS_GLS': 0.016754582939501671, 'CoM1Geval_15_MS_GLS': np.nan,
                         'CoM1Geval_3_MS_MS': 0.09164573968561901, 'CoM1Reval_3_MS_GLS': np.nan, 'rY_MS': np.nan,
                         'CoM1Yeval_10_GLS_GLS': 0.25112234249828386, 'CoM1Reval_7_MS_GLS': np.nan,
                         'GtailErr_GLS': 8.544150132061338, 'CoM0Rtd_GLS_MS': np.nan, 'aG_GLS': 82.883493662022488,
                         'MoI1Geval_10_GLS_GLS': 0.15050928865254934, 'mGaussHead_MS': 4.8707030003504572,
                         'MoI0Gtd_GLS_GLS': 15.656252461624346, 'MoI1RavgRes_GLS_GLS': -0.017127697138771342,
                         'MoI1Yscale_GLS_MS': np.nan, 'mRmG_GLS': 3.1325446087053859, 'CoM1Reval_22_GLS_MS': np.nan,
                         'CoM0Rtd_MS_MS': 7.5678238360403114, 'MoI1Geval_7_GLS_MS': np.nan,
                         'CoM0RavgRes_MS_GLS': np.nan, 'MoI1Rtd_GLS_MS': np.nan, 'fracG_MS': np.nan, 'Rtail_MS': np.nan,
                         'aG_MS': 374457612.18355203, 'MoI0GavgRes_MS_GLS': np.nan, 'residsEarlyY_MS': np.nan,
                         'CoM1Rscale_MS_MS': 1.0616295235532929, 'maxG_MS': 644044963.95032918,
                         'MoI0Rscale_GLS_MS': np.nan, 'sG_MS': 2.4378943824536652, 'MoI1Yeval_14_GLS_MS': np.nan,
                         'MoI0Rtd_GLS_GLS': 16.117403552870137, 'movement_MS': 0.72310756972111556,
                         'CoM1Rtd_GLS_GLS': 11.891775196635912, 'CoM0Rscale_GLS_MS': np.nan, 'MoI1Ytd_GLS_MS': np.nan,
                         'mRErr_MS': 0.47874209165381965, 'CoM1Gscale_MS_MS': 1.0773769735651011,
                         'MoI0Ytd_GLS_MS': np.nan, 'maxR_GLS': 1520448949.27, 'tailHead_MS': 1.1974548330712793,
                         'RtailErr_MS': np.nan, 'bR_MS': 652569363.40276909,
                         'MoI0YstdTail_GLS_GLS': 0.024237087911371771, 'fracY_MS': np.nan, 'CoM0Rscale_MS_GLS': np.nan,
                         'bGErr_GLS': np.nan, 'sSigHead_GLS': np.nan, 'CoM1Yeval_10_GLS_MS': np.nan,
                         'CoM1Geval_1_MS_GLS': np.nan, 'CoM0Gtd_MS_GLS': np.nan, 'CoM0Rscale_MS_MS': 1.1133331686013683,
                         'CoM1Yeval_22_GLS_MS': np.nan, 'MoI0Rscale_GLS_GLS': 1.0173531007792325,
                         'movement_GLS': 0.61299999999999999, 'CoM0Rtd_MS_GLS': np.nan, 'YtailErr_MS': np.nan,
                         'CoM0Yscale_GLS_MS': np.nan, 'CoM1Rscale_GLS_MS': np.nan,
                         'MoI1Reval_3_MS_MS': 0.081392654520657096, 'MoI1Reval_7_MS_GLS': np.nan, 'rGErr_GLS': 0.0,
                         'MoI0RavgRes_GLS_MS': np.nan, 'bG_MS': 561281842.79112184,
                         'CoM1Reval_10_GLS_GLS': 0.074167986144949793, 'devHead_MS': 14.475366780286052,
                         'CoM1Geval_3_MS_GLS': np.nan, 'residsLateR_MS': np.nan, 'CoM0RstdTail_GLS_MS': np.nan,
                         'mGErr_GLS': 9674.7295726784359, 'CoM1RstdTail_GLS_MS': np.nan, 'maxG_GLS': 1221864348.2490001,
                         'aR_MS': 567686857.67064536, 'MoI0YavgRes_GLS_GLS': -0.0074429715939374501,
                         'CoM0RstdTail_GLS_GLS': 0.051091936938660459, 'tMove_MS': 12.337154624662887,
                         'rYErr_MS': np.nan, 'CoM1GavgRes_MS_GLS': np.nan, 'MoI1Gscale_MS_GLS': np.nan,
                         'MoI1Reval_1_MS_MS': 0.0, 'CoM1YstdTail_GLS_GLS': 0.057511488939971274,
                         'MoI1Geval_14_GLS_MS': np.nan, 'residsLateY_MS': np.nan,
                         'CoM1Reval_7_MS_MS': 0.35316439241493691, 'mYmG_MS': np.nan, 'mY_GLS': 3.5573740956125084,
                         'MoI0GavgRes_GLS_MS': np.nan, 'MoI1Rscale_GLS_MS': np.nan, 'MoI1Reval_10_GLS_MS': np.nan,
                         'sGaussHead_MS': 3.4314956638065213, 'maxR_MS': 712540755.07840121, 'Gtail_MS': np.nan,
                         'sY_MS': np.nan, 'CoM1Gscale_MS_GLS': np.nan, 'MoI0Rtd_MS_MS': 11.157473065062517,
                         'mYmG_GLS': 3.5688820050737489, 'MoI1YavgRes_GLS_MS': np.nan, 'MoI1GavgRes_GLS_MS': np.nan,
                         'CoM1Rscale_MS_GLS': np.nan, 'MoI1RavgRes_MS_MS': -0.0077515702005317884,
                         'sGErr_GLS': 5315.7161028440887, 'Rtail_GLS': 137.95010020040081, 'sG_GLS': 1.4308109436045646,
                         'CoM1RstdTail_MS_MS': 0.14442330992136362, 'MoI1YstdTail_GLS_GLS': 0.018391669601371452,
                         'MoI0GstdTail_MS_MS': 0.032200733266257844, 'tScale_MS': 0.99911728336716554,
                         'MoI1Geval_1_MS_MS': 0.0, 'mG_GLS': -0.15314749403976621, 'bGErr_MS': 14500900.533183785,
                         'sGErr_MS': 5.1540717428128131, 'sR_GLS': 1.6954278567099046, 'CoM1Geval_7_MS_GLS': np.nan,
                         'MoI1Gscale_GLS_GLS': 0.98803784450945042, 'MoI1Rtd_MS_MS': 9.1618805343078833,
                         'mR_GLS': 3.0448841769218675, 'CoM0Rscale_GLS_GLS': 0.84210922953591472,
                         'MoI1Ytd_GLS_GLS': 15.790160614019424, 'CoM0Gscale_MS_MS': 1.142648366346291,
                         'CoM1RstdTail_MS_GLS': np.nan, 'residsLateG_MS': np.nan, 'CoM1Yscale_GLS_MS': np.nan,
                         'MoI0Gscale_GLS_MS': np.nan, 'CoM1Geval_7_MS_MS': 0.076833757389910759,
                         'CoM0Gscale_MS_GLS': np.nan, 'MoI1Geval_3_MS_MS': 0.26233653871158014,
                         'CoM1Reval_14_GLS_MS': np.nan, 'scaleLength_GLS': np.nan, 'Gtail_GLS': 71.470340681362714,
                         'MoI1YstdTail_GLS_MS': np.nan, 'sYErr_MS': np.nan, 'Ytail_MS': np.nan,
                         'MoI1Reval_14_GLS_MS': np.nan, 'aR_GLS': 108.96598326187959, 'MoI1Geval_7_MS_GLS': np.nan,
                         'CoM1YstdTail_GLS_MS': np.nan, 'MoI1GstdTail_MS_MS': 0.03006624352876662,
                         'fracG_GLS': 0.29684063258305876, 'mY_MS': np.nan, 'CoM0YstdTail_GLS_MS': np.nan,
                         'sY_GLS': 1.5716073972864846, 'CoM1Reval_14_GLS_GLS': 0.051167769205195772,
                         'MoI0YavgRes_GLS_MS': np.nan, 'fracY_GLS': 0.29489363259761014, 'CoM1Rtd_GLS_MS': np.nan,
                         'RtailErr_GLS': 13.602661752104453, 'CoM1Reval_7_GLS_GLS': 0.074187975227683681,
                         'aY_MS': np.nan, 'mR_MS': 1.6190216508854736, 'MoI1Reval_7_GLS_GLS': 0.16174012747092675,
                         'sRErr_GLS': 0.25791129763697007, 'mGaussHead_GLS': np.nan, 'devHead_GLS': np.nan,
                         'CoM1Gtd_MS_MS': 8.0078013909738139, 'CoM1Reval_10_GLS_MS': np.nan, 'mSigHead_GLS': np.nan,
                         'CoM0RstdTail_MS_GLS': np.nan, 'aGErr_MS': 23133412.289029256, 'fracR_MS': np.nan,
                         'MoI1YavgRes_GLS_GLS': -0.021179023485202921, 'mYmR_GLS': 0.42062668921160634,
                         'maxHead_MS': 1.6935765863653833, 'rR_MS': np.nan,
                         'CoM0RavgRes_GLS_GLS': 0.0029208141568767623, 'CoM1Yeval_7_GLS_MS': np.nan,
                         'scaleHead_GLS': np.nan, 'fracR_GLS': 0.39549283351524012,
                         'MoI1Reval_7_MS_MS': 0.11143903702511544, 'Ytail_GLS': 123.06657142857144,
                         'MoI1Geval_22_GLS_MS': np.nan, 'MoI0RstdTail_GLS_GLS': 0.037263901064747515,
                         'maxHead_GLS': np.nan, 'bRErr_MS': 35828272.665544622, 'rGErr_MS': 0.0,
                         'CoM1RavgRes_GLS_GLS': -0.0049455951580268933, 'MoI1Rscale_MS_GLS': np.nan,
                         'CoM1Yscale_GLS_GLS': 0.95031767290611879, 'bRErr_GLS': np.nan,
                         'CoM0RavgRes_MS_MS': -0.064035664126380759, 'MoI1Yscale_GLS_GLS': 0.9957512635461756,
                         'CoM0YavgRes_GLS_MS': np.nan, 'MoI1Yeval_7_GLS_MS': np.nan,
                         'MoI0Rscale_MS_MS': 0.96002757841642061, 'aGaussHead_MS': 0.47105900561102404,
                         'mYmR_MS': np.nan, 'CoM1Reval_3_MS_MS': 0.35343240273861337, 'CoM1Reval_1_MS_MS': 0.0,
                         'CoM1Ytd_GLS_GLS': 13.462541539296893, 'CoM1Reval_15_MS_MS': 0.34179017033137715,
                         'MoI1RstdTail_MS_GLS': np.nan, 'MoI1GavgRes_GLS_GLS': -0.030397242406518692,
                         'MoI1Geval_7_MS_MS': 0.29038788744512983, 'sSigHead_MS': 3.8470296262439221,
                         'MoI0Ytd_GLS_GLS': 16.24729674606899, 'MoI0YstdTail_GLS_MS': np.nan,
                         'MoI0RstdTail_MS_GLS': np.nan, 'MoI1Reval_15_MS_MS': 0.12159360389318613,
                         'CoM1RstdTail_GLS_GLS': 0.036734034475911885, 'CoM1RavgRes_MS_GLS': np.nan, 'rR_GLS': 1.0,
                         'MoI0Gtd_GLS_MS': np.nan, 'scaleLength_MS': 0.99710457179414957,
                         'residsEarlyR_GLS': 5.5415650242059842, 'MoI0Yscale_GLS_GLS': 1.0165083235051728,
                         'CoM0GavgRes_MS_MS': -0.068342497106730599, 'MoI1Reval_1_MS_GLS': np.nan,
                         'MoI1Geval_15_MS_MS': 0.27842724537392233, 'MoI1RstdTail_MS_MS': 0.034560079131718789,
                         'MoI1Reval_22_GLS_MS': np.nan, 'CoM1Ytd_GLS_MS': np.nan, 'MoI0Gscale_MS_MS': 1.069248951185322,
                         'YtailErr_GLS': 12.033746729345083, 'MoI0Rtd_MS_GLS': np.nan,
                         'CoM1Yeval_22_GLS_GLS': 0.10928107090720258, 'MoI1Yeval_10_GLS_GLS': 0.19233481258259472,
                         'CoM1Rscale_GLS_GLS': 0.79211809252104526, 'mG_MS': 4.1945861374594466,
                         'CoM0Ytd_GLS_GLS': 13.715495835487877, 'CoM0Yscale_GLS_GLS': 0.8932984850565574,
                         'MoI0GstdTail_GLS_GLS': 0.046080766567529623, 'rYErr_GLS': 0.0,
                         'CoM1GstdTail_MS_MS': 0.054054218412318319, 'aGErr_GLS': 4.4376695936598258,
                         'MoI1Geval_22_GLS_GLS': 0.22183626029912831, 'aRErr_GLS': 6.1842020713737034,
                         'MoI1Gtd_MS_MS': 11.573942933662986, 'aSigHead_GLS': np.nan,
                         'CoM1Yeval_14_GLS_GLS': 0.20736834139672225, 'CoM1Reval_1_MS_GLS': np.nan,
                         'tailLength_GLS': np.nan, 'mYErr_GLS': 0.202785969556528,
                         'CoM1Yeval_7_GLS_GLS': 0.2829103821808453, 'MoI1Yeval_22_GLS_GLS': 0.23212227571881977,
                         'CoM0YavgRes_GLS_GLS': 0.0050252905400519519, 'MoI0GavgRes_GLS_GLS': -0.05688431812117558,
                         'residsEarlyR_MS': np.nan, 'MoI0GstdTail_MS_GLS': np.nan, 'CoM0GstdTail_MS_GLS': np.nan,
                         'aYErr_MS': np.nan, 'MoI1RavgRes_GLS_MS': np.nan, 'tScale_GLS': 0.88857582042638272,
                         'CoM0GstdTail_MS_MS': 0.056881702923441711, 'sGaussHead_GLS': np.nan,
                         'MoI1Reval_10_GLS_GLS': 0.21795063918913996, 'mRmG_MS': -2.6392870963339612,
                         'aGaussHead_GLS': np.nan, 'CoM1GavgRes_MS_MS': -0.053169315455193746,
                         'residsEarlyG_GLS': 5.0664123053959331, 'CoM1Reval_7_GLS_MS': np.nan, 'rG_GLS': 1.0,
                         'MoI1Rtd_GLS_GLS': 15.941637140101728, 'devLength_GLS': np.nan,
                         'residsLateR_GLS': 48.415513706471387, 'aRErr_MS': 69044352.506834105,
                         'aSigHead_MS': 1.2820285259015392, 'MoI1GstdTail_MS_GLS': np.nan,
                         'CoM1Rtd_MS_MS': 7.8580611867023444, 'residsEarlyG_MS': np.nan,
                         'MoI0RavgRes_GLS_GLS': 0.0047263943402612099, 'MoI1Rscale_MS_MS': 1.0102760603353982,
                         'MoI1Reval_14_GLS_GLS': 0.23818475790484767, 'CoM1RavgRes_GLS_MS': np.nan,
                         'residsLateG_GLS': 2.3122898751147454, 'MoI0RstdTail_GLS_MS': np.nan,
                         'CoM0Gtd_MS_MS': 9.137077594465703, 'bYErr_MS': np.nan, 'CoM1Gtd_MS_GLS': np.nan,
                         'MoI1RavgRes_MS_GLS': np.nan, 'bG_GLS': np.nan, 'MoI1Geval_14_GLS_GLS': 0.21354670669762288,
                         'mYErr_MS': np.nan, 'MoI1Gscale_MS_MS': 1.0765541507431347, 'MoI1Geval_10_GLS_MS': np.nan,
                         'MoI1Yeval_14_GLS_GLS': 0.22195113943238579, 'MoI0Rtd_GLS_MS': np.nan,
                         'devLength_MS': 17.283217053063044, 'MoI1GstdTail_GLS_MS': np.nan,
                         'CoM0Rtd_GLS_GLS': 13.994681330902619, 'sR_MS': 2.7141462856148788,
                         'CoM1Yeval_14_GLS_MS': np.nan, 'residsLateY_GLS': 51.376579672990744,
                         'CoM1YavgRes_GLS_MS': np.nan, 'CoM0RstdTail_MS_MS': 0.12924853776448209,
                         'MoI1GavgRes_MS_MS': -0.0058619235428514424, 'MoI1Gtd_GLS_GLS': 15.485271248886118,
                         'rRErr_MS': 0.0, 'tailLength_MS': 186.91495593370414, 'MoI0GstdTail_GLS_MS': np.nan,
                         'CoM1Reval_15_MS_GLS': np.nan, 'mGErr_MS': 5.075682615954423, 'MoI1Geval_15_MS_GLS': np.nan,
                         'MoI1Reval_22_GLS_GLS': 0.24098998814683659}
        self.rnaiStd = {'bYErr_GLS': np.nan, 'MoI1Geval_7_GLS_GLS': 0.013315668475568435,
                        'residsEarlyY_GLS': 19.710681330847539, 'bY_GLS': np.nan,
                        'MoI0RstdTail_MS_MS': 0.022101249015179916, 'MoI1Gtd_GLS_MS': np.nan,
                        'MoI0Rscale_MS_GLS': np.nan, 'mRErr_GLS': 24.478729294172101,
                        'CoM0YstdTail_GLS_GLS': 0.023474606497852581, 'MoI1Reval_15_MS_GLS': np.nan,
                        'MoI1Rscale_GLS_GLS': 0.10003382048117804, 'MoI1Yeval_22_GLS_MS': np.nan,
                        'scaleHead_MS': 0.23805743414189515, 'bR_GLS': np.nan, 'rY_GLS': 0.0, 'tailHead_GLS': np.nan,
                        'MoI0RavgRes_MS_MS': 0.025039036638333602, 'tMove_GLS': 3.5001300738140317,
                        'MoI0Gtd_MS_MS': 1.2078131723153733, 'MoI1Gtd_MS_GLS': np.nan, 'MoI0Gtd_MS_GLS': np.nan,
                        'sYErr_GLS': 0.10328221884890527, 'MoI1Yeval_10_GLS_MS': np.nan, 'MoI1Reval_3_MS_GLS': np.nan,
                        'MoI1Geval_1_MS_GLS': np.nan, 'CoM0RavgRes_GLS_MS': np.nan, 'MoI1Rtd_MS_GLS': np.nan,
                        'CoM0Ytd_GLS_MS': np.nan, 'MoI1RstdTail_GLS_GLS': 0.011978732498279337,
                        'CoM1Reval_22_GLS_GLS': 0.057226745482188091, 'MoI1Yeval_7_GLS_GLS': 0.020608674825176861,
                        'MoI0Gscale_GLS_GLS': 0.10840188648693681, 'rRErr_GLS': 0.0, 'MoI1Reval_7_GLS_MS': np.nan,
                        'CoM1Geval_15_MS_MS': 0.039772809378581304, 'CoM1RavgRes_MS_MS': 0.071192672774666607,
                        'sRErr_MS': 0.42776590959482952, 'MoI0Gscale_MS_GLS': np.nan, 'MoI1Gscale_GLS_MS': np.nan,
                        'MoI1Geval_3_MS_GLS': np.nan, 'GtailErr_MS': np.nan, 'CoM1Rtd_MS_GLS': np.nan,
                        'CoM0GavgRes_MS_GLS': np.nan, 'CoM1GstdTail_MS_GLS': np.nan, 'MoI0Yscale_GLS_MS': np.nan,
                        'aY_GLS': 16.923422534384763, 'aYErr_GLS': 2.809441188463949, 'MoI1RstdTail_GLS_MS': np.nan,
                        'rG_MS': np.nan, 'MoI0RavgRes_MS_GLS': np.nan, 'CoM1Geval_1_MS_MS': 0.0, 'bY_MS': np.nan,
                        'MoI1GstdTail_GLS_GLS': 0.027403630115529754, 'MoI0GavgRes_MS_MS': 0.036824902659726072,
                        'MoI1GavgRes_MS_GLS': np.nan, 'mSigHead_MS': 3.9233046280411821,
                        'CoM1YavgRes_GLS_GLS': 0.055910321832127928, 'CoM1Geval_15_MS_GLS': np.nan,
                        'CoM1Geval_3_MS_MS': 0.038925529357094625, 'CoM1Reval_3_MS_GLS': np.nan, 'rY_MS': np.nan,
                        'CoM1Yeval_10_GLS_GLS': 0.046313366515791904, 'CoM1Reval_7_MS_GLS': np.nan,
                        'GtailErr_GLS': 4.8341823607870271, 'CoM0Rtd_GLS_MS': np.nan, 'aG_GLS': 17.486297317640243,
                        'MoI1Geval_10_GLS_GLS': 0.023550990937677489, 'mGaussHead_MS': 3.6117650393366922,
                        'MoI0Gtd_GLS_GLS': 2.2005542757572685, 'MoI1RavgRes_GLS_GLS': 0.028029432344795154,
                        'MoI1Yscale_GLS_MS': np.nan, 'mRmG_GLS': 1.8527397615635328, 'CoM1Reval_22_GLS_MS': np.nan,
                        'CoM0Rtd_MS_MS': 2.0126434735949674, 'MoI1Geval_7_GLS_MS': np.nan, 'CoM0RavgRes_MS_GLS': np.nan,
                        'MoI1Rtd_GLS_MS': np.nan, 'fracG_MS': np.nan, 'Rtail_MS': np.nan, 'aG_MS': 130141233.79015468,
                        'MoI0GavgRes_MS_GLS': np.nan, 'residsEarlyY_MS': np.nan,
                        'CoM1Rscale_MS_MS': 0.23996892736737191, 'maxG_MS': 174041536.16820961,
                        'MoI0Rscale_GLS_MS': np.nan, 'sG_MS': 0.81956991772772325, 'MoI1Yeval_14_GLS_MS': np.nan,
                        'MoI0Rtd_GLS_GLS': 2.3488226791892357, 'movement_MS': 0.44746286140096309,
                        'CoM1Rtd_GLS_GLS': 3.1628731438361628, 'CoM0Rscale_GLS_MS': np.nan, 'MoI1Ytd_GLS_MS': np.nan,
                        'mRErr_MS': 4.415408013934381, 'CoM1Gscale_MS_MS': 0.3632420970987551, 'MoI0Ytd_GLS_MS': np.nan,
                        'maxR_GLS': 308498175.09665745, 'tailHead_MS': 0.36741683481391191, 'RtailErr_MS': np.nan,
                        'bR_MS': 151449487.43374357, 'MoI0YstdTail_GLS_GLS': 0.0196579656651154, 'fracY_MS': np.nan,
                        'CoM0Rscale_MS_GLS': np.nan, 'bGErr_GLS': np.nan, 'sSigHead_GLS': np.nan,
                        'CoM1Yeval_10_GLS_MS': np.nan, 'CoM1Geval_1_MS_GLS': np.nan, 'CoM0Gtd_MS_GLS': np.nan,
                        'CoM0Rscale_MS_MS': 0.2452595953329596, 'CoM1Yeval_22_GLS_MS': np.nan,
                        'MoI0Rscale_GLS_GLS': 0.11404103259796251, 'movement_GLS': 0.48655010019524192,
                        'CoM0Rtd_MS_GLS': np.nan, 'YtailErr_MS': np.nan, 'CoM0Yscale_GLS_MS': np.nan,
                        'CoM1Rscale_GLS_MS': np.nan, 'MoI1Reval_3_MS_MS': 0.013234513732946353,
                        'MoI1Reval_7_MS_GLS': np.nan, 'rGErr_GLS': 0.0, 'MoI0RavgRes_GLS_MS': np.nan,
                        'bG_MS': 169716688.87625971, 'CoM1Reval_10_GLS_GLS': 0.036646133592636132,
                        'devHead_MS': 4.9862250275165998, 'CoM1Geval_3_MS_GLS': np.nan, 'residsLateR_MS': np.nan,
                        'CoM0RstdTail_GLS_MS': np.nan, 'mGErr_GLS': 214833.01881636577, 'CoM1RstdTail_GLS_MS': np.nan,
                        'maxG_GLS': 254717970.42858875, 'aR_MS': 148582659.77618688,
                        'MoI0YavgRes_GLS_GLS': 0.050152919977780298, 'CoM0RstdTail_GLS_GLS': 0.019898515814076468,
                        'tMove_MS': 1.1244055236233221, 'rYErr_MS': np.nan, 'CoM1GavgRes_MS_GLS': np.nan,
                        'MoI1Gscale_MS_GLS': np.nan, 'MoI1Reval_1_MS_MS': 0.0,
                        'CoM1YstdTail_GLS_GLS': 0.022549850753315104, 'MoI1Geval_14_GLS_MS': np.nan,
                        'residsLateY_MS': np.nan, 'CoM1Reval_7_MS_MS': 0.054236999828280563, 'mYmG_MS': np.nan,
                        'mY_GLS': 1.2493955812642508, 'MoI0GavgRes_GLS_MS': np.nan, 'MoI1Rscale_GLS_MS': np.nan,
                        'MoI1Reval_10_GLS_MS': np.nan, 'sGaussHead_MS': 1.9106406448896742,
                        'maxR_MS': 164276434.91912168, 'Gtail_MS': np.nan, 'sY_MS': np.nan, 'CoM1Gscale_MS_GLS': np.nan,
                        'MoI0Rtd_MS_MS': 1.5221191842980717, 'mYmG_GLS': 1.3176758876091448,
                        'MoI1YavgRes_GLS_MS': np.nan, 'MoI1GavgRes_GLS_MS': np.nan, 'CoM1Rscale_MS_GLS': np.nan,
                        'MoI1RavgRes_MS_MS': 0.026807271434306713, 'sGErr_GLS': 118311.42700690309,
                        'Rtail_GLS': 37.627623511087762, 'sG_GLS': 0.5484943947645371,
                        'CoM1RstdTail_MS_MS': 0.069275926720849673, 'MoI1YstdTail_GLS_GLS': 0.010615438634428423,
                        'MoI0GstdTail_MS_MS': 0.016331871410916906, 'tScale_MS': 0.11177310553694865,
                        'MoI1Geval_1_MS_MS': 0.0, 'mG_GLS': 1.5385198690033333, 'bGErr_MS': 21401777.298232034,
                        'sGErr_MS': 108.68137219420564, 'sR_GLS': 0.52374810360995494, 'CoM1Geval_7_MS_GLS': np.nan,
                        'MoI1Gscale_GLS_GLS': 0.11382362625532613, 'MoI1Rtd_MS_MS': 2.1805641445588857,
                        'mR_GLS': 0.69162376999836783, 'CoM0Rscale_GLS_GLS': 0.15029711172122709,
                        'MoI1Ytd_GLS_GLS': 2.3771893984538597, 'CoM0Gscale_MS_MS': 0.35938622674824189,
                        'CoM1RstdTail_MS_GLS': np.nan, 'residsLateG_MS': np.nan, 'CoM1Yscale_GLS_MS': np.nan,
                        'MoI0Gscale_GLS_MS': np.nan, 'CoM1Geval_7_MS_MS': 0.023998067267413322,
                        'CoM0Gscale_MS_GLS': np.nan, 'MoI1Geval_3_MS_MS': 0.030758806945028139,
                        'CoM1Reval_14_GLS_MS': np.nan, 'scaleLength_GLS': np.nan, 'Gtail_GLS': 19.258251285340187,
                        'MoI1YstdTail_GLS_MS': np.nan, 'sYErr_MS': np.nan, 'Ytail_MS': np.nan,
                        'MoI1Reval_14_GLS_MS': np.nan, 'aR_GLS': 19.68948542316668, 'MoI1Geval_7_MS_GLS': np.nan,
                        'CoM1YstdTail_GLS_MS': np.nan, 'MoI1GstdTail_MS_MS': 0.014974525495670892,
                        'fracG_GLS': 0.053267912585910009, 'mY_MS': np.nan, 'CoM0YstdTail_GLS_MS': np.nan,
                        'sY_GLS': 0.50185041417563359, 'CoM1Reval_14_GLS_GLS': 0.031495185342358303,
                        'MoI0YavgRes_GLS_MS': np.nan, 'fracY_GLS': 0.067664088775902603, 'CoM1Rtd_GLS_MS': np.nan,
                        'RtailErr_GLS': 7.7160298665543259, 'CoM1Reval_7_GLS_GLS': 0.032803536317039958,
                        'aY_MS': np.nan, 'mR_MS': 1.9947248119277994, 'MoI1Reval_7_GLS_GLS': 0.022023601052849985,
                        'sRErr_GLS': 1.6641338441699942, 'mGaussHead_GLS': np.nan, 'devHead_GLS': np.nan,
                        'CoM1Gtd_MS_MS': 1.0413083391277533, 'CoM1Reval_10_GLS_MS': np.nan, 'mSigHead_GLS': np.nan,
                        'CoM0RstdTail_MS_GLS': np.nan, 'aGErr_MS': 103737568.99064626, 'fracR_MS': np.nan,
                        'MoI1YavgRes_GLS_GLS': 0.034244320165720578, 'mYmR_GLS': 0.85505694079174377,
                        'maxHead_MS': 0.4076940571596101, 'rR_MS': np.nan, 'CoM0RavgRes_GLS_GLS': 0.052835269949788134,
                        'CoM1Yeval_7_GLS_MS': np.nan, 'scaleHead_GLS': np.nan, 'fracR_GLS': 0.072259291077765769,
                        'MoI1Reval_7_MS_MS': 0.019999309326649423, 'Ytail_GLS': 39.71859145312694,
                        'MoI1Geval_22_GLS_MS': np.nan, 'MoI0RstdTail_GLS_GLS': 0.022515787703298592,
                        'maxHead_GLS': np.nan, 'bRErr_MS': 484181352.01865989, 'rGErr_MS': 0.0,
                        'CoM1RavgRes_GLS_GLS': 0.034918553949632801, 'MoI1Rscale_MS_GLS': np.nan,
                        'CoM1Yscale_GLS_GLS': 0.15847439322619858, 'bRErr_GLS': np.nan,
                        'CoM0RavgRes_MS_MS': 0.097011225448848354, 'MoI1Yscale_GLS_GLS': 0.11645893834111688,
                        'CoM0YavgRes_GLS_MS': np.nan, 'MoI1Yeval_7_GLS_MS': np.nan,
                        'MoI0Rscale_MS_MS': 0.15509582818600129, 'aGaussHead_MS': 0.36537163358313979,
                        'mYmR_MS': np.nan, 'CoM1Reval_3_MS_MS': 0.036741712669728449, 'CoM1Reval_1_MS_MS': 0.0,
                        'CoM1Ytd_GLS_GLS': 3.278607810959659, 'CoM1Reval_15_MS_MS': 0.08714580653537228,
                        'MoI1RstdTail_MS_GLS': np.nan, 'MoI1GavgRes_GLS_GLS': 0.043444246644953451,
                        'MoI1Geval_7_MS_MS': 0.027778082017463814, 'sSigHead_MS': 4.3626864393813154,
                        'MoI0Ytd_GLS_GLS': 1.9869588694640588, 'MoI0YstdTail_GLS_MS': np.nan,
                        'MoI0RstdTail_MS_GLS': np.nan, 'MoI1Reval_15_MS_MS': 0.028425329997792251,
                        'CoM1RstdTail_GLS_GLS': 0.01487956537325784, 'CoM1RavgRes_MS_GLS': np.nan, 'rR_GLS': 0.0,
                        'MoI0Gtd_GLS_MS': np.nan, 'scaleLength_MS': 0.0695241625083549,
                        'residsEarlyR_GLS': 6.7699949400837163, 'MoI0Yscale_GLS_GLS': 0.10016935388971838,
                        'CoM0GavgRes_MS_MS': 0.059342705493761243, 'MoI1Reval_1_MS_GLS': np.nan,
                        'MoI1Geval_15_MS_MS': 0.039151174550097416, 'MoI1RstdTail_MS_MS': 0.01676494980287031,
                        'MoI1Reval_22_GLS_MS': np.nan, 'CoM1Ytd_GLS_MS': np.nan,
                        'MoI0Gscale_MS_MS': 0.18911186009322264, 'YtailErr_GLS': 6.698893514920174,
                        'MoI0Rtd_MS_GLS': np.nan, 'CoM1Yeval_22_GLS_GLS': 0.068549381255150671,
                        'MoI1Yeval_10_GLS_GLS': 0.026146541348139363, 'CoM1Rscale_GLS_GLS': 0.19068532622958839,
                        'mG_MS': 1.3176724873439174, 'CoM0Ytd_GLS_GLS': 3.1501612019727543,
                        'CoM0Yscale_GLS_GLS': 0.17940058759017369, 'MoI0GstdTail_GLS_GLS': 0.024675344095370023,
                        'rYErr_GLS': 0.0, 'CoM1GstdTail_MS_MS': 0.032829344303273042, 'aGErr_GLS': 57.803084713457089,
                        'MoI1Geval_22_GLS_GLS': 0.062365002492625034, 'aRErr_GLS': 78.526420014853031,
                        'MoI1Gtd_MS_MS': 1.2169134902382301, 'aSigHead_GLS': np.nan,
                        'CoM1Yeval_14_GLS_GLS': 0.060132511756171667, 'CoM1Reval_1_MS_GLS': np.nan,
                        'tailLength_GLS': np.nan, 'mYErr_GLS': 0.14817061251560831,
                        'CoM1Yeval_7_GLS_GLS': 0.04590304953575431, 'MoI1Yeval_22_GLS_GLS': 0.047982381256088162,
                        'CoM0YavgRes_GLS_GLS': 0.0547892678021695, 'MoI0GavgRes_GLS_GLS': 0.074249997052172181,
                        'residsEarlyR_MS': np.nan, 'MoI0GstdTail_MS_GLS': np.nan, 'CoM0GstdTail_MS_GLS': np.nan,
                        'aYErr_MS': np.nan, 'MoI1RavgRes_GLS_MS': np.nan, 'tScale_GLS': 0.1894896746182769,
                        'CoM0GstdTail_MS_MS': 0.049642809663482092, 'sGaussHead_GLS': np.nan,
                        'MoI1Reval_10_GLS_GLS': 0.43558916529820435, 'mRmG_MS': 1.8386562156477773,
                        'aGaussHead_GLS': np.nan, 'CoM1GavgRes_MS_MS': 0.038634034795517198,
                        'residsEarlyG_GLS': 1.6206358713816791, 'CoM1Reval_7_GLS_MS': np.nan, 'rG_GLS': 0.0,
                        'MoI1Rtd_GLS_GLS': 2.4246447835490361, 'devLength_GLS': np.nan,
                        'residsLateR_GLS': 35.930996797890181, 'aRErr_MS': 1028627385.3988408,
                        'aSigHead_MS': 0.31018922142655631, 'MoI1GstdTail_MS_GLS': np.nan,
                        'CoM1Rtd_MS_MS': 1.6864354713295786, 'residsEarlyG_MS': np.nan,
                        'MoI0RavgRes_GLS_GLS': 0.065036734448042435, 'MoI1Rscale_MS_MS': 0.1865851099698104,
                        'MoI1Reval_14_GLS_GLS': 0.25259644191420427, 'CoM1RavgRes_GLS_MS': np.nan,
                        'residsLateG_GLS': 25.523854217455344, 'MoI0RstdTail_GLS_MS': np.nan,
                        'CoM0Gtd_MS_MS': 1.6006575652002757, 'bYErr_MS': np.nan, 'CoM1Gtd_MS_GLS': np.nan,
                        'MoI1RavgRes_MS_GLS': np.nan, 'bG_GLS': np.nan, 'MoI1Geval_14_GLS_GLS': 0.048749536511729259,
                        'mYErr_MS': np.nan, 'MoI1Gscale_MS_MS': 0.23321863739934817, 'MoI1Geval_10_GLS_MS': np.nan,
                        'MoI1Yeval_14_GLS_GLS': 0.036017736924950354, 'MoI0Rtd_GLS_MS': np.nan,
                        'devLength_MS': 4.233778371000759, 'MoI1GstdTail_GLS_MS': np.nan,
                        'CoM0Rtd_GLS_GLS': 3.8002281803294502, 'sR_MS': 0.52386508654930297,
                        'CoM1Yeval_14_GLS_MS': np.nan, 'residsLateY_GLS': 171.84001610299094,
                        'CoM1YavgRes_GLS_MS': np.nan, 'CoM0RstdTail_MS_MS': 0.066331074061997009,
                        'MoI1GavgRes_MS_MS': 0.032807946911581579, 'MoI1Gtd_GLS_GLS': 2.4700695090340368,
                        'rRErr_MS': 0.0, 'tailLength_MS': 16.380751565953645, 'MoI0GstdTail_GLS_MS': np.nan,
                        'CoM1Reval_15_MS_GLS': np.nan, 'mGErr_MS': 104.73578695267206, 'MoI1Geval_15_MS_GLS': np.nan,
                        'MoI1Reval_22_GLS_GLS': 0.04444563641805712}
        #
        # KStests = {'bYErr_GLS': 0, 'MoI1Geval_7_GLS_GLS': 0.13646445975466304, 'residsEarlyY_GLS': 0.10458251246349637,
        #            'bY_GLS': 0, 'MoI0RstdTail_MS_MS': 0.24025469168900804, 'MoI1Gtd_GLS_MS': 0, 'MoI0Rscale_MS_GLS': 0,
        #            'mRErr_GLS': 0.15470632662101136, 'CoM0YstdTail_GLS_GLS': 0.18888992392485637,
        #            'MoI1Reval_15_MS_GLS': 0, 'MoI1Rscale_GLS_GLS': 0.14645691777143868, 'MoI1Yeval_22_GLS_MS': 0,
        #            'scaleHead_MS': 0.10157858613589565, 'bR_GLS': 0, 'rY_GLS': 0.0, 'tailHead_GLS': 0,
        #            'MoI0RavgRes_MS_MS': 0.069477211796246607, 'tMove_GLS': 0.19637633907778296,
        #            'MoI0Gtd_MS_MS': 0.09711796246648792, 'MoI1Gtd_MS_GLS': 0, 'MoI0Gtd_MS_GLS': 0,
        #            'sYErr_GLS': 0.15428795074812773, 'MoI1Yeval_10_GLS_MS': 0, 'MoI1Reval_3_MS_GLS': 0,
        #            'MoI1Geval_1_MS_GLS': 0, 'CoM0RavgRes_GLS_MS': 0, 'MoI1Rtd_MS_GLS': 0, 'CoM0Ytd_GLS_MS': 0,
        #            'MoI1RstdTail_GLS_GLS': 0.30117373078714488, 'CoM1Reval_22_GLS_GLS': 0.089120140077201704,
        #            'MoI1Yeval_7_GLS_GLS': 0.1980418545185384, 'MoI0Gscale_GLS_GLS': 0.1265167011730087,
        #            'rRErr_GLS': 0.0, 'MoI1Reval_7_GLS_MS': 0, 'CoM1Geval_15_MS_MS': 0.11700901633787542,
        #            'CoM1RavgRes_MS_MS': 0.076823056300268089, 'sRErr_MS': 0.10270963498849062, 'MoI0Gscale_MS_GLS': 0,
        #            'MoI1Gscale_GLS_MS': 0, 'MoI1Geval_3_MS_GLS': 0, 'GtailErr_MS': 0, 'CoM1Rtd_MS_GLS': 0,
        #            'CoM0GavgRes_MS_GLS': 0, 'CoM1GstdTail_MS_GLS': 0, 'MoI0Yscale_GLS_MS': 0,
        #            'aY_GLS': 0.093155863884634482, 'aYErr_GLS': 0.17397813857990851, 'MoI1RstdTail_GLS_MS': 0,
        #            'rG_MS': 0, 'MoI0RavgRes_MS_GLS': 0, 'CoM1Geval_1_MS_MS': 0, 'bY_MS': 0,
        #            'MoI1GstdTail_GLS_GLS': 0.3021114733736997, 'MoI0GavgRes_MS_MS': 0.10276139410187668,
        #            'MoI1GavgRes_MS_GLS': 0, 'mSigHead_MS': 0.18992887795216518,
        #            'CoM1YavgRes_GLS_GLS': 0.18128861977953736, 'CoM1Geval_15_MS_GLS': 0,
        #            'CoM1Geval_3_MS_MS': 0.19878846397002198, 'CoM1Reval_3_MS_GLS': 0, 'rY_MS': 0,
        #            'CoM1Yeval_10_GLS_GLS': 0.14584045247394795, 'CoM1Reval_7_MS_GLS': 0,
        #            'GtailErr_GLS': 0.28653306662868394, 'CoM0Rtd_GLS_MS': 0, 'aG_GLS': 0.26074944745830819,
        #            'MoI1Geval_10_GLS_GLS': 0.10094267025779147, 'mGaussHead_MS': 0.19689036966523876,
        #            'MoI0Gtd_GLS_GLS': 0.2840301338142493, 'MoI1RavgRes_GLS_GLS': 0.20542462350566681,
        #            'MoI1Yscale_GLS_MS': 0, 'mRmG_GLS': 0.23901544000553904, 'CoM1Reval_22_GLS_MS': 0,
        #            'CoM0Rtd_MS_MS': 0.11138069705093834, 'MoI1Geval_7_GLS_MS': 0, 'CoM0RavgRes_MS_GLS': 0,
        #            'MoI1Rtd_GLS_MS': 0, 'fracG_MS': 0, 'Rtail_MS': 0, 'aG_MS': 0.25491053677932407,
        #            'MoI0GavgRes_MS_GLS': 0, 'residsEarlyY_MS': 0, 'CoM1Rscale_MS_MS': 0.087426273458445003,
        #            'maxG_MS': 0.1816929133858268, 'MoI0Rscale_GLS_MS': 0, 'sG_MS': 0.1644207723035952,
        #            'MoI1Yeval_14_GLS_MS': 0, 'MoI0Rtd_GLS_GLS': 0.18449666783584706, 'movement_MS': 0.23064304461942256,
        #            'CoM1Rtd_GLS_GLS': 0.15209381586339835, 'CoM0Rscale_GLS_MS': 0, 'MoI1Ytd_GLS_MS': 0,
        #            'mRErr_MS': 0.075938543754175025, 'CoM1Gscale_MS_MS': 0.074785522788203718, 'MoI0Ytd_GLS_MS': 0,
        #            'maxR_GLS': 0.20327353121497183, 'tailHead_MS': 0.088064516129032322, 'RtailErr_MS': 0,
        #            'bR_MS': 0.12456047587574355, 'MoI0YstdTail_GLS_GLS': 0.23913077264998711, 'fracY_MS': 0,
        #            'CoM0Rscale_MS_GLS': 0, 'bGErr_GLS': 0, 'sSigHead_GLS': 0, 'CoM1Yeval_10_GLS_MS': 0,
        #            'CoM1Geval_1_MS_GLS': 0, 'CoM0Gtd_MS_GLS': 0, 'CoM0Rscale_MS_MS': 0.12010723860589811,
        #            'CoM1Yeval_22_GLS_MS': 0, 'MoI0Rscale_GLS_GLS': 0.096134565421851614,
        #            'movement_GLS': 0.31638021746755524, 'CoM0Rtd_MS_GLS': 0, 'YtailErr_MS': 0, 'CoM0Yscale_GLS_MS': 0,
        #            'CoM1Rscale_GLS_MS': 0, 'MoI1Reval_3_MS_MS': 0.13007404901229536, 'MoI1Reval_7_MS_GLS': 0,
        #            'rGErr_GLS': 0.0, 'MoI0RavgRes_GLS_MS': 0, 'bG_MS': 0.14943159286186383,
        #            'CoM1Reval_10_GLS_GLS': 0.090075510680555526, 'devHead_MS': 0.14427613941018766,
        #            'CoM1Geval_3_MS_GLS': 0, 'residsLateR_MS': 0, 'CoM0RstdTail_GLS_MS': 0,
        #            'mGErr_GLS': 0.25534994453913368, 'CoM1RstdTail_GLS_MS': 0, 'maxG_GLS': 0.24553651603633955,
        #            'aR_MS': 0.13876902713434813, 'MoI0YavgRes_GLS_GLS': 0.11326314611561306,
        #            'CoM0RstdTail_GLS_GLS': 0.18568850427257508, 'tMove_MS': 0.054926931106471732, 'rYErr_MS': 0,
        #            'CoM1GavgRes_MS_GLS': 0, 'MoI1Gscale_MS_GLS': 0, 'MoI1Reval_1_MS_MS': 0,
        #            'CoM1YstdTail_GLS_GLS': 0.14775966464834653, 'MoI1Geval_14_GLS_MS': 0, 'residsLateY_MS': 0,
        #            'CoM1Reval_7_MS_MS': 0.091376146788990809, 'mYmG_MS': 0, 'mY_GLS': 0.11338875093081749,
        #            'MoI0GavgRes_GLS_MS': 0, 'MoI1Rscale_GLS_MS': 0, 'MoI1Reval_10_GLS_MS': 0,
        #            'sGaussHead_MS': 0.13317757009345793, 'maxR_MS': 0.13110236220472438, 'Gtail_MS': 0, 'sY_MS': 0,
        #            'CoM1Gscale_MS_GLS': 0, 'MoI0Rtd_MS_MS': 0.064088471849865936, 'mYmG_GLS': 0.20599008834302951,
        #            'MoI1YavgRes_GLS_MS': 0, 'MoI1GavgRes_GLS_MS': 0, 'CoM1Rscale_MS_GLS': 0,
        #            'MoI1RavgRes_MS_MS': 0.06075067024128683, 'sGErr_GLS': 0.20850321286730267,
        #            'Rtail_GLS': 0.29610326457963826, 'sG_GLS': 0.20148649616734726,
        #            'CoM1RstdTail_MS_MS': 0.23317694369973191, 'MoI1YstdTail_GLS_GLS': 0.25856368698419013,
        #            'MoI0GstdTail_MS_MS': 0.1026005361930295, 'tScale_MS': 0.2277479892761394, 'MoI1Geval_1_MS_MS': 0,
        #            'mG_GLS': 0.17372185614085256, 'bGErr_MS': 0.18306780776826859, 'sGErr_MS': 0.18129733289430361,
        #            'sR_GLS': 0.24513617403034527, 'CoM1Geval_7_MS_GLS': 0, 'MoI1Gscale_GLS_GLS': 0.1731390631450539,
        #            'MoI1Rtd_MS_MS': 0.17454423592493301, 'mR_GLS': 0.14975055259415448,
        #            'CoM0Rscale_GLS_GLS': 0.14567160102184296, 'MoI1Ytd_GLS_GLS': 0.22685721203241838,
        #            'CoM0Gscale_MS_MS': 0.087788203753351196, 'CoM1RstdTail_MS_GLS': 0, 'residsLateG_MS': 0,
        #            'CoM1Yscale_GLS_MS': 0, 'MoI0Gscale_GLS_MS': 0, 'CoM1Geval_7_MS_MS': 0.091228367831693302,
        #            'CoM0Gscale_MS_GLS': 0, 'MoI1Geval_3_MS_MS': 0.08154304171408211, 'CoM1Reval_14_GLS_MS': 0,
        #            'scaleLength_GLS': 0, 'Gtail_GLS': 0.21408006088944759, 'MoI1YstdTail_GLS_MS': 0, 'sYErr_MS': 0,
        #            'Ytail_MS': 0, 'MoI1Reval_14_GLS_MS': 0, 'aR_GLS': 0.074505949269214999, 'MoI1Geval_7_MS_GLS': 0,
        #            'CoM1YstdTail_GLS_MS': 0, 'MoI1GstdTail_MS_MS': 0.13804289544235926,
        #            'fracG_GLS': 0.25104813688636773, 'mY_MS': 0, 'CoM0YstdTail_GLS_MS': 0, 'sY_GLS': 0.2205208911540596,
        #            'CoM1Reval_14_GLS_GLS': 0.14325637101603561, 'MoI0YavgRes_GLS_MS': 0,
        #            'fracY_GLS': 0.10986703988803359, 'CoM1Rtd_GLS_MS': 0, 'RtailErr_GLS': 0.26106504331886277,
        #            'CoM1Reval_7_GLS_GLS': 0.10114055686011406, 'aY_MS': 0, 'mR_MS': 0.1356921766757832,
        #            'MoI1Reval_7_GLS_GLS': 0.12890425674768358, 'sRErr_GLS': 0.21930727716689735, 'mGaussHead_GLS': 0,
        #            'devHead_GLS': 0, 'CoM1Gtd_MS_MS': 0.091085790884718507, 'CoM1Reval_10_GLS_MS': 0, 'mSigHead_GLS': 0,
        #            'CoM0RstdTail_MS_GLS': 0, 'aGErr_MS': 0.21627263747118863, 'fracR_MS': 0,
        #            'MoI1YavgRes_GLS_GLS': 0.22646850452421197, 'mYmR_GLS': 0.087373756132809732,
        #            'maxHead_MS': 0.13647453083109917, 'rR_MS': 0, 'CoM0RavgRes_GLS_GLS': 0.18202940326834127,
        #            'CoM1Yeval_7_GLS_MS': 0, 'scaleHead_GLS': 0, 'fracR_GLS': 0.11145550140267657,
        #            'MoI1Reval_7_MS_MS': 0.19143730886850152, 'Ytail_GLS': 0.33175702864486667, 'MoI1Geval_22_GLS_MS': 0,
        #            'MoI0RstdTail_GLS_GLS': 0.21333643844123584, 'maxHead_GLS': 0, 'bRErr_MS': 0.23508218277449044,
        #            'rGErr_MS': 0, 'CoM1RavgRes_GLS_GLS': 0.1396149666200901, 'MoI1Rscale_MS_GLS': 0,
        #            'CoM1Yscale_GLS_GLS': 0.17371020259930531, 'bRErr_GLS': 0, 'CoM0RavgRes_MS_MS': 0.14017426273458444,
        #            'MoI1Yscale_GLS_GLS': 0.11340222185663773, 'CoM0YavgRes_GLS_MS': 0, 'MoI1Yeval_7_GLS_MS': 0,
        #            'MoI0Rscale_MS_MS': 0.15982573726541555, 'aGaussHead_MS': 0.076959952798879017, 'mYmR_MS': 0,
        #            'CoM1Reval_3_MS_MS': 0.10278433513521829, 'CoM1Reval_1_MS_MS': 0,
        #            'CoM1Ytd_GLS_GLS': 0.22152140376144996, 'CoM1Reval_15_MS_MS': 0.085321298087255526,
        #            'MoI1RstdTail_MS_GLS': 0, 'MoI1GavgRes_GLS_GLS': 0.21256016146561094,
        #            'MoI1Geval_7_MS_MS': 0.10565999321343739, 'sSigHead_MS': 0.27130124777183595,
        #            'MoI0Ytd_GLS_GLS': 0.1701157488600491, 'MoI0YstdTail_GLS_MS': 0, 'MoI0RstdTail_MS_GLS': 0,
        #            'MoI1Reval_15_MS_MS': 0.071442306604680594, 'CoM1RstdTail_GLS_GLS': 0.10137556280080734,
        #            'CoM1RavgRes_MS_GLS': 0, 'rR_GLS': 0.0, 'MoI0Gtd_GLS_MS': 0, 'scaleLength_MS': 0.092761394101876671,
        #            'residsEarlyR_GLS': 0.11045605840898343, 'MoI0Yscale_GLS_GLS': 0.066565682589248304,
        #            'CoM0GavgRes_MS_MS': 0.093431635388739942, 'MoI1Reval_1_MS_GLS': 0,
        #            'MoI1Geval_15_MS_MS': 0.11467019185139987, 'MoI1RstdTail_MS_MS': 0.17012064343163541,
        #            'MoI1Reval_22_GLS_MS': 0, 'CoM1Ytd_GLS_MS': 0, 'MoI0Gscale_MS_MS': 0.14529490616621982,
        #            'YtailErr_GLS': 0.30614267047319432, 'MoI0Rtd_MS_GLS': 0,
        #            'CoM1Yeval_22_GLS_GLS': 0.066943813055788171, 'MoI1Yeval_10_GLS_GLS': 0.1151468513530845,
        #            'CoM1Rscale_GLS_GLS': 0.11162672311842142, 'mG_MS': 0.07942110836569008,
        #            'CoM0Ytd_GLS_GLS': 0.17214577713766013, 'CoM0Yscale_GLS_GLS': 0.099080279237528823,
        #            'MoI0GstdTail_GLS_GLS': 0.27660922217047046, 'rYErr_GLS': 0.0,
        #            'CoM1GstdTail_MS_MS': 0.1378686327077748, 'aGErr_GLS': 0.098094710314161726,
        #            'MoI1Geval_22_GLS_GLS': 0.12922917331631026, 'aRErr_GLS': 0.19950230959453785,
        #            'MoI1Gtd_MS_MS': 0.11175603217158178, 'aSigHead_GLS': 0,
        #            'CoM1Yeval_14_GLS_GLS': 0.087585031187794704, 'CoM1Reval_1_MS_GLS': 0, 'tailLength_GLS': 0,
        #            'mYErr_GLS': 0.16750925469884492, 'CoM1Yeval_7_GLS_GLS': 0.13687203157401834,
        #            'MoI1Yeval_22_GLS_GLS': 0.10236909476379052, 'CoM0YavgRes_GLS_GLS': 0.17849091755938518,
        #            'MoI0GavgRes_GLS_GLS': 0.29033379909951873, 'residsEarlyR_MS': 0, 'MoI0GstdTail_MS_GLS': 0,
        #            'CoM0GstdTail_MS_GLS': 0, 'aYErr_MS': 0, 'MoI1RavgRes_GLS_MS': 0, 'tScale_GLS': 0.28117133252421911,
        #            'CoM0GstdTail_MS_MS': 0.096233243967828425, 'sGaussHead_GLS': 0,
        #            'MoI1Reval_10_GLS_GLS': 0.061729185610184034, 'mRmG_MS': 0.12839160839160846, 'aGaussHead_GLS': 0,
        #            'CoM1GavgRes_MS_MS': 0.093780160857908856, 'residsEarlyG_GLS': 0.1089710842532452,
        #            'CoM1Reval_7_GLS_MS': 0, 'rG_GLS': 0.0, 'MoI1Rtd_GLS_GLS': 0.21458081778478594, 'devLength_GLS': 0,
        #            'residsLateR_GLS': 0.15239325048473146, 'aRErr_MS': 0.19556213017751478,
        #            'aSigHead_MS': 0.1124805424654029, 'MoI1GstdTail_MS_GLS': 0, 'CoM1Rtd_MS_MS': 0.10829758713136728,
        #            'residsEarlyG_MS': 0, 'MoI0RavgRes_GLS_GLS': 0.17198571650364847,
        #            'MoI1Rscale_MS_MS': 0.10143431635388744, 'MoI1Reval_14_GLS_GLS': 0.15855790920861196,
        #            'CoM1RavgRes_GLS_MS': 0, 'residsLateG_GLS': 0.21522963690669539, 'MoI0RstdTail_GLS_MS': 0,
        #            'CoM0Gtd_MS_MS': 0.20697050938337802, 'bYErr_MS': 0, 'CoM1Gtd_MS_GLS': 0, 'MoI1RavgRes_MS_GLS': 0,
        #            'bG_GLS': 0, 'MoI1Geval_14_GLS_GLS': 0.18993976688481545, 'mYErr_MS': 0,
        #            'MoI1Gscale_MS_MS': 0.12117962466487936, 'MoI1Geval_10_GLS_MS': 0,
        #            'MoI1Yeval_14_GLS_GLS': 0.23141576204151637, 'MoI0Rtd_GLS_MS': 0,
        #            'devLength_MS': 0.14345844504021449, 'MoI1GstdTail_GLS_MS': 0,
        #            'CoM0Rtd_GLS_GLS': 0.19915074046367834, 'sR_MS': 0.082999999999999963, 'CoM1Yeval_14_GLS_MS': 0,
        #            'residsLateY_GLS': 0.19931797924492264, 'CoM1YavgRes_GLS_MS': 0,
        #            'CoM0RstdTail_MS_MS': 0.18664879356568365, 'MoI1GavgRes_MS_MS': 0.13412868632707775,
        #            'MoI1Gtd_GLS_GLS': 0.28168349562178152, 'rRErr_MS': 0, 'tailLength_MS': 0.11105987611837578,
        #            'MoI0GstdTail_GLS_MS': 0, 'CoM1Reval_15_MS_GLS': 0, 'mGErr_MS': 0.14289875173370326,
        #            'MoI1Geval_15_MS_GLS': 0, 'MoI1Reval_22_GLS_GLS': 0.16525694222789661}
        #         ktot=0.
        #         for pN in PARAM_NAMES:
        #             if pN not in PARAM_NAMES_USE_GLS: KStests[pN+'_GLS']=0
        #             else: ktot+= KStests[pN+'_GLS']
        #             if pN not in PARAM_NAMES_USE_MS: KStests[pN+'_MS']=0
        #             else: ktot+= KStests[pN+'_MS']
        #         ktot/=len(PARAM_NAMES_USE_GLS)+len(PARAM_NAMES_USE_MS)#makes average weight 1
        #         for key in KStests:
        #             KStests[key]/=ktot
        #         if weights is None:
        #             self.weights = KStests

        if weights is None:
            self.weights = {}
            for key in self.controlOrigin:
                self.weights[key] = 1.
        else:
            self.weights = weights

        if origin is None:
            self.origin = self.controlOrigin
        # if origin is None: self.origin = rnaiMean
        else:
            self.origin = origin
        if norms is None:
            self.paramNorms = self.rnaiStd
        # if norms is None: self.paramNorms = controlStd
        else:
            self.paramNorms = norms

    def setNDPosition(self):  # paramsWT and paramsNorm (i.e. std for WT or RNAi) are dictionaries
        if self.origin is None:
            self.setOrigin()

        for i in range(len(PARAM_NAMES)):  # iterates through the index of parameter names to analyze
            pN = PARAM_NAMES[i]  # pN is the string value associated with index i
            if pN in self.paramsUseGLS and pN in PARAM_NAMES_USE_GLS:  # checks to see if parameter name is in paramsUseGLS(good usable parameter RNAiclass dictionary)
                if self.paramNorms[pN + '_GLS'] == 0:
                    self.posND[i] = self.weights[pN + '_GLS'] * (self.paramsGLS[pN] - self.origin[pN + '_GLS'])
                else:
                    self.posND[i] = self.weights[pN + '_GLS'] * (self.paramsGLS[pN] - self.origin[pN + '_GLS']) / \
                                    self.paramNorms[
                                        pN + '_GLS']  # calculates normalized value for each parameter and stores it in array posND
                    self.stdV[i] = self.controlStd[pN + '_GLS'] / self.paramNorms[pN + '_GLS']
            if pN in self.paramsUseMS and pN in PARAM_NAMES_USE_MS:  # checks to see if parameter name is in paramsUseGLS(good usable parameter RNAiclass dictionary)
                if self.paramNorms[pN + '_MS'] == 0:
                    self.posND[i + len(PARAM_NAMES)] = self.weights[pN + '_MS'] * (
                        self.paramsMS[pN] - self.origin[pN + '_MS'])
                else:
                    self.posND[i + len(PARAM_NAMES)] = self.weights[pN + '_MS'] * (
                        self.paramsMS[pN] - self.origin[pN + '_MS']) / self.paramNorms[
                                                           pN + '_MS']  # calculates normalized value for each parameter and stores it in array posND
                    self.stdV[i + len(PARAM_NAMES)] = self.controlStd[pN + '_MS'] / self.paramNorms[pN + '_MS']
        self.NDFlag = True

    #         self.posND[self.stdV>0.5]=0

    def getDims(self):  # returns the total number of considered dimensions
        return len([pN for pN in self.paramsUseGLS if pN in PARAM_NAMES_USE_GLS] +
                   [pN for pN in self.paramsUseMS if pN in PARAM_NAMES_USE_MS])

    def getOverlapDims(self, RNAi):  #
        n = 0  # overlap
        pNGLS = [pN for pN in self.paramsUseGLS if pN in PARAM_NAMES_USE_GLS]
        pNMS = [pN for pN in self.paramsUseMS if pN in PARAM_NAMES_USE_MS]
        pNGLSRNAi = [pN for pN in RNAi.paramsUseGLS if pN in PARAM_NAMES_USE_GLS]
        pNMSRNAi = [pN for pN in RNAi.paramsUseMS if pN in PARAM_NAMES_USE_MS]
        for pN in pNGLS:  # counts the number of overlapping parameters between RNAi_1 and RNAi_2 in GLS
            if pN in pNGLSRNAi:
                n += 1  # adds 1 to the value of n
        for pN in pNMS:  # counts the number of overlapping parameters between RNAi_1 and RNAi_2 in MS
            if pN in pNMSRNAi:
                n += 1  # adds 1 to the value of n
        return n

    def getDist2Zero(self):
        if not self.NDFlag:
            self.setNDPosition()
        return getNDDistance(self.posND, 'control', self.pNamesUse)

    def getDistance(self, RNAi):
        if not self.NDFlag:
            self.setNDPosition()
        if not RNAi.NDFlag:
            RNAi.setNDPosition()
        return getNDDistance(self.posND, RNAi.posND, self.pNamesUse)

    def getMiddle(self, RNAi):
        if self.NDFlag:
            return getNDAvgVector([self.posND, RNAi.posND])
        else:
            printLog('WARNING: origin is not set yet')
            return np.nan

    def getPAD(self, RNAi):
        """
        calculates PAD between two RNAi conditions
        :param RNAi:
        :return:
        """
        dist = self.getDistance(RNAi)
        if np.isinf(dist):
            return 0.
        pad = 1. - dist / (self.getDist2Zero() + RNAi.getDist2Zero())
        if pad > 0:
            return pad
        else:
            return 0

            #     def getPAD(self, RNAi): #calculates PAD between two RNAi conditions
            #         if not self.NDFlag: self.setNDPosition()
            #         if not RNAi.NDFlag: RNAi.setNDPosition()
            #         avgV = embdFunc.getNDAvgVector((np.abs(self.posND),np.abs(RNAi.posND)))
            #         dist = embdFunc.getNDDistance(self.posND/avgV/2, RNAi.posND/avgV/2, self.pNamesUse)
            #         return 1.-dist/np.sqrt(len(PARAM_NAMES_USE_GLS)+len(PARAM_NAMES_USE_MS))


    def get_cosine_dist(self, RNAi):
        """
        takes in two RNAi objects and calculates the cosine distance between the two conditions (vectors in
        n-dimensional space)- cosine of 0 is 1 (maximally similar), a cosine of 90deg is 0.
        :param RNAi: RNAi object
        :return: cosine distance
        """
        v1 = self.posND
        v2 = RNAi.posND

        N = len(self.pNamesUse[0]) + len(self.pNamesUse[1])  # total number of parameters considered GLS + MS
        vOverlap = np.invert(np.isnan(
            v1 * v2))  # overlapping vector, checks whether a parameter (coordinate) is not nan for each considered group (RNAi)
        if sum(vOverlap) == 0:
            return np.inf  # if there are no overlapping coordinates, returns infinity
        dotp_overlap = np.dot(v1[vOverlap],v2[vOverlap])  # finds the dot product for overlapping coordinates within the
        # vector (this is equiv. to x1*x2 + y1*y2 + z1*z2...)
        i1 = np.sum(np.invert(np.isnan(v1))) - np.sum(
            vOverlap)  # sums all coordinates that are not equal to nan for v1 and subtracts overlap
        i2 = np.sum(np.invert(np.isnan(v2))) - np.sum(
            vOverlap)  # sums all coordinates that are not equal to nan for v2 and subtracts overlap
        weight = 1/ (1. * N / (N - i1 - i2))  # uses same weighting as used for distance calculation, but inverts it
        # since higher values have opposite meaning in euc distance vs. cosine dis
        print("i1 = {0}, i2 = {1}, overlap = {2}, weight = {3}".format(i1, i2, np.sum(vOverlap), weight))
        # cosine_dist = dotp_overlap / (self.getDist2Zero() * RNAi.getDist2Zero())
        cosine_dist = dotp_overlap / (np.linalg.norm(v1[vOverlap]) * np.linalg.norm(v2[vOverlap]))

        adjusted_cos_dist = cosine_dist * weight
        if adjusted_cos_dist > 0:
            return adjusted_cos_dist
        else:
            return 0

    def get_norm_dist(self, r2):
        """
        calculates PAD and normalized distance between two RNAi conditions
        :param r2: RNAiClass object
        :return: PAD and distance to r2 corrected by number of dimensions
        """
        return np.exp(-self.getDistance(r2) / len(PARAM_NAMES_USE_GLS + PARAM_NAMES_USE_MS))

    def getSimParams(self, RNAi):
        if not self.NDFlag: self.setNDPosition()
        if not RNAi.NDFlag: RNAi.setNDPosition()
        diff = np.abs(self.posND - RNAi.posND)
        inds = np.argsort(diff)
        names = np.array([pN + '_GLS' for pN in PARAM_NAMES] + [pN + '_MS' for pN in PARAM_NAMES])
        names = names[inds]
        posOrdered = self.posND[inds]
        bestPars, bestParVals = [], []
        i = 0
        while len(bestPars) < 10 and i < names.size:
            pN = names[i]
            if pN[-3:] == 'GLS' and pN[:-4] in PARAM_NAMES_USE_GLS:
                bestPars.append(pN)
                bestParVals.append(posOrdered[i])
            elif pN[-3:] == '_MS' and pN[:-3] in PARAM_NAMES_USE_MS:
                bestPars.append(pN)
                bestParVals.append(posOrdered[i])
            i += 1
        return bestPars, bestParVals

    def getDiffParams(self, RNAi):
        if not self.NDFlag: self.setNDPosition()
        if not RNAi.NDFlag: RNAi.setNDPosition()
        diff = np.abs(self.posND - RNAi.posND)
        inds = np.argsort(diff)[::-1]
        names = np.array([pN + '_GLS' for pN in PARAM_NAMES] + [pN + '_MS' for pN in PARAM_NAMES])
        names = names[inds]
        posOrdered = self.posND[inds]
        worstPars, worstParVals = [], []
        i = 0
        while len(worstPars) < 10 and i < names.size:
            pN = names[i]
            if pN[-3:] == 'GLS' and pN[:-4] in PARAM_NAMES_USE_GLS:
                worstPars.append(pN)
                worstParVals.append(posOrdered[i])
            elif pN[-3:] == '_MS' and pN[:-3] in PARAM_NAMES_USE_MS:
                worstPars.append(pN)
                worstParVals.append(posOrdered[i])
            i += 1
        return worstPars, worstParVals

    def getPrev(self):
        printLog('GLS: {0} out of {1} not moving'.format(len(self.GLSNoMoveEmbs),
                                                         len(self.GLSNoMoveEmbs) + len(self.GLSMoveEmbs)))
        printLog('MS: {0} out of {1} not moving'.format(len(self.MSNoMoveEmbs),
                                                        len(self.MSNoMoveEmbs) + len(self.MSMoveEmbs)))
        return len(self.GLSNoMoveEmbs), len(self.MSNoMoveEmbs)

    def getEmbsUsed(self):
        self.refresh_params_data()
        MSembLabels = ''
        GLSembLabels = ''
        for emb in self.MSUseEmbs:
            MSembLabels += emb.label + ','
        print 'MS embs = {0}'.format(MSembLabels[:-1])
        for emb in self.GLSUseEmbs:
            GLSembLabels += emb.label + ','
        print 'GLS embs = {0}'.format(GLSembLabels[:-1])

    def getRepEmb(self, RNAi=None):
        '''
        finds the representative embryo for a seed (self) condition. Finds the EMBRYO VECTOR that is closest to the seed
        RNAi VECTOR (ND position). If RNAi argument is enabled, it finds the closest seed embryo to the queried RNAi ND
        position. I.e. if seed is EMBD0126 and RNAi is EMBD0133, it finds the best example embryo from EMBD0126 that
        resembles the representative embryo for EMBD0133.
        :param RNAi: seed, test RNAi object
        :return: nearest GLS embryo object, nearest MS embryo object
        '''
        if RNAi is None: RNAi = self
        ''' returns one GLS and one MS embryo objects that represent the rnai the best '''
        if len(self.MSUseEmbs) == 0: self.refresh_params_data()
        if len(RNAi.MSUseEmbs) == 0: RNAi.refresh_params_data()

        if not self.NDFlag: self.setNDPosition()
        if not RNAi.NDFlag: RNAi.setNDPosition()
        # calculate distances for each GLS embryo to RNAi
        distGLS = []
        for emb in self.GLSUseEmbs:
            vR = RNAi.posND[:len(PARAM_NAMES)]  # vector RNAi
            vE = np.ones_like(vR) * np.nan     #  vector embryo
            for i in range(len(PARAM_NAMES)):  # iterates through the index of parameter names to analyze
                pN = PARAM_NAMES[i]  # pN is the string value associated with index i
                if pN in self.paramsUseGLS and pN in PARAM_NAMES_USE_GLS:  # checks to see if parameter name is in paramsUseGLS(good usable parameter RNAiclass dictionary)
                    if self.paramNorms[pN + '_GLS'] == 0:
                        vE[i] = self.weights[pN + '_GLS'] * (emb.params[pN] - self.origin[pN + '_GLS'])
                    else:
                        vE[i] = self.weights[pN + '_GLS'] * (emb.params[pN] - self.origin[pN + '_GLS']) / \
                                self.paramNorms[
                                    pN + '_GLS']  # calculates normalized value for each parameter and stores it in array posND
            distGLS.append(getNDDistance(vR, vE, (self.pNamesUse[0], [])))

        # calculate distances for each MS embryo to RNAi
        distMS = []
        for emb in self.MSUseEmbs:
            vR = RNAi.posND[len(PARAM_NAMES):]
            vE = np.ones_like(vR) * np.nan
            for i in range(len(PARAM_NAMES)):  # iterates through the index of parameter names to analyze
                pN = PARAM_NAMES[i]  # pN is the string value associated with index i
                if pN in self.paramsUseMS and pN in PARAM_NAMES_USE_MS:  # checks to see if parameter name is in paramsUseGLS(good usable parameter RNAiclass dictionary)
                    if self.paramNorms[pN + '_MS'] == 0:
                        vE[i] = self.weights[pN + '_MS'] * (emb.params[pN] - self.origin[pN + '_MS'])
                    else:
                        vE[i] = self.weights[pN + '_MS'] * (emb.params[pN] - self.origin[pN + '_MS']) / self.paramNorms[
                            pN + '_MS']  # calculates normalized value for each parameter and stores it in array posND
            distMS.append(getNDDistance(vR, vE, ([], self.pNamesUse[1])))
        return self.GLSUseEmbs[np.argmin(distGLS)], self.MSUseEmbs[np.argmin(distMS)]

    def plotSigmFitAll(self):
        fig1 = myFigure()  # creates a blank fig
        for emb in self.GLSUseEmbs[:4]:  # iterates through embryo object in GLS folder
            fig1 = emb.showSigmFit(fig=fig1, chs=[0, 1,
                                                  2])  # plot asymmetric sigmoidal for spots. assign it to fig so all embs will be plotted on the same plot. lead with emb. to refer to function showsigmfit within emb object within Embryos class
        fig1.title('GLS sigma fit')
        #         fig1.setBGColor('k')
        #         fig1.save(FOLDER_IN+'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/GLS_SpotsSigmaFit_All_{0}.svg'.format(self.label))
        fig2 = myFigure()
        for emb in self.MSUseEmbs[:4]:  # iterates through embryo object in MS folder
            fig2 = emb.showSigmFit(fig=fig2, chs=[0, 1])  # plot asymmetric sigmoidal for total intensity.
        x, y, yerr = self.MSUseEmbs[0].avgCurves['tIntG']
        fig2.errorbar(x[::10], y[::10], yerr[::10], color='g')
        # fig2.errorbar(x[::10], y[::10], yerr[::10], color='greenyellow')

        x, y, yerr = self.MSUseEmbs[0].avgCurves['tIntR']
        # fig2.errorbar(x[::10], y[::10], yerr[::10], color='mistyrose')
        fig2.errorbar(x[::10], y[::10], yerr[::10], color='r')

        fig2.title('MS sigma fit')

        # fig2.setBGColor('k')
        fig2.save(FOLDER_IN+'Automated_analysis/Plots_Images/MS_SigmaFit_All_{0}.svg'.format(self.label))
        # fig1.setBGColor('k')
        fig1.save(FOLDER_IN+'Automated_analysis/Plots_Images/GLS_SigmaFit_All_{0}.svg'.format(self.label))
        return fig1, fig2

    def plotSigmFitAvg(self):
        def sigmoidalGLS(x, a, m, s, r):
            return a * (1. - (1. + np.exp((x - m) / s)) ** (-r))

        def sigmoidalMS(x, a, b, m, s, r):
            return a - a * (1. + np.exp((x - m) / s)) ** (-r)

        avgR = AvgCurve()
        avgG = AvgCurve()
        avgY = AvgCurve()
        for emb in self.GLSUseEmbs:  # iterates through embryo object in GLS folder
            printLog('{0} time alignment = {1}'.format(emb.label, emb.time))
            if emb.t0 != 0:
                x = emb.time - emb.t0
                avgG.add(x, emb.green)
                avgR.add(x, emb.red)
                avgY.add(x, emb.yellow)
        fig1 = myFigure()  # creates a blank fig
        x, y, yerr = avgG.getAvg(2, every=5)
        fig1.errorbar(x, y, yerr, color='green')
        fig1.plot(x, sigmoidalGLS(x, self.paramsGLS['aG'], self.paramsGLS['mG'], self.paramsGLS['sG'], 1), color='w')
        x, y, yerr = avgR.getAvg(2, every=5)
        fig1.errorbar(x, y, yerr, color='red')
        fig1.plot(x, sigmoidalGLS(x, self.paramsGLS['aR'], self.paramsGLS['mR'], self.paramsGLS['sR'], 1), color='w')
        x, y, yerr = avgY.getAvg(2, every=5)
        fig1.errorbar(x, y, yerr, color='yellow')
        fig1.plot(x, sigmoidalGLS(x, self.paramsGLS['aY'], self.paramsGLS['mY'], self.paramsGLS['sY'], 1), color='w')
        fig1.title('GLS sigma fit')
        fig1.xlim([-5, 20])
        fig1.ylim([-50, 350])
        #         fig1.setBGColor('k')
        #         fig1.save(FOLDER_IN+'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/GLS_SpotsSigmaFit_Avg_{0}.svg'.format(self.label))

        avgR = AvgCurve()
        avgG = AvgCurve()
        for emb in self.MSUseEmbs:  # iterates through embryo object in GLS folder
            if emb.t0 != 0:
                x = emb.time - emb.t0
                avgG.add(x, emb.tInt[0])
                avgR.add(x, emb.tInt[1])
        fig1 = myFigure()  # creates a blank fig
        x, y, yerr = avgG.getAvg(2, every=5)
        fig1.errorbar(x, y, yerr, color='green')
        fig1.plot(x, sigmoidalMS(x, self.paramsMS['aG'], self.paramsMS['bG'], self.paramsGLS['mG'], self.paramsMS['sG'],
                                 1), color='g')
        x, y, yerr = avgR.getAvg(2, every=5)
        fig1.errorbar(x, y, yerr, color='red')
        fig1.plot(x, sigmoidalMS(x, self.paramsMS['aR'], self.paramsMS['bR'], self.paramsGLS['mR'], self.paramsMS['sR'],
                                 1), color='r')
        x, y, yerr = self.MSUseEmbs[0].avgCurves['tIntG']
        fig1.errorbar(x[::10], y[::10], yerr[::10], color='greenyellow')
        x, y, yerr = self.MSUseEmbs[0].avgCurves['tIntR']
        fig1.errorbar(x[::10], y[::10], yerr[::10], color='mistyrose')
        fig1.title('MS sigma fit')
        fig1.xlim([-15, 15])

    #         fig1.setBGColor('k')
    #         fig1.save(FOLDER_IN+'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/MS_totIntSigmaFit_Avg_{0}.svg'.format(self.label))
        fig1.show()

    def plotCentMassAll(self):
        fig1_0 = myFigure()
        fig1_1 = myFigure()
        fig2_0 = myFigure()
        fig2_1 = myFigure()
        for emb in self.GLSUseEmbs[:5]:
            fig1_0 = emb.showDistCentMass(0, fig=fig1_0, chs=[1, 2])
            fig1_1 = emb.showDistCentMass(1, fig=fig1_1, chs=[1, 2])
        for emb in self.MSUseEmbs:
            fig2_0 = emb.showDistCentMass(0, fig=fig2_0, chs=[1])
            fig2_1 = emb.showDistCentMass(1, fig=fig2_1, chs=[1])
        fig1_0.title('center of mass- GLS end-on')
        fig1_0.ylim([0, 1])
        x, y, yerr = self.GLSUseEmbs[0].avgCurves['CoM0R']
        fig1_0.errorbar(x[::10], y[::10], yerr[::10], color='mistyrose')
        x, y, yerr = self.GLSUseEmbs[0].avgCurves['CoM0Y']
        fig1_0.errorbar(x[::10], y[::10], yerr[::10], color='gold')
        #         fig1_0.setBGColor('k')
        #         fig1_0.save(FOLDER_IN+'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/CoM_All_GLS_endOn_{0}.svg'.format(self.label))
        fig1_1.title('center of mass- GLS AP')
        fig1_1.ylim([0, 1])
        x, y, yerr = self.GLSUseEmbs[0].avgCurves['CoM1R']
        fig1_1.errorbar(x[::10], y[::10], yerr[::10], color='mistyrose')
        x, y, yerr = self.GLSUseEmbs[0].avgCurves['CoM1Y']
        fig1_1.errorbar(x[::10], y[::10], yerr[::10], color='gold')
        #         fig1_1.setBGColor('k')
        #         fig1_1.save(FOLDER_IN+'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/CoM_All_GLS_AP_{0}.svg'.format(self.label))
        fig2_0.title('center of mass- MS end-on')
        fig2_0.ylim([0, 1])
        #         fig2_0.setBGColor('k')
        x, y, yerr = self.MSUseEmbs[0].avgCurves['CoM0']
        fig2_0.errorbar(x[::10], y[::10], yerr[::10], color='mistyrose')
        #         fig2_0.save(FOLDER_IN+'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/CoM_All_MS_endOn_{0}.svg'.format(self.label))
        fig2_1.title('center of mass- MS AP')
        fig2_1.ylim([0, 1])
        #         fig2_1.setBGColor('k')
        x, y, yerr = self.MSUseEmbs[0].avgCurves['CoM1']
        fig2_1.errorbar(x[::10], y[::10], yerr[::10], color='mistyrose')

    #         fig2_1.save(FOLDER_IN+'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/CoM_All_MS_AP_{0}.svg'.format(self.label))

    def plotCentMassAvg(self):
        fig1_0 = myFigure()
        fig1_1 = myFigure()
        fig2_0 = myFigure()
        fig2_1 = myFigure()
        avgCM0R = AvgCurve()
        avgCM1R = AvgCurve()
        avgCM0Y = AvgCurve()
        avgCM1Y = AvgCurve()
        for emb in self.GLSUseEmbs[:1]:
            x = np.arange(31) - emb.t0
            y = np.array([emb.getDistFromGreen(1, slide, 0) for slide in range(31)])
            avgCM0R.add(x, y)
            dist1 = np.array([emb.getDistFromGreen(1, slide, 1) for slide in range(31)])
            dist2 = np.array([emb.getDistFromGreen(1, slide, 2) for slide in range(31)])
            y = (dist1 + dist2) / 2.
            avgCM1R.add(x, y)
            y = np.array([emb.getDistFromGreen(2, slide, 0) for slide in range(31)])
            avgCM0Y.add(x, y)
            dist1 = np.array([emb.getDistFromGreen(2, slide, 1) for slide in range(31)])
            dist2 = np.array([emb.getDistFromGreen(2, slide, 2) for slide in range(31)])
            y = (dist1 + dist2) / 2.
            avgCM1Y.add(x, y)

        avgCM0R_MS = AvgCurve()
        avgCM1R_MS = AvgCurve()
        for emb in self.MSUseEmbs:
            x = np.arange(31) - emb.t0
            y = np.array([emb.getDistFromGreen(1, slide, 0) for slide in range(31)])
            avgCM0R_MS.add(x, y)
            dist1 = np.array([emb.getDistFromGreen(1, slide, 1) for slide in range(31)])
            dist2 = np.array([emb.getDistFromGreen(1, slide, 2) for slide in range(31)])
            y = (dist1 + dist2) / 2.
            avgCM1R_MS.add(x, y)

        x, y, yerr = avgCM0R.getAvg(2, every=5)
        fig1_0.errorbar(x, y, yerr, color='red')
        x, y, yerr = avgCM0Y.getAvg(2, every=5)
        fig1_0.errorbar(x, y, yerr, color='yellow')
        #         x = [self.paramsGLS['mR']+self.paramsGLS['sR'], self.paramsGLS['mR']+4*self.paramsGLS['sR'], self.paramsGLS['mR']+6*self.paramsGLS['sR']]
        #         y = [self.paramsGLS['R2Gax0v1'], self.paramsGLS['R2Gax0v2'], self.paramsGLS['R2Gax0v3']]
        #         fig1_0.scatter(x, y, color='k')
        #         x = [self.paramsGLS['mY']+self.paramsGLS['sY'], self.paramsGLS['mY']+4*self.paramsGLS['sY'], self.paramsGLS['mY']+6*self.paramsGLS['sY']]
        #         y = [self.paramsGLS['Y2Gax0v1'], self.paramsGLS['Y2Gax0v2'], self.paramsGLS['Y2Gax0v3']]
        #         fig1_0.scatter(x, y, color='k')
        fig1_0.title('center of mass- GLS end-on')
        fig1_0.ylim([0, 1])
        fig1_0.xlim([-5, 20])
        #         fig1_0.setBGColor('k')
        fig1_0.save(FOLDER_IN + 'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/CoM_Avg_GLS_endOn_{0}.svg'.format(
            self.label))

        x, y, yerr = avgCM1R.getAvg(2, every=5)
        fig1_1.errorbar(x, y, yerr, color='red')
        x, y, yerr = avgCM1Y.getAvg(2, every=5)
        fig1_1.errorbar(x, y, yerr, color='yellow')
        #         x = [self.paramsGLS['mR']+self.paramsGLS['sR'], self.paramsGLS['mR']+4*self.paramsGLS['sR'], self.paramsGLS['mR']+6*self.paramsGLS['sR']]
        #         y = [self.paramsGLS['R2Gax1v1'], self.paramsGLS['R2Gax1v2'], self.paramsGLS['R2Gax1v3']]
        #         fig1_0.scatter(x, y, color='k')
        #         x = [self.paramsGLS['mY']+self.paramsGLS['sY'], self.paramsGLS['mY']+4*self.paramsGLS['sY'], self.paramsGLS['mY']+6*self.paramsGLS['sY']]
        #         y = [self.paramsGLS['Y2Gax1v1'], self.paramsGLS['Y2Gax1v2'], self.paramsGLS['Y2Gax1v3']]
        #         fig1_0.scatter(x, y, color='k')
        fig1_1.title('center of mass- GLS AP')
        fig1_1.ylim([0, 1])
        fig1_1.xlim([-5, 20])
        #         fig1_1.setBGColor('k')
        fig1_1.save(
            FOLDER_IN + 'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/CoM_Avg_GLS_AP_{0}.svg'.format(self.label))

        x, y, yerr = avgCM0R_MS.getAvg(2, every=5)
        fig2_0.errorbar(x, y, yerr, color='red')
        #         x = [self.paramsMS['mR'], self.paramsMS['mR']+self.paramsMS['sR'], self.paramsMS['mR']+2*self.paramsMS['sR']]
        #         y = [self.paramsMS['R2Gax0v1'], self.paramsMS['R2Gax0v2'], self.paramsMS['R2Gax0v3']]
        #         fig1_0.scatter(x, y, color='k')
        fig2_0.title('center of mass- MS end-on')
        fig2_0.ylim([0, 1])
        fig2_0.xlim([-15, 15])
        #         fig2_0.setBGColor('k')
        x, y, yerr = self.MSUseEmbs[0].avgCurves['CoM0']
        fig2_0.errorbar(x[::10], y[::10], yerr[::10], color='mistyrose')
        #         fig2_0.save(FOLDER_IN+'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/CoM_Avg_MS_endOn_{0}.svg'.format(self.label))

        x, y, yerr = avgCM1R_MS.getAvg(2, every=5)
        fig2_1.errorbar(x, y, yerr, color='red')
        #         x = [self.paramsMS['mR'], self.paramsMS['mR']+self.paramsMS['sR'], self.paramsMS['mR']+2*self.paramsMS['sR']]
        #         y = [self.paramsMS['R2Gax1v1'], self.paramsMS['R2Gax1v2'], self.paramsMS['R2Gax1v3']]
        #         fig1_0.scatter(x, y, color='k')
        fig2_1.title('center of mass- MS AP')
        fig2_1.ylim([0, 1])
        fig2_1.xlim([-15, 15])
        #         fig2_1.setBGColor('k')
        x, y, yerr = self.MSUseEmbs[0].avgCurves['CoM1']
        fig2_1.errorbar(x[::10], y[::10], yerr[::10], color='mistyrose')

    #         fig2_1.save(FOLDER_IN+'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/CoM_Avg_MS_AP_{0}.svg'.format(self.label))

    def plotMOIAll(self):
        fig1_0 = myFigure()
        fig1_1 = myFigure()
        fig2_0 = myFigure()
        fig2_1 = myFigure()
        for emb in self.GLSUseEmbs:
            fig1_0 = emb.showMOI(0, fig=fig1_0, chs=range(3),
                                 setColor=True)  # plots MOI for GLS. If set color arg is False, uses autocolor (would use this if plotting lots of embryos on one plot)
            fig1_1 = emb.showMOI(1, fig=fig1_1, chs=range(3), setColor=True)
        for emb in self.MSUseEmbs:
            fig2_0 = emb.showMOI(0, fig=fig2_0, chs=range(2),
                                 setColor=True)  # plots MOI for GLS. If set color arg is False, uses autocolor (would use this if plotting lots of embryos on one plot)
            fig2_1 = emb.showMOI(1, fig=fig2_1, chs=range(2), setColor=True)
        fig1_0.title('Moment of Inertia- GLS end-on')
        #         x, y, yerr = self.GLSUseEmbs[0].avgCurves['MoI0G']
        #         fig1_0.errorbar(x[::10], y[::10], yerr[::10], color = 'springgreen')
        #         x, y, yerr = self.GLSUseEmbs[0].avgCurves['MoI0R']
        #         fig1_0.errorbar(x[::10], y[::10], yerr[::10], color = 'mistyrose')
        #         x, y, yerr = self.GLSUseEmbs[0].avgCurves['MoI0Y']
        #         fig1_0.errorbar(x[::10], y[::10], yerr[::10], color = 'gold')
        fig1_0.legend(4)
        #         fig1_0.setBGColor('k')
        #         fig1_0.save(FOLDER_IN+'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/MoI_All_GLS_endOn_{0}.svg'.format(self.label))
        fig1_1.title('Moment of Inertia- GLS AP')
        #         x, y, yerr = self.GLSUseEmbs[0].avgCurves['MoI1G']
        #         fig1_1.errorbar(x[::10], y[::10], yerr[::10], color = 'springgreen')
        #         x, y, yerr = self.GLSUseEmbs[0].avgCurves['MoI1R']
        #         fig1_1.errorbar(x[::10], y[::10], yerr[::10], color = 'mistyrose')
        #         x, y, yerr = self.GLSUseEmbs[0].avgCurves['MoI1Y']
        #         fig1_1.errorbar(x[::10], y[::10], yerr[::10], color = 'gold')
        #         fig1_1.setBGColor('k')
        fig1_1.legend(4)
        #         fig1_1.save(FOLDER_IN+'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/MoI_All_GLS_AP_{0}.svg'.format(self.label))
        fig2_0.title('Moment of Inertia- MS end-on')
        #         fig2_0.setBGColor('k')
        #         x, y, yerr = self.MSUseEmbs[0].avgCurves['MoI0G']
        #         fig2_0.errorbar(x[::10], y[::10], yerr[::10], color = 'greenyellow')
        #         x, y, yerr = self.MSUseEmbs[0].avgCurves['MoI0R']
        #         fig2_0.errorbar(x[::10], y[::10], yerr[::10], color = 'mistyrose')
        fig2_0.xlim((-5, 30))
        fig2_0.legend(4)
        #         fig2_0.save(FOLDER_IN+'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/MoI_All_MS_endOn_{0}.svg'.format(self.label))
        fig2_1.title('Moment of Inertia- MS AP')
        #         fig2_1.setBGColor('k')
        #         x, y, yerr = self.MSUseEmbs[0].avgCurves['MoI1G']
        #         fig2_1.errorbar(x[::10], y[::10], yerr[::10], color = 'greenyellow')
        #         x, y, yerr = self.MSUseEmbs[0].avgCurves['MoI1R']
        #         fig2_1.errorbar(x[::10], y[::10], yerr[::10], color = 'mistyrose')
        fig2_1.legend(4)
        #         fig2_1.save(FOLDER_IN+'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/MoI_All_MS_AP_{0}.svg'.format(self.label))
        fig2_1.show()
        fig2_0.show()

        # fig2_0.close()
        # fig2_1.close()

        fig1_0.show()
        fig1_1.show()

    def plotMOIAvg(self):  # FIXME emb.align no longer exists
        fig1_0 = myFigure()
        fig1_1 = myFigure()
        fig2_0 = myFigure()
        fig2_1 = myFigure()
        avgMI0G_GLS = AvgCurve()
        avgMI1G_GLS = AvgCurve()
        avgMI0R_GLS = AvgCurve()
        avgMI1R_GLS = AvgCurve()
        avgMI0Y_GLS = AvgCurve()
        avgMI1Y_GLS = AvgCurve()
        for emb in self.GLSUseEmbs:
            x = np.arange(31) - emb.t0
            y = emb.getAllMoments(0, 0)  # retrieves MOI values for axis 0
            avgMI0G_GLS.add(x, y)
            MOI1 = emb.getAllMoments(0, 1)
            MOI2 = emb.getAllMoments(0, 2)
            y = (MOI1 + MOI2) / 2.  # averages MOI values for axis 1,2
            avgMI1G_GLS.add(x, y)

            y = emb.getAllMoments(1, 0)  # retrieves MOI values for axis 0
            avgMI0R_GLS.add(x, y)
            MOI1 = emb.getAllMoments(1, 1)
            MOI2 = emb.getAllMoments(1, 2)
            y = (MOI1 + MOI2) / 2.  # averages MOI values for axis 1,2
            avgMI1R_GLS.add(x, y)

            y = emb.getAllMoments(2, 0)  # retrieves MOI values for axis 0
            avgMI0Y_GLS.add(x, y)
            MOI1 = emb.getAllMoments(2, 1)
            MOI2 = emb.getAllMoments(2, 2)
            y = (MOI1 + MOI2) / 2.  # averages MOI values for axis 1,2
            avgMI1Y_GLS.add(x, y)

        avgMI0G_MS = AvgCurve()
        avgMI1G_MS = AvgCurve()
        avgMI0R_MS = AvgCurve()
        avgMI1R_MS = AvgCurve()
        for emb in self.MSUseEmbs:
            x = np.arange(31) - emb.t0
            y = emb.getAllMoments(0, 0)  # retrieves MOI values for axis 0
            avgMI0G_MS.add(x, y)
            MOI1 = emb.getAllMoments(0, 1)
            MOI2 = emb.getAllMoments(0, 2)
            y = (MOI1 + MOI2) / 2.  # averages MOI values for axis 1,2
            avgMI1G_MS.add(x, y)

            y = emb.getAllMoments(1, 0)  # retrieves MOI values for axis 0
            avgMI0R_MS.add(x, y)
            MOI1 = emb.getAllMoments(1, 1)
            MOI2 = emb.getAllMoments(1, 2)
            y = (MOI1 + MOI2) / 2.  # averages MOI values for axis 1,2
            avgMI1R_MS.add(x, y)

        x, y, yerr = avgMI0G_GLS.getAvg(2, every=5)
        fig1_0.errorbar(x, y, yerr, color='green')
        x, y, yerr = avgMI0R_GLS.getAvg(2, every=5)
        fig1_0.errorbar(x, y, yerr, color='red')
        x, y, yerr = avgMI0Y_GLS.getAvg(2, every=5)
        fig1_0.errorbar(x, y, yerr, color='yellow')
        #         x = [self.paramsGLS['mG']+self.paramsGLS['sG'], self.paramsGLS['mG']+4*self.paramsGLS['sG'], self.paramsGLS['mG']+6*self.paramsGLS['sG']]
        #         y = [self.paramsGLS['miGv1a0'], self.paramsGLS['miGv2a0'], self.paramsGLS['miGv3a0']]
        #         fig1_0.scatter(x-self.paramsGLS['mG'], y, color='k')
        #         x = [self.paramsGLS['mR']+self.paramsGLS['sR'], self.paramsGLS['mR']+4*self.paramsGLS['sR'], self.paramsGLS['mR']+6*self.paramsGLS['sR']]
        #         y = [self.paramsGLS['miRv1a0'], self.paramsGLS['miRv2a0'], self.paramsGLS['miRv3a0']]
        #         fig1_0.scatter(x-self.paramsGLS['mG'], y, color='k')
        #         x = [self.paramsGLS['mY']+self.paramsGLS['sY'], self.paramsGLS['mY']+4*self.paramsGLS['sY'], self.paramsGLS['mY']+6*self.paramsGLS['sY']]
        #         y = [self.paramsGLS['miYv1a0'], self.paramsGLS['miYv2a0'], self.paramsGLS['miYv3a0']]
        #         fig1_0.scatter(x-self.paramsGLS['mG'], y, color='k')
        fig1_0.title('Moment of Inertia- GLS end-on')
        fig1_0.ylim([0, 1])
        fig1_0.xlim([-5, 20])
        #         fig1_0.setBGColor('k')
        fig1_0.save(FOLDER_IN + 'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/MoI_Avg_GLS_endOn_{0}.svg'.format(
            self.label))

        x, y, yerr = avgMI1G_GLS.getAvg(2, every=5)
        fig1_1.errorbar(x, y, yerr, color='green')
        x, y, yerr = avgMI1R_GLS.getAvg(2, every=5)
        fig1_1.errorbar(x, y, yerr, color='red')
        x, y, yerr = avgMI1Y_GLS.getAvg(2, every=5)
        fig1_1.errorbar(x, y, yerr, color='yellow')
        #         x = [self.paramsGLS['mG']+self.paramsGLS['sG'], self.paramsGLS['mG']+4*self.paramsGLS['sG'], self.paramsGLS['mG']+6*self.paramsGLS['sG']]
        #         y = [self.paramsGLS['miGv1a1'], self.paramsGLS['miGv2a1'], self.paramsGLS['miGv3a1']]
        #         fig1_1.scatter(x-self.paramsGLS['mG'], y, color='k')
        #         x = [self.paramsGLS['mR']+self.paramsGLS['sR'], self.paramsGLS['mR']+4*self.paramsGLS['sR'], self.paramsGLS['mR']+6*self.paramsGLS['sR']]
        #         y = [self.paramsGLS['miRv1a1'], self.paramsGLS['miRv2a1'], self.paramsGLS['miRv3a1']]
        #         fig1_1.scatter(x-self.paramsGLS['mG'], y, color='k')
        #         x = [self.paramsGLS['mY']+self.paramsGLS['sY'], self.paramsGLS['mY']+4*self.paramsGLS['sY'], self.paramsGLS['mY']+6*self.paramsGLS['sY']]
        #         y = [self.paramsGLS['miYv1a1'], self.paramsGLS['miYv2a1'], self.paramsGLS['miYv3a1']]
        #         fig1_0.scatter(x-self.paramsGLS['mG'], y, color='k')
        fig1_1.title('Moment of Inertia- GLS AP')
        fig1_1.ylim([0, 1])
        fig1_1.xlim([-5, 20])
        #         fig1_1.setBGColor('k')
        fig1_1.save(
            FOLDER_IN + 'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/MoI_Avg_GLS_AP_{0}.svg'.format(self.label))

        x, y, yerr = avgMI0G_MS.getAvg(2, every=5)
        fig2_0.errorbar(x, y, yerr, color='green')
        x, y, yerr = avgMI0R_MS.getAvg(2, every=5)
        fig2_0.errorbar(x, y, yerr, color='red')
        #         x = [self.paramsMS['mG'], self.paramsMS['mG']+self.paramsMS['sG'], self.paramsMS['mG']+2*self.paramsMS['sG']]
        #         y = [self.paramsMS['miGv1a0'], self.paramsMS['miGv2a0'], self.paramsMS['miGv3a0']]
        #         fig2_0.scatter(x-self.paramsMS['mG'], y, color='k')
        #         x = [self.paramsMS['mR'], self.paramsMS['mR']+self.paramsMS['sR'], self.paramsMS['mR']+2*self.paramsMS['sR']]
        #         y = [self.paramsMS['miRv1a0'], self.paramsMS['miRv2a0'], self.paramsMS['miRv3a0']]
        #         fig2_0.scatter(x-self.paramsMS['mG'], y, color='k')
        fig2_0.title('Moment of Inertia- MS end-on')
        fig2_0.ylim([0, 1])
        fig2_0.xlim([-15, 15])
        #         fig2_0.setBGColor('k')
        fig2_0.save(FOLDER_IN + 'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/MoI_Avg_MS_endOn_{0}.svg'.format(
            self.label))

        x, y, yerr = avgMI1G_MS.getAvg(2, every=5)
        fig2_1.errorbar(x, y, yerr, color='green')
        x, y, yerr = avgMI1R_MS.getAvg(2, every=5)
        fig2_1.errorbar(x, y, yerr, color='red')
        #         x = [self.paramsMS['mG'], self.paramsMS['mG']+self.paramsMS['sG'], self.paramsMS['mG']+2*self.paramsMS['sG']]
        #         y = [self.paramsMS['miGv1a1'], self.paramsMS['miGv2a1'], self.paramsMS['miGv3a1']]
        #         fig2_1.scatter(x-self.paramsMS['mG'], y, color='k')
        #         x = [self.paramsMS['mR'], self.paramsMS['mR']+self.paramsMS['sR'], self.paramsMS['mR']+2*self.paramsMS['sR']]
        #         y = [self.paramsMS['miRv1a1'], self.paramsMS['miRv2a1'], self.paramsMS['miRv3a1']]
        #         fig2_1.scatter(x-self.paramsMS['mG'], y, color='k')
        fig2_1.title('Moment of Inertia- MS AP')
        fig2_1.ylim([0, 1])
        fig2_1.xlim([-15, 15])
        #         fig2_1.setBGColor('k')
        fig2_1.save(
            FOLDER_IN + 'Automated_analysis/Plots_Images/4LabMeeting_Aug2016/MoI_Avg_MS_AP_{0}.svg'.format(self.label))

    def showMaxProjMS(self, t):
        self.getEmbryos()
        self.asignUse()
        #         fig = plt.figure()
        for emb in self.MSUseEmbs[:3]:
            figE = emb.showMaxProj(t)
            #             fig.axes.append(figE.ax)
            #         return fig

    def showMaxProjGLS(self, t):
        self.getEmbryos()
        self.asignUse()
        for emb in self.GLSUseEmbs[:3]:
            figE = emb.showMaxProj(t)

    def makeTalignMovies(self):  # , fileName
        '''makes composite movie of all embryos considered (GLS and MS) for a given RNAi. Representative embryo for each strain is marked with italics'''
        from fileLookup import FILENAME_WEB_CONTENT

        space = 0  # variable to provide spacing between images

        self.getEmbryos()
        self.asignUse()
        allE = [self.GLSUseEmbs, self.MSUseEmbs]

        repEmbG, repEmbM = self.getRepEmb()  # a list of representative embryo object [GLS,MS]
        repEmbGLS = repEmbG.folder.split('/')[4]
        repEmbMS = repEmbM.folder.split('/')[4]
        for i in allE:  # iterates through GLS embs, then through MS embs
            #             for emb in self.GLSUseEmbs[:2]:#[:1]
            maxProjEmbs = []
            shape = []
            name = []
            for emb in i:  # [:1]
                embProjTime = []  # projections for each timepoint for a given embryo
                name.append(emb.folder.split('/')[4])
                for t in range(-2,
                               27):  # goes through all timepoints-- here time has been converted  (nearest slide used for dimensionless time)- displays time points starting 3 back from t0.
                    if emb.strain == 'MS':
                        embProj, (xG, yG), (xR, yR) = emb.getMaxProj(
                            t - DT0_GLS_2_MS)  # gets maximum projection for the slide (t). To display similar stages, we shift t by delta t0 value for MS
                    else:
                        embProj, (xG, yG), (xR, yR) = emb.getMaxProj(t)  # gets maximum projection for the slide (t)
                    embProj = (255 * embProj).astype(np.uint8)  # converts to 8bit

                    embProjIm = cv2.cvtColor(embProj, cv2.COLOR_BGR2RGB)  # Convert the image to RGB (OpenCV uses BGR)
                    pil_im3 = Image.fromarray(embProjIm)  # Pass the image to PIL
                    draw = ImageDraw.Draw(pil_im3)  #
                    fontPil3 = ImageFont.truetype("arial.ttf", size=18, index=0)  # use a truetype font
                    fontPil4 = ImageFont.truetype("ariali.ttf", size=18, index=0)  # use a truetype font
                    embName = emb.folder.split('/')[4]
                    if emb.folder.split('/')[3] == 'GLS':
                        if emb.folder.split('/')[
                            4] == repEmbGLS:  # this if statement is in place to determine which image
                            # is the representative embryo
                            draw.text((5, 5), '**' + embName + '**', font=fontPil4,
                                      fill=(255, 255, 255))  # Draw the text
                        else:
                            draw.text((5, 5), embName, font=fontPil3, fill=(255, 255, 255))  # Draw the text
                            #                     pil_im3.show()
                    else:
                        #                         if emb.folder.split('/')[4]=='Emb1':
                        #                             draw.text((len(embProj)+40,10), str((t)*20)+'m', font=fontPil3, fill=(255,255,255)) # Note- time displayed is GLS normalized time
                        if emb.folder.split('/')[4] == repEmbMS:
                            draw.text((5, 5), '**' + embName + '**', font=fontPil4,
                                      fill=(255, 255, 255))  # Draw the text
                        else:
                            draw.text((5, 5), embName, font=fontPil3, fill=(255, 255, 255))  # Draw the text

                    embProj = cv2.cvtColor(np.array(pil_im3), cv2.COLOR_RGB2BGR)

                    embProjTime.append(embProj)  # adds to list for each emb
                maxProjEmbs.append(
                    embProjTime)  # appends each list for each embryo to maxProjEmbs for each RNAi for each strain
                embProjTime = np.array(embProjTime)  # converts to array (this is mostly for next line)
                shape.append(
                    embProjTime.shape)  # appends shape of each emb (this is needed for layout of all embs in final output)
                emb.image = None  # deletes images to save on memory
            shape = np.array(shape)
            time = max(shape[:, 0])  # finds number of timepoints for the embryo with the longest timeseries
            width = max(shape[:, 2])  # finds width of the widest image
            height = sum(shape[:, 1]) + (len(
                maxProjEmbs) - 1) * space  # finds the sum of all heights for each RNAi and adds in extra space for in between images
            offset = 0  # len(maxProjEmbs)+2 #not sure what is reasonable offset
            image = np.zeros((time, height + offset, int(width + 0.5 * offset),
                              3))  # creates a black background to place max projects

            yTot = 0  # used to define spacing as each embryo is added
            for i in range(len(maxProjEmbs)):  # compiles all timepoints and stacks movies for each strain
                maxP = maxProjEmbs[i]
                t = shape[i][0]
                x = shape[i][2]  # defines length of projection
                y = shape[i][1]  # defines height of projection
                #                 eName= name[i]#i.e. emb1
                image[:, yTot:yTot + y, int(0.25 * offset):x,
                :] = maxP  # redefines pixel values for specified coordinates [yRange,xRange] on black background to maxP image values
                #                 for timept in range(len(maxP)): #this loop is for writing embryo number labels on all images
                #
                #                     if emb.folder.split('/')[3]=='GLS': #this if statement is in place to determine which image is the representative embryo
                #                         pass
                # #                         if eName==repEmbGLS:
                # #                             font= cv2.FONT_HERSHEY_DUPLEX
                # #                             cv2.putText(image[timept],eName,(10,yTot+20), font, 0.6,(100,255,255),1)#cv2.putText(img, text, (x,y), fontFace, fontScale, color[, thickness[, lineType[, bottomLeftOrigin]]])
                # #                         else:
                # #                             font= cv2.FONT_HERSHEY_DUPLEX
                # #                             cv2.putText(image[timept],eName,(10,yTot+20), font, 0.6,(255,255,255),1)#cv2.putText(img, text, (x,y), fontFace, fontScale, color[, thickness[, lineType[, bottomLeftOrigin]]])
                #                     else:
                #                         if eName==repEmbMS:
                #                             font= cv2.FONT_HERSHEY_DUPLEX
                #                             cv2.putText(image[timept],eName,(10,yTot+20), font, 0.6,(100,255,255),1)#cv2.putText(img, text, (x,y), fontFace, fontScale, color[, thickness[, lineType[, bottomLeftOrigin]]])
                #                         else:
                #                             font= cv2.FONT_HERSHEY_DUPLEX
                #                             cv2.putText(image[timept],eName,(10,yTot+20), font, 0.6,(255,255,255),1)#cv2.putText(img, text, (x,y), fontFace, fontScale, color[, thickness[, lineType[, bottomLeftOrigin]]])

                yTot += y + space  # adds y to yTot so that next image is properly positioned below previous image

            image = np.array(image).astype(np.uint8)  # required because black background is not uint8
            for tp in image[:1]:
                fig = myFigure()
                fig.imshow(tp, colorbar=False)  # plots image onto figure
                fig.noAxis()
            # fig.show()

            #             saveImagesMulti(image, 'C:\\Users\\Becky\\Documents\\ODLAB\\EMBD_project\\Automated_Analysis\\{0}_{1}.tif'.format('EMBD{0:04}'.format(self.RNAi), emb.folder.split('/')[3]))

            if emb.folder.split('/')[-3] == 'GLS':
                imageGLS = image
            else:
                imageMS = image

        widthComp = max(imageGLS.shape[2], imageMS.shape[2])  # finds the width of the widest set of images (GLS or MS)
        heightComp = max(imageGLS.shape[1],
                         imageMS.shape[1])  # finds the length of the widest set of images (GLS or MS)

        hOffset = 50
        imageComp = np.zeros((time, heightComp + hOffset, int(2 * widthComp + 0.5 * offset),
                              3))  # creates a black background to place max projects
        imageComp = np.array(imageComp).astype(np.uint8)  # makes black background 8 bit
        #         for tmpt in range(len(imageComp)):#adds embd number and c.e.ID label to the top of the movie for each timepoint
        #             font2= cv2.FONT_HERSHEY_TRIPLEX+cv2.FONT_ITALIC
        #             font3= cv2.FONT_HERSHEY_TRIPLEX
        #             cv2.putText(imageComp[tmpt],'EMBD{0:04}   {1}'.format(self.RNAi,getGeneName(self.RNAi), getGeneCommonName(self.RNAi)),(15,20), font3, 0.6,(255,255,255),1)#adds the EMBD#, CE ID, and common gene name to the header of the compiled movie file
        #             cv2.putText(imageComp[tmpt],'{0}'.format(getGeneCommonName(self.RNAi)),(widthComp,20), font2, 0.6,(255,255,255),1)#adds the EMBD#, CE ID, and common gene name to the header of the compiled movie file
        #                 cv2.putText(imageComp[tmpt],getGeneName(self.RNAi),(widthComp,20), font2, 0.7,(255,255,255),1)

        text_to_show = 'EMBD{0:04}'.format(self.RNAi)
        if self.RNAi == 1903:
            text_to_show2 = 'hmg-3/hmg-4'
        elif getGeneName(self.RNAi)== getGeneCommonName(self.RNAi):
            text_to_show2 = '{0}'.format(getGeneCommonName(self.RNAi))
        else:
            text_to_show2 = '{1}({0})'.format(getGeneName(self.RNAi), getGeneCommonName(self.RNAi))
        cv2_im_rgb = cv2.cvtColor(imageComp[0], cv2.COLOR_BGR2RGB)  # Convert the image to RGB (OpenCV uses BGR)
        #             pil_im= Image.frombytes("L",imageComp.size, imageComp.tobytes())
        pil_im = Image.fromarray(cv2_im_rgb)  # Pass the image to PIL
        draw = ImageDraw.Draw(pil_im)  #
        fontPil = ImageFont.truetype("arial.ttf", size=22, index=0)  # use a truetype font
        fontPil2 = ImageFont.truetype("ariali.ttf", size=22,
                                      index=0)  # use a truetype font (arial italic is 'ariali' and arial bold is 'arialbd')
        draw.text((15, 10), text_to_show, font=fontPil, fill=(255, 255, 255))  # Draw the text
        draw.text((widthComp - 80, 10), text_to_show2, font=fontPil2, fill=(255, 255, 255))  # Draw the text
        #         pil_im.show()
        #         cv2_im_processed = cv2.cvtColor(np.array(pil_im), cv2.COLOR_RGB2BGR)#Get back the image to OpenCV

        cv2_im_processedAll = []
        for tmpt in range(-2, 27):
            im_pil_t = copy(pil_im)  # makes a copy of imageComp, so old timestamp is not there for the next timept
            draw = ImageDraw.Draw(im_pil_t)
            #             draw.text((2*widthComp-60,10), str((tmpt-1)*20)+'m', font=fontPil, fill=(0,0,0)) # Note- time displayed is GLS normalized time
            draw.text((2 * widthComp - 60, 10), str((tmpt) * 20) + 'm', font=fontPil,
                      fill=(255, 255, 255))  # Note- time displayed is GLS normalized time
            cv2_im_processed = cv2.cvtColor(np.array(im_pil_t), cv2.COLOR_RGB2BGR)  # Get back the image to OpenCV
            cv2_im_processedAll.append(cv2_im_processed)

        cv2_im_processedAll = np.array(cv2_im_processedAll)
        imageComp[:, :cv2_im_processedAll.shape[1], :cv2_im_processedAll.shape[2], :] = cv2_im_processedAll
        imageComp[:, hOffset:imageGLS.shape[1] + hOffset, :imageGLS.shape[2],
        :] = imageGLS  # re-assigns pixel values for GLS composite on black background to left edge
        imageComp[:, hOffset:imageMS.shape[1] + hOffset, widthComp:widthComp + imageMS.shape[2],
        :] = imageMS  # re-assigns pixel values for MS images to the right of GLS images on black background
        for tp in imageComp[:1]:
            fig = myFigure()
            fig.imshow(tp, colorbar=False)  # plots image onto figure
            fig.noAxis()
        # fig.show()

        saveImagesMulti(imageComp, FOLDER_TIME_ALIGN_MOVIES + '{0}_Comp.tif'.format('EMBD{0:04}'.format(self.RNAi)))
        saveFlag = True  # saves to web content folder
        if saveFlag:
            folder = FILENAME_WEB_CONTENT + 'EMBD{n:04}/'.format(n=self.RNAi)
            if not os.path.exists(folder):
                os.makedirs(folder)
            saveImagesMulti(imageComp, folder + '{0}_Comp.tif'.format('EMBD{0:04}'.format(self.RNAi)))

    def get_human_gene_from_ortholist_embd(self, ce_id):
        from fileLookup import FILENAME_EMBD_ORTHOLIST
        import pandas as pd

        file = FILENAME_EMBD_ORTHOLIST
        data = pd.read_csv(file, ',')
        self_data = data[data['Locus ID'] == ce_id]
        self_data_sort = self_data.sort_values('No. of Databases', ascending=False)
        h_gene_list = self_data_sort['HGNC Symbol'].values.tolist()
        h_gene_str = ", ".join([str(elem) for elem in h_gene_list])
        print(h_gene_str)
        return (h_gene_str)

    def fix_head_params_RNAi(self):
        '''This function was created to specifically update two headInt params (devHead and tailHead) for a given
        RNAi condition.'''

        self.paramsUseMS = []
        self.getEmbryos()
        self.prev_GLS, self.prev_MS = self.getPrev()
        self.asignUse()
        # self.calcAverageParams()

        PARAMS_TO_FIX = ['tailHead','devHead']
        for pN in PARAMS_TO_FIX:
            goodParamsMS = []  # empty list to add good MS params to
            for emb in self.MSUseEmbs:
                if emb.checkGoodParam(pN): goodParamsMS.append(emb.params[pN])  # Appends value for each good parameter
            if len(goodParamsMS) > 0:  # checks to be sure there are values in the list
                self.paramsUseMS.append(pN)  # adds parameter names to paramsUseMS list
                self.paramsMS[pN] = np.median(goodParamsMS)  # adds median good parameters to dictionary self.params
                # self.paramsMSstd[pN] = np.std(goodParamsMS)  # adds std for RNAi to dictionary
                # self.paramsEmbsMS[pN] = goodParamsMS

        # self.paramsMS = make_dict_float
        # self.save_params_data()

        conn, curs = db_utils_embryos.initialize_db()
        columns = {}
        update_set = ['tailHead','devHead']
        for p in update_set:
        # for key in self.paramsMS:
            columns[p + '_MS'] = self.paramsMS[p]
        # columns['prev_MS'] = self.prev_MS
        db_utils_embryos.update_row(self.id, columns, curs, 'genes')
        conn.commit()
        conn.close()




        # conn, curs = db_utils_embryos.initialize_db()
        # columns = {}
        # for key in self.paramsGLS:
        #     columns[key + '_GLS'] = self.paramsGLS[key]
        # for key in self.paramsMS:
        #     columns[key + '_MS'] = self.paramsMS[key]
        # columns['prev_GLS'] = self.prev_GLS
        # columns['prev_MS'] = self.prev_MS
        # columns['version'] = self.__version__
        # db_utils_embryos.update_row(self.id, columns, curs, 'genes')
        # conn.commit()
        # conn.close()


class ControlClass(RNAiClass):
    def __init__(self, verb=True):
        RNAiClass.__init__(self, 0, verb)

    def asignUse(self):
        self.GLSUseEmbs = self.GLSMoveEmbs
        self.MSUseEmbs = self.MSMoveEmbs

    def defineFolders(self):
        folder = FOLDER_IN + 'cropped/EMBD0000/GLS/'
        badFolders = ['20150129T130903/', '20140813T133845/']
        self.GLSfolder = [folder + d + '/' for d in os.listdir(folder) if d[0] == '2' and d not in badFolders][:20]
        folder = FOLDER_IN + 'cropped/EMBD0000/MS/'
        badFolders = ['20140514T130911/', '20140515T131046/']
        self.MSfolder = [folder + d + '/' for d in os.listdir(folder) if d[0] == '2' and d not in badFolders][:20]

    def getEmbryos(self): #note- need to define folders firstv
        embs = []


        for f in self.GLSfolder:
            try:
                embs += loadEmbs(f)  # , inds = range(1))
            except:
                pass

        # for emb in embs:
            # emb.updateParams()
        # embs = embdFunc.updateEmbsMulti(embs)
        #         embT = []
        #         for emb in embs:
        #             if emb.spotsLoaded: embT.append(emb)
        self.GLSMoveEmbs, self.GLSNoMoveEmbs = self.splitByMovement(
            embs)  # generates two lists- moving and not moving GLS embs
        embs = []

        for f in self.MSfolder:
            try:
                embs += loadEmbs(f)  # , inds = range(1))
            except:
                pass

        # for emb in embs:
            # emb.updateParams()
        # embs = embdFunc.updateEmbsMulti(embs)
        self.MSMoveEmbs, self.MSNoMoveEmbs = self.splitByMovement(
            embs)  # generates two lists- moving and not moving MS embs


def get_control_origin():
    con, cur = db_utils_embryos.initialize_db()
    c_mean = {}
    c_std = {}
    for pn in PARAM_NAMES:
        print(pn)
        pars = db_utils_embryos.get_all_params(cur, pn, 'GLS', 0)
        if not np.isnan(pars).all():
            c_mean[pn + '_GLS'] = np.nanmean(pars)
            c_std[pn + '_GLS'] = np.nanstd(pars)
        else:
            c_mean[pn + '_GLS'] = np.nan
            c_std[pn + '_GLS'] = np.nan
        pars = db_utils_embryos.get_all_params(cur, pn, 'MS', 0)
        if not np.isnan(pars).all():
            c_mean[pn + '_MS'] = np.nanmean(pars)
            c_std[pn + '_MS'] = np.nanstd(pars)
        else:
            c_mean[pn + '_MS'] = np.nan
            c_std[pn + '_MS'] = np.nan
    return c_mean, c_std


def write_controlAvg_mysql():
    """ pulls average parameter values from control embryos (uses origin dictionary) and writes average values for each
    parameter to genes mySQL table for embd_id 0. Note this is not in db_utils_genes because there is an circular import
    issue with ControlClass when this function is placed in that file."""

    conn, curs = db_utils_embryos.initialize_db()
    table = "genes"
    id = 2627  # this is the unique identifier row number id for EMBD0000
    columns={}  # creates dictionary of column labels and values for update_row function to use
    c = ControlClass()
    control_ave_dict = c.origin   #origin is a dictionary of ave values for parameters from control embryos
    for key in control_ave_dict:
        columns[key] = control_ave_dict[key]  # reads param name and value from c.origin and assigns column name and
        # value to dictionary to populate mySQL
        # print key, control_ave_dict[key]
    db_utils_embryos.update_row(id, columns, curs, table)  # updates row containing control values
    conn.commit()
    conn.close()




if __name__ == '__main__':
    '''run a list of RNAs'''
    skipped = ['EMBD1537', 'EMBD1538', 'EMBD1539', 'EMBD1540', 'EMBD1541', 'EMBD1543', 'EMBD1544', 'EMBD1546', 'EMBD1547', 'EMBD1552', 'EMBD1553', 'EMBD1560', 'EMBD1567', 'EMBD1570', 'EMBD1572', 'EMBD1573', 'EMBD1576', 'EMBD1577', 'EMBD1578', 'EMBD1580', 'EMBD1581', 'EMBD1586', 'EMBD1588', 'EMBD1590', 'EMBD1592', 'EMBD1594', 'EMBD1595', 'EMBD1596', 'EMBD1599', 'EMBD1603', 'EMBD1604', 'EMBD1605', 'EMBD1606', 'EMBD1607', 'EMBD1609', 'EMBD1610', 'EMBD1611', 'EMBD1616', 'EMBD1618', 'EMBD1619', 'EMBD1621', 'EMBD1622', 'EMBD1623', 'EMBD1625', 'EMBD1626', 'EMBD1627', 'EMBD1630', 'EMBD1631', 'EMBD1632', 'EMBD1633', 'EMBD1634', 'EMBD1635', 'EMBD1636', 'EMBD1639', 'EMBD1640', 'EMBD1641', 'EMBD1642', 'EMBD1644', 'EMBD1645', 'EMBD1646', 'EMBD1647', 'EMBD1648', 'EMBD1649', 'EMBD1650', 'EMBD1652', 'EMBD1653', 'EMBD1654', 'EMBD1656', 'EMBD1658', 'EMBD1661', 'EMBD1662', 'EMBD1663', 'EMBD1664', 'EMBD1665', 'EMBD1666', 'EMBD1667', 'EMBD1670', 'EMBD1671', 'EMBD1672', 'EMBD1673', 'EMBD1674', 'EMBD1675', 'EMBD1677', 'EMBD1678', 'EMBD1679', 'EMBD1680', 'EMBD1683', 'EMBD1684', 'EMBD1685', 'EMBD1687', 'EMBD1688', 'EMBD1689', 'EMBD1690', 'EMBD1694', 'EMBD1696', 'EMBD1697', 'EMBD1699', 'EMBD1701', 'EMBD1702', 'EMBD1703', 'EMBD1706', 'EMBD1707', 'EMBD1708', 'EMBD1710', 'EMBD1711', 'EMBD1712', 'EMBD1713', 'EMBD1714', 'EMBD1716', 'EMBD1718', 'EMBD1719', 'EMBD1720', 'EMBD1721', 'EMBD1722', 'EMBD1724', 'EMBD1725', 'EMBD1726', 'EMBD1727', 'EMBD1728']
    #
    # runlist = [f for f in range(1766,1767)]
    # from fileLookup import redo_batch2
    # runlist = redo_batch2

    #
    runlist = [1906,1907,1908,1909]
    # for i in runlist:
    #     print(i)
    #     if 'EMBD{0}'.format(i) in skipped:
    #         print("skipped {0}".format(i))
    #         pass
    #     else:
    #         print('generating', i)
    #         r = RNAiClass(i)
    #         r.refresh_embs()
    #         r.refresh_params_data()



    '''get PADs for RNAi conditions outside the first 500'''
    # query = 3043
    # r = RNAiClass(query)
    # # r2_list = [i for i in range(504)]+[1903, 1904, 1905, 1360,703,885,1406,674,952,866]
    # r2_list = [i for i in range(504)]
    #
    # all_pads = []
    # for r2 in r2_list:
    #     test = RNAiClass(r2)
    #     print([query, r2, r.getPAD(test)])
    #     all_pads.append([query, r2, r.getPAD(test)])
    # import pandas as pd
    # df = pd.DataFrame(all_pads)
    # df_sorted = df.sort_values([2], ascending=False)
    #
    # # print(df)
    # print(df_sorted.head(n=40))

    # query = 3043
    # r = RNAiClass(query)
    # r2 = 1600
    # test = RNAiClass(r2)
    # print([query, r2, r.getPAD(test)])

    '''get representative embryo for query'''
    # r = RNAiClass(1527)
    # r2 = RNAiClass(57)
    # print(r.getPAD(r2))
    # print(r.getPAD(r2))
    # GLS, MS = r.getRepEmb(r2)

    '''make time aligned movies'''
    from fileLookup import SKIPPED
    # #
    # run_list = [i for i in [3045]]
    run_list = [1906,1907,1908,1909]

    for j in run_list:
        if 'EMBD{0:04}'.format(j) not in SKIPPED:
            print("STARTING {0} ---------------------").format(j)
            r = RNAiClass(j)
            r.makeTalignMovies()

    # r = RNAiClass(3020)

    # print(r.getDist2Zero())

    # r1 = RNAiClass(76)
    # r.getEmbryos()
    # r.asignUse()

    # r1.getEmbryos()
    # r1.asignUse()
    # r.setNDPosition()
    # r1.setNDPosition()
    #
    # cos_dist = r.get_cosine_dist(r1)
    # print("cosine distance = {0}".format(cos_dist))

    # print(r.getDist2Zero())
    # gls_spots, ms_tInt = r.plotSigmFitAll()
    # r.plotSigmFitAvg()

    # r.plotMOIAvg()

    # gls_spots.show()
    # ms_tInt.show()



    # print(len(r.pNamesUse[0]) + len(r.pNamesUse[1]))
    # print(len(r1.pNamesUse[0]) + len(r1.pNamesUse[1]))
    # print("440 to 52 dist is {0}".format(r.getDistance(r1)))
    # print("52 to 440 dist is {0}".format(r1.getDistance(r)))

    # c = ControlClass()
    # r = RNAiClass(5)
    # print c.paramsGLS
    # print c.paramsGLS['maxG']
    # print np.mean(c.paramsGLS['maxG'])
    # print np.std(c.paramsGLS['maxG'])
    # print get_control_origin()

    # write_controlAvg_mysql()

    # r = RNAiClass(67)
    # # r.refresh_embs()
    # r.refresh_params_data()

    # r_list = [91,180,205,205,256,256,281,284,336,336,336,336,336,441,441,441,441,441,441,441,444]
    # r_list = np.unique(r_list)
    # for i in r_list:
    #     print('generating', i)
    #     r = RNAiClass(i)
    #     # r.refresh_embs()
    #     r.refresh_params_data()

    # rList = []
    # tscale_gls = []
    # tscale_ms = []
    # tscale_avg = []
    # for r in range(504):
    #     rnai = RNAiClass(r)
    #     if rnai.prev == 0:
    #         rList.append(r)
    #         tscale_gls.append(rnai.paramsGLS['tScale'])
    #         tscale_ms.append(rnai.paramsMS['tScale'])
    #         tscale_avg.append(np.mean((rnai.paramsGLS['tScale'], rnai.paramsMS['tScale'])))
    # ind_min_gls = np.nanargmin(tscale_gls)
    # ind_min_ms = np.nanargmin(tscale_ms)
    # ind_min_avg = np.nanargmin(tscale_avg)
    # # ind_min_gls = np.nanargmax(tscale_gls)
    # # ind_min_ms = np.nanargmax(tscale_ms)
    # # ind_min_avg = np.nanargmax(tscale_avg)
    # print('GLS: {gls} tScale = {t_gls}; MS: {ms} tScale = {t_ms}; Both: {avg} tScale = {t_avg};'.format(gls=rList[ind_min_gls],
    #                                                                                      t_gls=tscale_gls[ind_min_gls],
    #                                                                                      ms=rList[ind_min_ms],
    #                                                                                      t_ms=tscale_ms[ind_min_ms],
    #                                                                                      avg=rList[ind_min_avg],
    #                                                                                      t_avg=tscale_gls[ind_min_avg]))

    #     r1.refresh_params_data()
    #     print(r1.paramsGLS['aG'])
    #     print(r1.paramsGLS['aR'])
    #     print(r1.paramsGLS['aY'])
    #     r1.refresh_params_data()
    #     r1 = RNAiClass(441)
    #     r2.refresh_params_data()
    #     print(r2.paramsGLS['aG'])
    #     print(r2.paramsGLS['aR'])
    #     print(r2.paramsGLS['aY'])
    #     print(r1.getPAD(r2))
    #     print(r1.getPAD(r2))
    #     embGLS, embMS = r1.getRepEmb(r2)
    #     r1.makeTalignMovies()


    '''make time aligned movies'''
    # for i in range(55,56):
    #     print("starting {0}".format(i))
    #     r=RNAiClass(i)
    #     r.makeTalignMovies()

    '''get PAD or representative embryo'''
    #     r1 = RNAiClass(32)
    #     r2 = RNAiClass(441)
    #     print(r1.getPAD(r2))
    #     embGLS, embMS = r1.getRepEmb(r2)
    #     print('representative embryos: {lg} GLS, {lm} MS'.format(lg = embGLS.label, lm = embMS.label))
    #     print('embryos used in GLS:',[emb.label for emb in r1.GLSUseEmbs])
    #     print('embryos used in MS:',[emb.label for emb in r1.MSUseEmbs])

    '''test'''
# image = np.zeros(252,428,3)#creates a black background to place max projects
#     image = 'C:\Users\Becky\Pictures\dylan\IMAG0160.jpg'
# #     image= np.zeros(255,255)
#     img= cv2.imread(image,cv2.IMREAD_GRAYSCALE)
#     cv2.line(img,(0,0), (150,50), (255, 255, 255), 15)
# 
#     cv2.imshow('test',img)
#     cv2.waitKey(0)
#     cv2.destroyAllWindows()

    '''get human gene from ortholist'''
    # r = RNAiClass(1)
    # h_genes = r.get_human_gene_from_ortholist_embd('K11C4.3')

    '''fix head params'''

    # for i in range (254,255): # comment out lines 202-204 in EMBRYOS to prevent embryos from being refreshed while this runs. DONT forget to put it back!
    #     print(i)
    #     r = RNAiClass(i)
    #     r.fix_head_params_RNAi()



    '''other'''
    # r = RNAiClass(364) # NOTE THAT CONTROL ORIGIN HAS BEEN CHANGED FOR tailHead and devHead as of 04/12/22!
    # r2 = RNAiClass(1903)
    # PAD = r.getPAD(r2)
    # dims = r.getDims()
    # dist = r.getDist2Zero()
    # dist_RNAi = r.getDistance(r2)
    # print((dist_RNAi, dist, dims, PAD))


    # r = RNAiClass(1904)  # NOTE THAT CONTROL ORIGIN HAS BEEN CHANGED FOR tailHead and devHead as of 04/12/22!
    # list = [364,703,885,408,130,1360,674,1406,866, 952, 1903, 1904]
    #
    #
    # for gene in list:
    #     r2 = RNAiClass(gene)
    #     PAD = r.getPAD(r2)
    #     print(gene, PAD)



