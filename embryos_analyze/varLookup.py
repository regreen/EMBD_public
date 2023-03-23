'''
Created on May 2, 2016

@author: Becky

Contains lookup values and files for other programs
'''
global PARAM_NAMES
import time
import json

nCores = 6  # number of parallel processes to run
FOLDER_IN = 'Z:/'
# MYSQL_DB_CREDENTIALS = credential file  # Analysis station

PARAM_NAMES = ['maxR', 'maxG', 'aG', 'aGErr', 'bG', 'bGErr', 'mG', 'mGErr', 'sG', 'sGErr', 'rG', 'rGErr', 'Gtail',
               'GtailErr', 'residsEarlyG', 'residsLateG', \
               'aR', 'aRErr', 'bR', 'bRErr', 'mR', 'mRErr', 'sR', 'sRErr', 'rR', 'rRErr', 'Rtail', 'RtailErr',
               'residsEarlyR', 'residsLateR', \
               'aY', 'aYErr', 'bY', 'bYErr', 'mY', 'mYErr', 'sY', 'sYErr', 'rY', 'rYErr', 'Ytail', 'YtailErr',
               'residsEarlyY', 'residsLateY', \
               'tMove', 'movement', 'mRmG', 'mYmG', 'mYmR', 'tScale'] + \
              ['fracR', 'fracG', 'fracY'] + \
              ['aSigHead', 'mSigHead', 'sSigHead', 'aGaussHead', 'mGaussHead', 'sGaussHead', 'maxHead', 'tailHead',
               'scaleHead', 'devHead'] + \
              ['scaleLength', 'devLength', 'tailLength'] + \
              ['CoM{0}{1}{2}_GLS'.format(ax, c, par) for c in ['R', 'Y'] for ax in range(2) for par in
               ['td', 'scale', 'avgRes', 'stdTail']] + \
              ['MoI{0}{1}{2}_GLS'.format(ax, c, par) for c in ['G', 'R', 'Y'] for ax in range(2) for par in
               ['td', 'scale', 'avgRes', 'stdTail']] + \
              ['CoM{0}{1}{2}_MS'.format(ax, c, par) for c in ['G', 'R'] for ax in range(2) for par in
               ['td', 'scale', 'avgRes', 'stdTail']] + \
              ['MoI{0}{1}{2}_MS'.format(ax, c, par) for c in ['G', 'R'] for ax in range(2) for par in
               ['td', 'scale', 'avgRes', 'stdTail']]

# movement and alignment validation file
FILENAME_ALIGNMENT_0000_MS = FOLDER_IN + 'EMBD_Alignment_Validation_EMBD0000_MS.csv'
FILENAME_ALIGNMENT_RNAI_MS = FOLDER_IN + 'EMBD_Alignment_Validation_RNAi_MS.csv'
FILENAME_ALIGNMENT_0000_GLS = FOLDER_IN + 'EMBD_Alignment_Validation_EMBD0000_GLS.csv'
FILENAME_ALIGNMENT_RNAI_GLS = FOLDER_IN + 'EMBD_Alignment_Validation_RNAi_GLS.csv'

# obtained from avgCalculator.py
FILENAME_AVG_CURVES_MS_ADJTH = FOLDER_IN + 'AVG_CURVES_0000_MS.pkl'  # (adjustable intensity threshold) pickle of the average curve for the distance between the green and red signal in control embryos of MS

FILENAME_AVG_CURVES_GLS = FOLDER_IN + 'AVG_CURVES_0000_GLS.pkl'  # pickle of the average curve for the distance between the green and red signal in control embryos of MS

AVG_CURVES_MS = ['tIntR', 'tIntG', 'CoM0R', 'CoM1R', 'CoM0G', 'CoM1G', 'MoI0G', 'MoI1G', 'MoI0R', 'MoI1R',
                 'lengthR', 'headInt']
AVG_CURVES_GLS = ['spotsG', 'spotsR', 'spotsY', 'CoM0R', 'CoM1R', 'CoM0Y', 'CoM1Y', \
                  'MoI0G', 'MoI1G', 'MoI0R', 'MoI1R', 'MoI0Y', 'MoI1Y']

# obtained from avgCalculator.py
RADIUS_AVG = 53.
LENGTH_AVG = 100.

AG_AVG_GLS = 90  # 90+/-7 #run022318
AR_AVG_GLS = 110  # 110+/-11
AY_AVG_GLS = 82  # 82+/-5
SG_AVG_GLS = 1.3  # 1.3+/-0.4
SR_AVG_GLS = 1.4  # 1.4+/-0.3
SY_AVG_GLS = 1.4  # 1.4+/-0.2
T0R_M_T0G_AVG_GLS = 3.0  # 3.0+/-0.8
T0Y_M_T0G_AVG_GLS = 3.7  # 3.7+/-0.7
T0Y_M_T0R_AVG_GLS = 0.7  # 0.7+/-0.6

GLS_MOVE_PRECISION = 4
# AG_AVG_GLS = 89  # old version of averages
# AR_AVG_GLS = 111
# AY_AVG_GLS = 83
# SG_AVG_GLS = 1.3
# SR_AVG_GLS = 1.5
# SY_AVG_GLS = 1.3
# T0R_M_T0G_AVG_GLS = 3.0
# T0Y_M_T0G_AVG_GLS = 3.7
# T0Y_M_T0R_AVG_GLS = 0.7
T0_2_MOVE_AVG_GLS = 17

AR_AVG_MS = 6.1e8  # 614736852+/-110817335 #run022318
AG_AVG_MS = 4.4e8  # 440091691+/-108409019
SR_AVG_MS = 2.6  # 2.6+/-0.5
SG_AVG_MS = 2.2  # 2.2+/-0.5
AVG_R_BASE_MS = 8.7e7  # 87318240+/-110817335
AVG_G_BASE_MS = 2.1e8  # 212369928+/-108409019
T0G_M_T0R_AVG_MS = 3  # 3+/-1
T0_2_MOVE_AVG_MS = 12  # 12+/-2

MS_MOVE_PRECISION = 4
# AR_AVG_MS = 5.9e8  # old version of averages
# AG_AVG_MS = 4.2e8
# SG_AVG_MS = 2.2
# SR_AVG_MS = 2.7
# AVG_G_BASE_MS = 2.0e8
# AVG_R_BASE_MS = 1.1e8
# T0G_M_T0R_AVG_MS = 3.
# T0_2_MOVE_AVG_MS = 12
MG_2_MOVE_AVG_MS = 8 #use to check if movement that is found makes sense
MG_2_MOVE_STD_MS = 2

DT0_GLS_2_MS = T0_2_MOVE_AVG_GLS - T0_2_MOVE_AVG_MS
# LOG_FILE = FOLDER_IN + 'LOG/logFile{0}.txt'.format(time.strftime('%m_%d_%Y'))

# important variables used in multiple places
nLambda_MS = 6.  # used in fitting cutoff in curves- determined empirically, see googleDoc lambdaOptimization_MS_022218
tauScale_MS = 0.0  # used in fitting time scaling in curves
nLambda_GLS = 3.  # used in fitting cutoff in curves- determined empirically, see googleDoc lambdaOptimGLS_022318
tauScale_GLS = 0.  # used in fitting time scaling in curves


def get_log_file():
    return FOLDER_IN + 'LOG/logFile{0}.txt'.format(time.strftime('%m_%d_%Y'))


def printLog(string):
    print(string)
    try:
        with open(get_log_file(), "a") as myfile:
            myfile.write(time.strftime("%a, %d %b %Y %H:%M:%S || ", time.localtime()) + json.dumps(string) + '\n')
    except Exception, e:
        print('Can not print to log file Error:{0}'.format(str(e)))
