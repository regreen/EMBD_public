'''
Created on Mar 7, 2016

@author: Becky
'''
import numpy as np
from varLookup import *
from params_use_assign import get_params_use
from RNAiClass import RNAiClass, ControlClass
import PADassist
from itertools import repeat
import multiprocessing as mp
from lmfit import Parameters, minimize, report_fit, Model, Parameter
import embdFunc

PARAM_NAMES_USE_GLS, PARAM_NAMES_USE_MS = get_params_use()


def getCSIS(input):
    distAll, i = input
    print('csi', i)
    csiRNAi = []
    for j in range(len(distAll[i])):
        n1 = len(distAll[
                     i]) - j  # number of genes with distance less than interacting gene j including this gene (count below, because small distance ~ high pcc)
        r1, r2, d, d20, nOv, pad = distAll[i][j]
        rcounted = [rt2 for rt1, rt2, dt, d20, nOv, pad in distAll[i][j:]]
        n2 = getNBelow(distAll, r2, d, rcounted)
        csiRNAi.append(1 - 1. * (n1 + n2) / len(distAll))
    return csiRNAi


def getCSI(distAll):
    pool = mp.Pool(processes=6)
    csiListAll = pool.map(getCSIS, zip(repeat(distAll), range(len(distAll))))
    pool.close()
    pool.join()
    return csiListAll


def getCSIPAD(distAll):
    csiListAll = []
    for i in range(len(distAll)):
        print('csi', i)
        csiRNAi = []
        for j in range(len(distAll[i])):
            padth = distAll[i][j][5]
            data = np.array(distAll[i])
            pads = data[:, 5].astype(np.float)
            n1 = np.sum(pads >= padth)
            rcounted = data[:, 1][np.where(pads >= padth)]
            n2 = getNAbovePAD(distAll, r2, d, rcounted)
            csiRNAi.append(1 - 1. * (n1 + n2) / len(distAll))
        csiListAll.append(csiRNAi)
    return csiListAll


def getNBelow(distAll, rna, dist, rcounted):
    for i in range(len(distAll)):
        if distAll[i][0][0] == rna:
            dists = [d for (r1, r2, d, d20, nOv, pad) in distAll[i] if d <= dist and not (r2 in rcounted)]
            return len(dists)
    return 0


def getNAbovePAD(distAll, rna, padth, rcounted):
    data = np.array(distAll)
    data = data[np.where(data[:, 0] == rna)[0]]
    np.delete(data, )
    for i in range(len(distAll)):
        if distAll[i][0][0] == rna:
            pads = [pad for (r1, r2, d, d20, nOv, pad) in distAll[i] if pad >= padth and not (r2 in rcounted)]
            return len(pads)
    return 0


def getPADS(input):
    RNAi, i = input
    print('distance', i)
    distRNAi = []
    for j in range(len(RNAi)):
        #         distRNAi.append((RNAi[i].label, RNAi[j].label, 1-RNAi[i].getPAD(RNAi[j]),\
        distRNAi.append((RNAi[i].label, RNAi[j].label, RNAi[i].getDistance(RNAi[j]), \
                         RNAi[i].getDist2Zero(), RNAi[i].getOverlapDims(RNAi[j]), RNAi[i].getPAD(RNAi[j])))
    distRNAi.sort(key=lambda x: x[2])
    d = distRNAi[::-1]  # sort given rnai list by distance to another rnai (large first)
    return d


def optimizePAD(rGroups, c):
    pGLS = {}
    pMS = {}
    pStd = {}
    sumAll = 0.
    for key in c.paramsGLS:
        if key in PARAM_NAMES_USE_GLS:
            pGLS[key] = c.paramsGLS[key]
            pStd[key + '_GLS'] = c.paramsGLSstd[key]
            sumAll += c.paramsGLSstd[key]
    for key in c.paramsMS:
        if key in PARAM_NAMES_USE_MS:
            pMS[key] = c.paramsMS[key]
            pStd[key + '_MS'] = c.paramsMSstd[key]
            sumAll += c.paramsMSstd[key]

    nParams = len(pStd)
    eqWeights = {}
    for key in pStd:
        eqWeights[key] = 1. / nParams

    # Use all RNAi mean and std instead of controls
    pMeanRNAi, pStd = getParamsNorms()

    #     pGLS, pMS = {}, {}
    #     for key in pMeanRNAi:
    #         if key[-3:]=='_MS' and key[:-3] in PARAM_NAMES_USE_MS: pMS[key[:-3]] = pMeanRNAi[key]
    #         elif key[-4:]=='_GLS' and key[:-4] in PARAM_NAMES_USE_GLS: pGLS[key[:-4]] = pMeanRNAi[key]

    def getAvgPAD(weights):
        pads = []
        dist20 = []
        for rList in rGroups:
            for r in rList:
                r.setNDPosWeights(pGLS, pMS, pStd, weights)
                dist20.append(r.getDist2Zero())
            for i in range(len(rList) - 1):
                for j in range(i + 1, len(rList)):
                    r1 = rList[i]
                    r2 = rList[j]
                    pads.append(r1.getPAD(r2))
        if sum(np.isnan(pads)) == len(pads):
            print('bad value')
            return 10000.
        else:
            #             print (np.nanmean(pads),np.nanmean(dist20))
            return (1. - np.nanmean(pads))  # /np.nanmean(dist20)

    params = Parameters()
    expr = '1.'
    for key in PARAM_NAMES_USE_GLS[1:]:
        if not np.isnan(pGLS[key]):
            params.add(key + '_GLS', 1. / nParams, vary=True, min=0, max=1)
            expr += '-{0}_GLS'.format(key)
    for key in pMS:
        if not np.isnan(pMS[key]):
            params.add(key + '_MS', 1. / nParams, vary=True, min=0, max=1)
            expr += '-{0}_MS'.format(key)
    params[PARAM_NAMES_USE_GLS[0] + '_GLS'] = Parameter(value=1. / nParams, expr=expr, min=0., max=1.)
    minzer = minimize(getAvgPAD, params=params, method='cobyla', options={'maxiter': 1000})
    res = minzer.params
    print('successful fit ', minzer.success)
    print(report_fit(res))
    print('weight sum={0}'.format(np.sum([res['{0}'.format(key)].value for key in params])))
    print('number of parameters={0}'.format(nParams))
    print('optimized weights PAD={0}'.format(1 - getAvgPAD(res)))
    print('equal weights PAD={0}'.format(1 - getAvgPAD(eqWeights)))


def getParamsNorms():
    RNAi = []
    for i in range(1, 504):
        try:
            r = RNAiClass(i, verb=True)
            RNAi.append(r)
        except Exception, e:
            printLog('!!!!!!!!!ERROR OCCURED RNAi {0} :{1}'.format(i, str(e)))

    pAllRNAi, pMean, pSTD = {}, {}, {}
    for key in PARAM_NAMES:
        pAllRNAi[key + '_GLS'] = []
        pAllRNAi[key + '_MS'] = []
        for r in RNAi:
            #             print(key, r.paramsEmbsGLS[key], r.paramsEmbsMS[key])
            pAllRNAi[key + '_GLS'].append(r.paramsGLS[key])
            pAllRNAi[key + '_MS'].append(r.paramsMS[key])
        pMean[key + '_GLS'] = np.nanmean(pAllRNAi[key + '_GLS'])
        pSTD[key + '_GLS'] = np.nanstd(pAllRNAi[key + '_GLS'])
        pMean[key + '_MS'] = np.nanmean(pAllRNAi[key + '_MS'])
        pSTD[key + '_MS'] = np.nanstd(pAllRNAi[key + '_MS'])
    return pMean, pSTD


def get_RNAi_object_list():
    """
    :return: a list of RNAi objects for the genes in range
    """
    RNAi = []
    for i in range(1, 504):
    # for i in range(1, 6):
        try:
            r = RNAiClass(i, verb=True)
            #             if r.__version__!=embdFunc.get_rnai_revision_version(): r.refreshParamsData()
            r.setNDPosition()
            RNAi.append(r)
        except Exception, e:
            printLog('!!!!!!!!!ERROR OCCURED RNAi {0} :{1}'.format(i, str(e)))
    return RNAi


def calc_pairwise_metrics(RNAi):
    """
    calculates pairwise metrics (distance, PAD, CSI) for a list of RNAi conditions
    :param RNAi: list of RNAi objects
    :return: distAll: which is list of lists. each list includes RNAi[i].label, RNAi[j].label, RNAi[i].getDistance(RNAi[j]),\
                         RNAi[i].getDist2Zero(), RNAi[i].getOverlapDims(RNAi[j]), RNAi[i].getPAD(RNAi[j]))
            csi: CSI values for i, j pair

    """
    print("check params used!")
    pool = mp.Pool(processes=6)
    distAll = pool.map(getPADS, zip(repeat(RNAi), range(len(RNAi))))
    pool.close()
    pool.join()
    csi = getCSI(distAll)
    return distAll, csi


def write_phenodistinfo_to_CSV(RNAi, distAll, csi):
    """
    writes into csv file RNAi pair, distance, distance to wt, dimensions, CSI and PAD values for each pair of RNAi conditions
    :param RNAi: list of RNAi objects
    :param distAll: list of lists. each list includes RNAi[i].label, RNAi[j].label, RNAi[i].getDistance(RNAi[j]),\
                         RNAi[i].getDist2Zero(), RNAi[i].getOverlapDims(RNAi[j]), RNAi[i].getPAD(RNAi[j]))
    :param csi:
    :return: writes values to specified csv file
    """
    fd = open('Z:\\phenoDist_rank_{d}.csv'.format(d=time.strftime('%m%d%Y')), 'w')
    # fd = open('Z:\\phenoDist_rank_ALL_PARAMS{d}.csv'.format(d=time.strftime('%m%d%Y')), 'w')

    pUsed = ''
    for p in PARAM_NAMES_USE_GLS:
        pUsed += '{0},'.format(p)
    fd.write('params used GLS=,' + pUsed + '\n')
    pUsed = ''
    for p in PARAM_NAMES_USE_MS:
        pUsed += '{0},'.format(p)
    fd.write('params used MS=,' + pUsed + '\n')
    fd.write('RNA1,RNA2,dist,dist2WT,dims,CSI,PAD\n')
    for i in range(len(RNAi)):
        for j in range(len(RNAi)):
            r1, r2, d, d20, nOv, pad = distAll[i][j]
            #             print(r1,r2,csi[i][j],d)
            fd.write('{0},{1},{2},{3},{4},{5},{6}\n'.format(r1, r2, d, d20, nOv, csi[i][j], pad))
    fd.close()


def write_phenodistinfo_to_mySQL(RNAi, distAll, csi):
    """
    populates pairwise distance, pad and csi data into mySQL. Checks for existing entry and if it is present, it updates
     it, if not it adds it as a new entry
    :param RNAi: list of RNAi objects
    :param distAll:  list of lists. each list includes RNAi[i].label, RNAi[j].label, RNAi[i].getDistance(RNAi[j]),\
                         RNAi[i].getDist2Zero(), RNAi[i].getOverlapDims(RNAi[j]), RNAi[i].getPAD(RNAi[j]))
    :param csi:
    :return:
    """
    from db_utils_embryos import initialize_db, insert_row, update_row
    conn, cursor = initialize_db()
    table = "pair_lookup"

    for i in range(len(RNAi)):
        print(i)
        for j in range(len(RNAi)):  # get all possible unique combinations of i and j
            rnai1 = RNAi[i]
            rnai2 = RNAi[j]
            # sql_rows = cursor.fetchall()
            r1, r2, d, d20, ov_dims, pad = distAll[i][j]
            cs = csi[i][j]
            colms = {'rnai1': int(r1[4:]), 'rnai2': int(r2[4:]), 'dist': d, 'dist2wt': d20, 'dims': ov_dims, 'CSI': cs,
                     'PAD': pad}
            if r1 != r2:
                sql = "SELECT id FROM {t} WHERE rnai1={r1} AND rnai2={r2}".format(t=table, r1=int(r1[4:]), r2=int(r2[4:]))
                cursor.execute(sql)
                vals = cursor.fetchall()
                if len(vals)>0:
                    update_row(vals[0][0], colms, cursor, table)
                elif len(vals) == 0:
                    insert_row(colms, cursor, table)
        conn.commit()
    conn.close()

def run_phenodistcalc_populate():
    '''run phenoDistCalc to populate csv and/or mySQL'''
    RNAi = get_RNAi_object_list()  # retrieves a list of RNAi objects for the genes in range (range is defined within this function)
    distAll, csi = calc_pairwise_metrics(RNAi)
    write_phenodistinfo_to_CSV(RNAi, distAll, csi)
    # write_phenodistinfo_to_mySQL(RNAi, distAll, csi)

def initialize():

    run_phenodistcalc_populate()

if __name__ == '__main__':
    #     c = ControlClass()
    #     rGroups = []
    # #     rGroups.append([RNAiClass(i) for i in [9,5,55,6,435,59]])
    #     rGroups.append([RNAiClass(i) for i in [63,19,118,417,77]])
    # #     rGroups.append([RNAiClass(i) for i in [7,52,10,57,217]])
    # #     rGroups.append([RNAiClass(i) for i in [31,32,255]])
    #     optimizePAD(rGroups, c)

    #     pMeanRNAi, pStd = getParamsNorms()
    #     print('pMeanRNAi = ',pMeanRNAi)
    #     print('pStd = ',pStd)
    initialize()

    # RNAi = []
    # for i in range(1, 504):
    #     try:
    #         r = RNAiClass(i, verb=True)
    #         #             if r.__version__!=embdFunc.get_rnai_revision_version(): r.refreshParamsData()
    #         r.setNDPosition()
    #         RNAi.append(r)
    #     except Exception, e:
    #         printLog('!!!!!!!!!ERROR OCCURED RNAi {0} :{1}'.format(i, str(e)))
    #
    # pool = mp.Pool(processes=6)
    # distAll = pool.map(getPADS, zip(repeat(RNAi), range(len(RNAi))))
    # pool.close()
    # pool.join()
    # csi = getCSI(distAll)
    # fd = open('Z:\\phenoDist_rank_{d}.csv'.format(d=time.strftime('%m%d%Y')), 'w')
    # pUsed = ''
    # for p in PARAM_NAMES_USE_GLS:
    #     pUsed += '{0},'.format(p)
    # fd.write('params used GLS=,' + pUsed + '\n')
    # pUsed = ''
    # for p in PARAM_NAMES_USE_MS:
    #     pUsed += '{0},'.format(p)
    # fd.write('params used MS=,' + pUsed + '\n')
    # fd.write('RNA1,RNA2,dist,dist2WT,dims,CSI,PAD\n')
    # for i in range(len(RNAi)):
    #     for j in range(len(RNAi)):
    #         r1, r2, d, d20, nOv, pad = distAll[i][j]
    #         #             print(r1,r2,csi[i][j],d)
    #         fd.write('{0},{1},{2},{3},{4},{5},{6}\n'.format(r1, r2, d, d20, nOv, csi[i][j], pad))
