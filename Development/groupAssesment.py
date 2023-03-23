'''
Created on Mar 3, 2017

@author: Admin
'''

import myFunc
import numpy as np
from varLookup import printLog
from RNAiClass import RNAiClass, ControlClass
import PADassist
from myFigure import myFigure
import embdFunc
import dataValidation
import sklearn
from sklearn.model_selection import train_test_split
import multiprocessing as mp
from itertools import repeat
import json
from sklearn.manifold import TSNE
from params_use_assign import get_params_use

PARAM_NAMES_USE_GLS, PARAM_NAMES_USE_MS = get_params_use()


def getDist2WT(rnaList):
    dist = []
    for i in rnaList:
        r = RNAiClass(i)
        r.setNDPosition()
        d = r.getDist2Zero()
        dist.append(d)
    return np.array(dist)


def getDist2EO(rnaList):
    dist = []
    for i in range(len(rnaList) - 1):
        for j in range(i + 1, len(rnaList)):
            if isinstance(rnaList[i], object):
                r1 = rnaList[i]
                r2 = rnaList[j]
            else:
                r1 = RNAiClass(rnaList[i])
                r2 = RNAiClass(rnaList[j])
            r1.setNDPosition()
            r2.setNDPosition()
            d = r1.getDistance(r2)
            if not np.isnan(d) and not np.isinf(d) and d > 0: dist.append(d)
    return np.array(dist)


def getDist2AVG(rnaList):
    v, rs = [], []
    for i in range(len(rnaList)):
        r = RNAiClass(rnaList[i])
        r.setNDPosition()
        v.append(r.posND)
        rs.append(r)
    avgV = embdFunc.getNDAvgVector(v)
    dist = np.array([embdFunc.getNDDistance(r.posND, avgV, r.pNamesUse) for r in rs])
    return dist

def getAvgVector(rnaList):
    '''
    Takes in a list of embd numbers and returns the average vector position for that list
    :param rnaList: list of embd numbers
    :return: average ND vector
    '''
    v, rs = [], []
    for i in range(len(rnaList)):
        r = RNAiClass(rnaList[i])
        r.setNDPosition()
        v.append(r.posND)
        rs.append(r)
    avgV = embdFunc.getNDAvgVector(v)
    return avgV

def getPAD2EO(rnaList):
    '''
    calculates PAD values for all possible gene combinations within a list- excludes redundant measures between gene pairs i.e. A to B, B to A
    INPUT: rnaList is a list of RNAiClass objects
    OUTPUT: numpy array of PAD values for 
    '''
    pad = []
    for i in range(len(rnaList) - 1):
        for j in range(i + 1,
                       len(rnaList)):  # checks all possible combinations of RNAiClass objects (genes) within rnaList
            p = rnaList[i].getPAD(rnaList[j])  # gets the pad value for each combination of RNAiClass objects (genes)
            if not np.isnan(p) and not np.isinf(p):  # checks if not nan or infinity
                # pad.append(p)  # appends PAD values
                pad.append([rnaList[i].label,rnaList[j].label,p])
                # if p < 0.3:
                #     print('small pad', rnaList[i].RNAi, rnaList[j].RNAi)
    # print np.array(pad)
    return np.array(pad)  # array of all PAD values for genes within provided rnaList


def get_all_PAD2EO(rna_list):
    '''
    calculates PAD values for all possible gene combinations within a list- includes redundant measures between gene pairs
    INPUT: rnaList is a list of RNAiClass objects
    OUTPUT: numpy array of PAD values for
    '''
    pad = []
    for i in range(len(rna_list)):
        for j in range(len(rna_list)):
            p = rna_list[i].getPAD(rna_list[j])  # gets the pad value for each combination of RNAiClass objects (genes)
            if not np.isnan(p) and not np.isinf(p):  # checks if not nan or infinity
                pad.append(p)  # appends PAD values
    print np.array(pad)
    return np.array(pad)  # array of all PAD values for genes within provided rnaList


def get_all_dist2EO(rna_list):
    '''
    calculates dist values for all possible gene combinations within a list- includes redundant measures between gene pairs
    INPUT: rnaList is a list of RNAiClass objects
    OUTPUT: numpy array of dist values for all gene pairs
    '''
    dist = []
    for i in range(len(rna_list)):
        for j in range(len(rna_list)):
            p = rna_list[i].getDistance(
                rna_list[j])  # gets the dist value for each combination of RNAiClass objects (genes)
            if not np.isnan(p) and not np.isinf(p):  # checks if not nan or infinity
                dist.append(p)  # appends dist values
    # print np.array(pad)
    return np.array(dist)  # array of all dist values for genes within provided rnaList


def getCSI2EO(rnaList, data):
    csi = []
    for i in range(len(rnaList) - 1):
        for j in range(i + 1, len(rnaList)):
            c = PADassist.getCSIFromFile(rnaList[i], rnaList[j], data)
            if not np.isnan(c) and not np.isinf(c) and c > 0 and c < 0.99:
                csi.append(c)
    return np.array(csi)


def plotDist2WTvsArrest(rnaiList, labels):
    fig = myFigure()
    fig.markerSize = 8
    allDists = []
    allDistsErrs = []
    for i in range(len(rnaiList)):
        dist = getDist2WT(rnaiList[i])
        x = np.arange(-0.25, 0.25, 0.5 / dist.size) + i
        fig.scatter(x, dist)
        fig.errorbar(i, np.mean(dist), np.std(dist), color='k')
        allDists.append(np.mean(dist))
        allDistsErrs.append(np.std(dist))
        print('dist2WT {l} outlier {r}'.format(l=labels[i], r=rnaiList[i][np.argmax(dist)]))
    fig.xticks(range(len(rnaiList)), labels)
    #     fig.xlim([-1,6])
    fig.title('Dist to WT')
    fig.ylim((0, None))
    return fig, allDists, allDistsErrs


def plotDistInGroupVsArrest(rnaiList, labels):
    allDists = []
    allDistsErr = []
    fig = myFigure()
    fig.markerSize = 8
    for i in range(len(rnaiList)):
        dist = getDist2AVG(rnaiList[i])
        x = np.arange(-0.25, 0.251, 0.5 / dist.size) + i
        if x.size > dist.size: x = x[:dist.size]
        fig.scatter(x, dist)
        #         print('{l} avgDist={d} std={e}'.format(l=labels[i], d=np.mean(dist), e=np.std(dist)))
        fig.errorbar(i, np.mean(dist), np.std(dist), color='k')
        allDists.append(np.mean(dist))
        allDistsErr.append(np.std(dist))
        print('{l} outlier {r}'.format(l=labels[i], r=rnaiList[i][np.argmax(dist)]))
    fig.xticks(range(len(rnaiList)), labels)
    #     fig.xlim([-1,6])
    fig.title('Dist to avg')
    fig.ylim((0, None))
    return fig, allDists, allDistsErr


def plotPADInGroup(rnaiList, labels, simGroups):
    use_distance_instead = False
    fig = myFigure()
    fig.markerSize = 8
    rnaiGroups = []
    sep = []
    allD = []
    allDc = []
    groupSpread = []
    for rs in rnaiList:
        rns = []
        for r in rs:
            rns.append(RNAiClass(r))
            rns[-1].setNDPosition()
        rnaiGroups.append(rns)
    for i in range(len(rnaiGroups)):

        if use_distance_instead:
            dist = getDist2EO(rnaiGroups[i])
        else:
            dist = getPAD2EO(rnaiGroups[i])

        x = np.arange(-0.25, 0.251, 0.5 / dist.size) + i
        x = x[:dist.size]
        fig.scatter(x, dist)
        distContr = []
        for m in range(len(rnaiGroups[i])):
            for j in range(len(rnaiGroups)):
                if int(rnaiGroups[j][0].label[-4:]) not in simGroups[i]:
                    for k in range(len(rnaiGroups[j])):
                        if i != j:
                            if use_distance_instead:
                                distContr.append(rnaiGroups[i][m].getDistance(rnaiGroups[j][k]))
                                if distContr[-1] == 1.: print('DIST=1', rnaiGroups[i][m].RNAi, rnaiGroups[j][k].RNAi)
                            else:
                                distContr.append(rnaiGroups[i][m].getPAD(rnaiGroups[j][k]))
                                if distContr[-1] == 1.: print('PAD=1', rnaiGroups[i][m].RNAi, rnaiGroups[j][k].RNAi)

                            #                     if i!=j: distContr.append(PADassist.getPADFromFile(rnaiList[i][m], rnaiList[j][k], data))
        distContr = np.array(distContr)
        x = np.arange(-0.25, 0.251, 0.5 / distContr.size) + i
        x = x[:distContr.size]
        fig.scatter(x, distContr, marker='+')
        print('{l} avgDist={d} std={e}'.format(l=labels[i], d=np.mean(dist), e=np.std(dist)))
        fig.errorbar(i, np.mean(dist), np.std(dist), color='k')
        fig.errorbar(i, np.mean(distContr), np.std(distContr), color='magenta')
        groupSpread.append(np.std(dist))
        sep.append(np.mean(dist) - np.median(distContr))
        allD.append(np.mean(dist))
        allDc.append(np.mean(distContr))
    fig.xticks(range(len(rnaiList)), labels)
    #     fig.xlim([-1,6])
    #     fig.ylim([-0.2,1.2])
    if use_distance_instead:
        fig.title('Dist separation={0}, {1}'.format(np.mean(sep), min(allD) - max(allDc)))
    else:
        fig.title('PAD separation={0}, {1}'.format(np.mean(sep), min(allD) - max(allDc)))
    # fig.title('PAD separation={0}'.format(min(allD)-max(allDc)))
    # fig.title('PAD separation={0}, group spread={1}'.format((np.mean(allD)-np.mean(allDc)), np.mean(groupSpread)))

    return fig


def plotCSIInGroupVsArrest(rnaiList, labels):
    fileName = 'Z:\\Automated_analysis\\phenoDistCalc\\phenoDist_06232017_CSIonPAD.csv'
    data = myFunc.readCSV(fileName, ',')[3:]
    data = np.array(data)
    fig = myFigure()
    fig.markerSize = 8
    for i in range(len(rnaiList)):
        dist = getCSI2EO(rnaiList[i], data)
        x = np.arange(-0.25, 0.251, 0.5 / dist.size) + i
        if x.size > dist.size: x = x[:dist.size]
        fig.scatter(x, dist)
        distContr = []
        for m in range(len(rnaiList[i])):
            for j in range(len(rnaiList)):
                for k in range(len(rnaiList[j])):
                    if i != j: distContr.append(PADassist.getCSIFromFile(rnaiList[i][m], rnaiList[j][k], data))
        distContr = np.array(distContr)
        x = np.arange(-0.25, 0.251, 0.5 / distContr.size) + i
        x = x[:distContr.size]
        fig.scatter(x, distContr, marker='+')
        print('{l} avgDist={d} std={e}'.format(l=labels[i], d=np.mean(dist), e=np.std(dist)))
        fig.errorbar(i, np.mean(dist), np.std(dist), color='k')
        fig.errorbar(i, np.mean(distContr), np.std(distContr), color='magenta')
    fig.xticks(range(len(rnaiList)), labels)
    fig.xlim([-1, 6])
    fig.ylim([-0.2, 1.2])
    fig.title('CSI')
    return fig


def getRNAiGroupsArrest():
    #     rnaiList = [[19,31,63,77,118,181,264,357], [9,397,357,181,53],[5,274,268,435,255,264], [7,52,54,390,447,398],[10,217,291,412,59,288,398],[177,195,239,330,61,241,433,462,493,390]]
    rnaiList = [[19, 31, 63, 77, 118, 181, 264], [9, 397, 357, 53], [5, 90, 274, 268, 435, 255], [7, 52, 54, 447],
                [10, 217, 291, 412, 15, 59, 288, 390], [177, 195, 239, 330, 61, 241, 433, 462, 493]]
    labels = ['CF', 'comma', 'rupture', '1.5', '2', '3']
    return rnaiList, labels


def getRNAiGroupsPhenoOld():
    ''' initial groups for first pass optimization '''
    rnaiList = [[10, 15, 271, 18, 28, 439, 54, 177, 291, 288, 219, 398, 7, 209], [21, 422, 387, 408],
                [181, 182, 184, 397, 185, 357, 13], [31, 32, 255, 388, 130], [34, 16, 356, 386],
                [63, 19, 77, 118, 441, 501], [5, 55, 6, 133, 114, 137, 359, 435], [9, 45, 141], [426, 1, 76, 197]]
    labels = ['two-fold', 'mixed fate', 'high red', 'epid def', 'exc red', 'mex', 'rupture', 'enclosure', 'no markers']
    return rnaiList, labels

def calculate_avg_vectors_man_groups():  # in progress
    '''
    Calculates the average vectors for all of the manual groups
    :return: prints average vectors to be stored and accessed when calculating nearest group
    '''

    rnaiList, labels, defect_label = get_manual_groups()

    avg_vect_list = []
    labels_list = []

    for i in range(len(rnaiList)):
        g_label = labels[i]
        v = getAvgVector(rnaiList[i])
        avg_vect_list.append(v)
        labels_list.append(g_label)

    # print(labels_list)
    # print(avg_vect_list)
    return labels_list, avg_vect_list


def get_manual_groups():
    '''
    stable list to access the 17 groups genes
    :return: rnaiList, labels, defect_label
    '''
    rnaiList = [[21, 364, 408, 130], [264, 386, 31, 32, 255, 388, 422, 118, 359], [417, 64, 115, 281, 77],
                [63, 19, 501, 77], [181, 182, 357], [184, 185, 154, 363, 117, 45, 447, 108], [34, 16, 31, 264],
                [9, 398, 435, 45, 52, 38],
                [95, 4, 57, 5, 277], [67, 90, 235, 403, 503, 261, 404, 25],
                [26, 453, 235, 327, 489, 379, 420, 225, 289, 261, 350],
                [10, 15, 217, 18, 28, 439, 177, 291, 209], [110, 142, 98, 101], # note that 91 and 186 were removed from dev delay because these were originally put in the 'test' group and 91 is bad
                [383, 495, 498, 414, 375, 396, 321], [387, 385, 422, 31], [320, 58, 125, 288],
                [426, 197, 1, 76]]  # 320 group is sim in gls but morph is dif. 387 is wnt group, 426 is no/low markers

    labels = ['eCFS_R', 'sect rupt', 'eMEX', 'lMEX', 'lCFS_R', 'hiNeur', 'dorsRupt', 'dorsBend', 'rupture', 'gRupt',
              'gMix', '2xCrunch', 'dev delay', 'wt', 'CFS_wnt', 'CFS_lGlY', 'no_markers']
    defect_label = ['cf', 'cf', 'cf', 'cf', 'cf', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'cf', 'cf', 'cf']
    return rnaiList,labels,defect_label

def getRNAiGroupsPheno():
    #     ''' carefully vetted extended groups generated'''
    #     rnaiList = [[426, 1, 76, 197], [21,  364, 408, 130],[264, 386,31,32,  255, 388, 422, 118, 359], [417, 64,115, 281, 77],\
    #                 [63, 19, 501, 77],[181, 182, 357],[184,185, 154, 363, 117, 45, 447, 108],[34, 16, 31, 264],[ 9, 398, 435, 45,  52,38],\
    #                 [ 95, 4, 57, 5, 277],[67,90,235,403, 503, 261, 404, 25],[26, 453, 235, 327, 489, 379, 420, 225, 289, 261, 350, 3],\
    #                 [10, 15, 217, 18 ,28, 439,  177, 291, 209]]
    #     simGroups = [[426,417,326],[21,264,417],[264,21,385,326,63,34],[417,63,426,21],[63,417,21,264],[181,184],[184,181,9],[34,264],[9,184],\
    #                  [95,184,9],[67,10,26],[26,67],[10,9]]
    #     labels = ['no mark','eCFS_R','sect rupt', 'eMEX', 'lMEX','lCFS_R','hiNeur','dorsRupt', 'dorsBend', 'rupture', 'gRupt', 'gMix', '2xCrunch']
    # REMOVED No Markers group
    # rnaiList = [[21, 364, 408, 130], [264, 386, 31, 32, 255, 388, 422, 118, 359], [417, 64, 115, 281, 77],
    #             [63, 19, 501, 77], [181, 182, 357], [184, 185, 154, 363, 117, 45, 447, 108], [34, 16, 31, 264],
    #             [9, 398, 435, 45, 52, 38],
    #             [95, 4, 57, 5, 277], [67, 90, 235, 403, 503, 261, 404, 25],
    #             [26, 453, 235, 327, 489, 379, 420, 225, 289, 261, 350],
    #             [10, 15, 217, 18, 28, 439, 177, 291, 209], [110, 142, 98, 101, 91, 186], # note that 91 is a bad member
    #             [383, 495, 498, 414, 375, 396, 321], [387, 385, 422, 31], [320, 58, 125, 288],
    #             [426, 197, 1, 76]]  # 320 group is sim in gls but morph is dif. 387 is wnt group, 426 is no/low markers

    rnaiList = [[21, 364, 408, 130], [264, 386, 31, 32, 255, 388, 422, 118, 359], [417, 64, 115, 281, 77],
                [63, 19, 501, 77], [181, 182, 357], [184, 185, 154, 363, 117, 45, 447, 108], [34, 16, 31, 264],
                [9, 398, 435, 45, 52, 38],
                [95, 4, 57, 5, 277], [67, 90, 235, 403, 503, 261, 404, 25],
                [26, 453, 235, 327, 489, 379, 420, 225, 289, 261, 350],
                [10, 15, 217, 18, 28, 439, 177, 291, 209], [110, 142, 98, 101], # note that 91 and 186 were removed from dev delay because these were originally put in the 'test' group and 91 is bad
                [383, 495, 498, 414, 375, 396, 321], [387, 385, 422, 31], [320, 58, 125, 288],
                [426, 197, 1, 76]]  # 320 group is sim in gls but morph is dif. 387 is wnt group, 426 is no/low markers

    labels = ['eCFS_R', 'sect rupt', 'eMEX', 'lMEX', 'lCFS_R', 'hiNeur', 'dorsRupt', 'dorsBend', 'rupture', 'gRupt',
              'gMix', '2xCrunch', 'dev delay', 'wt', 'CFS_wnt', 'CFS_lGlY', 'no_markers']
    defect_label = ['cf', 'cf', 'cf', 'cf', 'cf', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'cf', 'cf', 'cf']

    rlist_train = []
    rlist_test = []
    labels_test = []
    defect_label_test = []
    # np.random.seed(071076)
    np.random.seed(1)

    for j in range(len(rnaiList)):
        group = rnaiList[j]
        y = np.zeros_like(group)
        X_train, X_test, y_train, y_test = train_test_split(group, y, test_size=0.3)
        if len(X_train) > 2 and len(X_test) > 1:
            rlist_train.append(X_train)
            rlist_test.append(X_test)
            labels_test.append(labels[j])
            defect_label_test.append(defect_label[j])
        else:
            rlist_train.append(group)

    # simGroups = [[21, 264, 417], [264, 21, 385, 326, 63, 34], [417, 63, 21], [63, 417, 21, 264], [181, 184],
    #              [184, 181, 9], [34, 264], [9, 184],
    #              [95, 184, 9], [67, 10, 26], [26, 67, 110, 383], [10, 9], [110, 383, 26], [383, 110, 26], [387, 264,417, 320], [320, 264, 387, 26], [426, 21, 417]]


    simGroups_by_label = [['eCFS_R', 'sect rupt', 'eMEX'],
                          ['sect rupt', 'eCFS_R', 'CFS_wnt', 'no_markers', 'lMEX', 'dorsRupt'],
                          ['eMEX', 'lMEX', 'eCFS_R'], ['lMEX', 'eMEX', 'eCFS_R', 'sect rupt'], ['lCFS_R', 'hiNeur'],
                          ['hiNeur', 'lCFS_R', 'dorsBend'], ['dorsRupt', 'sect rupt'], ['dorsBend', 'hiNeur'],
                          ['rupture', 'hiNeur', 'dorsBend'], ['gRupt', '2xCrunch', 'gMix'],
                          ['gMix', 'gRupt', 'dev delay', 'wt'], ['2xCrunch', 'dorsBend'], ['dev delay', 'wt', 'gMix'],
                          ['wt', 'dev delay', 'gMix'], ['CFS_wnt', 'sect rupt', 'eMEX', 'CFS_lGlY'],
                          ['CFS_lGlY', 'sect rupt', 'CFS_wnt', 'gMix'], ['no_markers', 'eCFS_R', 'eMEX']]

    simGroups = []
    simGroups_train = []
    simGroups_test = []

    for sim_group in simGroups_by_label:  # loops through all groups within simGroups_by _label
        simGroups.append([])  # adds new sublist to simGroups list
        simGroups_train.append([])
        simGroups_test.append([])
        for sg_label in sim_group:  # checks the length of the simsGroup_by label, sg_label is name of
            # individual phenotypic group within the simGroup_by_label
            ind = np.where(np.array(labels) == sg_label)[0][0]
            simGroups[-1].append(rnaiList[ind][0])
            simGroups_train[-1].append(rlist_train[ind][0])
            ind = np.where(np.array(labels_test) == sg_label)[0]
            if ind.size > 0:
                simGroups_test[-1].append(rlist_test[ind[0]][0])

    return rlist_test, rlist_train, simGroups_test, simGroups_train, labels, labels_test, defect_label, defect_label_test, rnaiList, labels, simGroups
    # return rnaiList, labels, simGroups, defect_label


def get_rnai_groups_pheno_distinct():
    rnaiList = [[21, 364, 408, 130], [264, 386, 31, 32, 255, 388, 422, 118, 359], [417, 64, 115, 281, 77, 63, 19, 501],
                [181, 182, 357, 184, 185, 154, 363, 117, 45, 447, 108], [9, 398, 435, 45, 52, 38], [95, 4, 57, 5, 277],
                [26, 453, 235, 327, 489, 379, 420, 225, 289, 261, 350],
                [10, 15, 217, 18, 28, 439, 177, 291, 209], [110, 142, 98, 101, 91, 186],
                [383, 495, 498, 414, 375, 396, 321]]
    labels = ['eCFS_R', 'sect rupt', 'MEX', 'hiNeur', 'dorsBend', 'rupture', 'gMix', '2xCrunch', 'dev delay', 'wt']
    return rnaiList, labels, []

def retrieve_and_map_man_group_autodata():
    '''
    pulls all data associated with man group genes and maps conditions as either "non-morph" or "morph" (can alter mapping to continuous by group or other). Selects certain columns for consideration. Can be toggled to return test/train split data.
    # :return: dataframes for labels and considered data are returned
    '''
    from sklearn.feature_extraction import DictVectorizer
    import pandas as pd
    from db_utils_genes import get_sql_data_for_man_group_conditions
    data = get_sql_data_for_man_group_conditions()  # retrieves data from sql
    data_to_use = data[['man_group','maxG_GLS','aG_GLS','sY_GLS','fracR_GLS','fracG_GLS', 'CoM0Rtd_GLS_GLS','CoM1Rscale_GLS_GLS','MoI0GavgRes_GLS_GLS', 'MoI1Gscale_GLS_GLS','MoI1GstdTail_GLS_GLS','MoI0Rtd_GLS_GLS','MoI0RavgRes_GLS_GLS','MoI1Rtd_GLS_GLS', 'MoI0Ytd_GLS_GLS','MoI0YavgRes_GLS_GLS', 'MoI1Yscale_GLS_GLS','MoI1YavgRes_GLS_GLS', 'sG_MS', 'maxR_MS', 'aR_MS', 'sR_MS', 'maxHead_MS', 'tailHead_MS', 'devHead_MS', 'devLength_MS', 'tailLength_MS', 'CoM0Gtd_MS_MS', 'CoM0Gscale_MS_MS', 'CoM1Gscale_MS_MS', 'CoM1GavgRes_MS_MS', 'CoM0Rscale_MS_MS', 'CoM0RavgRes_MS_MS', 'CoM1Rtd_MS_MS', 'CoM1RavgRes_MS_MS', 'CoM1RstdTail_MS_MS', 'MoI0Gscale_MS_MS', 'MoI1Gscale_MS_MS', 'MoI1GstdTail_MS_MS', 'MoI1Rscale_MS_MS', 'MoI1RavgRes_MS_MS', 'MoI1RstdTail_MS_MS']]

    data[data == np.inf] = np.nan   #flags infinity values as NaN
    data.fillna(0, inplace=True) #replaces Nan values with zero... not a great idea, but just trying it

    labels = data[['man_group']]
    labels_mapping = {'WT-like': 0, 'dev delay': 1, 'gRupt, gMix': 2, 'gMix': 3, 'gRupt': 4, '2x Crunch': 5,
                      'rupture': 6, 'hiNeur': 7,
                      'hiNeur, dorsBend': 8, 'dorsBend': 9, 'dors rupt': 10, 'sect rupt, dors rupt': 11,
                      'sect rupt': 12,
                      'sect rupt, CFS_wnt': 13, 'sect rupt, dors rupt, CFS_wnt': 14, 'CFS_lGlY': 15, 'lCFS_R': 16,
                      'CFS_wnt': 17, 'lMEX': 18, 'eMEX, lMEX': 19, 'eMEX': 20, 'eCFS_R': 21, 'no_markers': 22}
    labels_mapping_G_M = {'WT-like': 0, 'dev delay': 0, 'gRupt, gMix': 2, 'gMix': 2, 'gRupt': 1, '2x Crunch': 1,
                      'rupture': 1, 'hiNeur': 3,
                      'hiNeur, dorsBend': 2, 'dorsBend': 1, 'dors rupt': 1, 'sect rupt, dors rupt': 2,
                      'sect rupt': 2,
                      'sect rupt, CFS_wnt': 2, 'sect rupt, dors rupt, CFS_wnt': 2, 'CFS_lGlY': 3, 'lCFS_R': 3,
                      'CFS_wnt': 3, 'lMEX': 3, 'eMEX, lMEX': 3, 'eMEX': 3, 'eCFS_R': 3, 'no_markers': 3}  # this maps to "0-WT, 1-Morph, 2-Transition, and 3-Germ defect classes"

    labels_mapping_morph_or_other = {'WT-like': 0, 'dev delay': 0, 'gRupt, gMix': 0, 'gMix': 0, 'gRupt': 0, '2x Crunch': 1,
                          'rupture': 1, 'hiNeur': 0,
                          'hiNeur, dorsBend': 1, 'dorsBend': 1, 'dors rupt': 1, 'sect rupt, dors rupt': 1,
                          'sect rupt': 1,
                          'sect rupt, CFS_wnt': 0, 'sect rupt, dors rupt, CFS_wnt': 0, 'CFS_lGlY': 0, 'lCFS_R': 0,
                          'CFS_wnt': 0, 'lMEX': 0, 'eMEX, lMEX': 0, 'eMEX': 0, 'eCFS_R': 0,
                          'no_markers': 0}  # this maps to "0-non morph, 1-Morph

    # labels['man_group_mapped_narrow'] = labels.man_group.map(labels_mapping)
    # labels['man_group_mapped_broad'] = labels.man_group.map(labels_mapping_G_M)
    labels['man_group_mapped_binary'] = labels.man_group.map(labels_mapping_morph_or_other)

    # data_to_use['man_group_mapped_narrow'] = data_to_use.man_group.map(labels_mapping)
    # data_to_use['man_group_mapped_broad'] = data_to_use.man_group.map(labels_mapping_G_M)
    data['man_group_mapped_binary'] = data_to_use.man_group.map(labels_mapping_morph_or_other)


    # data_to_use.to_csv('Z:\\Automated_analysis\\man_groups_sql.csv')
    # print("SAVED")

    print(labels.man_group_mapped_binary.value_counts())


    labels_mapped = labels[['man_group_mapped_binary']]
    considered_data = data[['man_group_mapped_binary','maxG_GLS','aG_GLS','sY_GLS','fracR_GLS','fracG_GLS', 'CoM0Rtd_GLS_GLS','CoM1Rscale_GLS_GLS','MoI0GavgRes_GLS_GLS', 'MoI1Gscale_GLS_GLS','MoI1GstdTail_GLS_GLS','MoI0Rtd_GLS_GLS','MoI0RavgRes_GLS_GLS','MoI1Rtd_GLS_GLS', 'MoI0Ytd_GLS_GLS','MoI0YavgRes_GLS_GLS', 'MoI1Yscale_GLS_GLS','MoI1YavgRes_GLS_GLS', 'sG_MS', 'maxR_MS', 'aR_MS', 'sR_MS', 'maxHead_MS', 'tailHead_MS', 'devHead_MS', 'devLength_MS', 'tailLength_MS', 'CoM0Gtd_MS_MS', 'CoM0Gscale_MS_MS', 'CoM1Gscale_MS_MS', 'CoM1GavgRes_MS_MS', 'CoM0Rscale_MS_MS', 'CoM0RavgRes_MS_MS', 'CoM1Rtd_MS_MS', 'CoM1RavgRes_MS_MS', 'CoM1RstdTail_MS_MS', 'MoI0Gscale_MS_MS', 'MoI1Gscale_MS_MS', 'MoI1GstdTail_MS_MS', 'MoI1Rscale_MS_MS', 'MoI1RavgRes_MS_MS', 'MoI1RstdTail_MS_MS']]

    # x_train, x_test, y_train, y_test = train_test_split(labels_mapped, considered_data, test_size=0.3)
    #
    # return(x_train, y_train, x_test, y_test)
    return (labels, considered_data)

def test_two_sample_ks_test (labels, considered_data):
    '''
    takes in a df with binary labels and performs two sample ks test on all considered params (one at a time). Prints values and generates a cum distribution plot for each.
    :param considered_data: dataframe where first column contains binary data labels
    :return:
    '''
    from scipy import stats
    from myFigure import myFigure
    import pandas as pd

    plotFlag = True
    saveFlag = True

    df = considered_data
    cols = df.columns  # list of param names
    non_morph_data = df[df['man_group_mapped_binary']==0] # data where man groups are either WT-like or CFS
    morph_only_data = df[df['man_group_mapped_binary']==1]  #data where man_groups are morph in nature (2-fold, rupture)
    output = []
    for i in range(len(cols)): #iterates through the columns to test the distributions of each param
        j = i+1  # first column is the label, so this enables skipping that column
        if j <42:
            param = cols[j]  #param being evaluated
        if j >1 and j<42:
            data1 = non_morph_data.iloc[:, j]  #selects the data
            data2 = morph_only_data.iloc[:, j]
            test = stats.ks_2samp(data1, data2)
            stat = np.around(test[0], decimals=3)
            pval = np.around(test[1], decimals=5)
            info = (param, stat, pval)
            output.append(info)
            print(info)
            if plotFlag:
                '''plot fig of distribution'''
                fig = myFigure()
                data1p = np.sort(data1)
                data2p = np.sort(data2)
                y = np.arange(data1p.size)
                y2 = np.arange(data2p.size)
                fig.plot(data1p, 1. * y / data1p.size, color='black', label='non-morph')
                fig.plot(data2p, 1. * y2 / data2p.size, color='red', label='morph')
                fig.title((param, pval))
                fig.legend(2)
                fig.save('Z:\\Automated_analysis\\ks_tests\\ks_{0}.svg'.format(param))
                # fig.show()
            else:
                pass
        else:
            pass
    if saveFlag:
        outputdf = pd.DataFrame(output)
        outputdf.to_csv('Z:\\Automated_analysis\\ks_tests\\ks_tests_nonmorph_morph.csv')
        print("SAVED")

def log_regression():
    from sklearn.linear_model import LogisticRegression
    from matplotlib import pyplot as plt
    from scipy.special import expit


    lab, considered_data = retrieve_and_map_man_group_autodata()
    df = considered_data
    # cols = df.columns  # list of param names
    # non_morph_data = df[df['man_group_mapped_binary'] == 0]  # data where man groups are either WT-like or CFS
    # morph_only_data = df[
    #     df['man_group_mapped_binary'] == 1]  # data where man_groups are morph in nature (2-fold, rupture)

    labels = df['man_group_mapped_binary']  # binary 0 non-morph, 1 morph
    # data = df['sG_MS'] # continuous
    data = df['MoI1RstdTail_MS_MS']
    X = np.array(data)

    X = X.reshape((94,1))  # requires reshape of array
    y = np.array(labels)
    y = y.reshape((94,1))
    clf = LogisticRegression(random_state=50).fit(X, y)

    plt.figure(1, figsize=(4, 3))
    plt.clf()
    plt.scatter(X.ravel(), y, color='black', zorder=20)
    X_test = np.linspace(-0.5, 1, 300)

    loss = expit(X_test * clf.coef_ + clf.intercept_).ravel()
    plt.plot(X_test, loss, color='red', linewidth=3)
    plt.axhline(.5, color='.5')

    plt.ylabel('y')
    plt.xlabel('X')
    plt.xticks(range(-1, 1))
    plt.yticks([0, 0.5, 1])
    plt.ylim(-.25, 1.25)
    plt.xlim(-1, 1)
    plt.legend(('Logistic Regression Model'),
               loc="lower right", fontsize='small')
    plt.tight_layout()
    plt.show()

    print(clf.coef_)
    print(clf.intercept_)
    print(clf.predict(X))
    print(clf.predict([1.7]))
    # print(clf.predict_proba(X))
    # print(clf.score(X,y))
    print("done")

'''towards data sci example'''
    # >> > from sklearn.datasets import load_iris
    # >> > from sklearn.linear_model import LogisticRegression
    # >> > X, y = load_iris(return_X_y=True)
    # >> > clf = LogisticRegression(random_state=0).fit(X, y)
    # >> > clf.predict(X[:2, :])
    # array([0, 0])
    # >> > clf.predict_proba(X[:2, :])
    # array([[9.8...e - 01, 1.8...e - 02, 1.4...e - 08],
    #        [9.7...e - 01, 2.8...e - 02, ...e - 08]])
    # >> > clf.score(X, y)
    # 0.97...




# def search_cv(x_train, y_train, x_test, y_test):
#     '''
#     search cross validation to find key parameters according to manual groups. IN PROGRESS!!
#     :param x_train:
#     :param y_train:
#     :param x_test:
#     :param y_test:
#     :return:
#     '''
#     from sklearn.model_selection import cross_val_score, GridSearchCV
#     from sklearn.ensemble import RandomForestRegressor, GradientBoostingClassifier
#     def evaluate(model, test_features, test_labels):
#         predictions = model.predict(test_features)
#         errors = abs(predictions - test_labels)
#         mape = 100 * np.mean(errors / test_labels)
#         accuracy = 100 - mape
#         print('Model Performance')
#         print('Average Error: {:0.4f} degrees.'.format(np.mean(errors)))
#         print('Accuracy = {:0.2f}%.'.format(accuracy))
#         return accuracy
#     # grid search
#     # model = GradientBoostingClassifier(n_estimators=30)
#     param_grid = {
#         'bootstrap': [True],
#         'max_depth': [80, 90, 100, 110],
#         'max_features': [2, 3],
#         'min_samples_leaf': [3, 4, 5],
#         'min_samples_split': [8, 10, 12],
#         'n_estimators': [100, 200, 300, 1000]}
#     # Create a based model
#     rf = RandomForestRegressor()
#     # Instantiate the grid search model
#     grid_search = GridSearchCV(estimator=rf, param_grid=param_grid,
#                                cv=2, n_jobs=-1, verbose=2)
#     # Fit the grid search to the data
#     grid_search.fit(y_train, x_train)
#     print("training complete!!!!!!!!!!!!!!!!")
#     grid_search.best_params_
#     print("search complete--------------------")
#     best_grid = grid_search.best_estimator_
#
#     # grid_accuracy = evaluate(best_grid, test_features, test_labels)
#     # parameters = {"kernel": ["linear", "rbf"], "C": [1, 2, 4], "gamma": [0.125, 0.25, 0.5, 1, 2, 4]}
#     # parameters = [
#     #     {"C": [1, 10, 100, 1000], "kernel": ["linear"]},
#     #     {"C": [1, 10, 100, 1000], "gamma": [0.001, 0.0001], "kernel": ["rbf"]},
#     # ]
#     # clf = GridSearchCV(model, param_grid=param_grid)
#     # grid_search = clf.fit(x_train, y_train)
#     print("Best score: %0.3f" % grid_search.best_score_)
#     print(grid_search.best_estimator_)
#
#     # best prarams
#     print('best params:', grid_search.best_params_)
#
#     print('-----grid search end------------')
#     print('on all train set')
#     #
#     # print('importance!!!')
#     # # get importance
#     # importance = grid_search.feature_importances_
#     # # summarize feature importance
#     # for i, v in enumerate(importance):
#     #     print('Feature: %0d, Score: %.5f' % (i, v))
#     # # plot feature importance
#     # best_estimator = RandomForestRegressor(bootstrap=True, criterion='mse', max_depth=110,
#     #        max_features=0.1, max_leaf_nodes=None, min_impurity_split=1e-07,
#     #        min_samples_leaf=3, min_samples_split=8,
#     #        min_weight_fraction_leaf=0.0, n_estimators=200, n_jobs=1,
#     #        oob_score=False, random_state=None, verbose=0, warm_start=False)
#     # scores = cross_val_score(estimator=best_estimator, X=x_train, y=y_train, cv=3, scoring='accuracy')
#
#     scores = cross_val_score(estimator=grid_search.best_estimator_, X=x_train, y=y_train, cv=3, scoring='accuracy')
#     print(scores.mean(), scores)
#     print('on test set')
#     # scores = cross_val_score(estimator=best_estimator, X=x_test, y=y_test, cv=3, scoring='accuracy')
#     scores = cross_val_score(estimator=grid_search.best_estimator_, X=x_test, y=y_test, cv=3, scoring='accuracy')
#     print(scores.mean(), scores)
#
# def testing_importances():
#
#     '''
#     IN PROGRESS
#     :return:
#     '''
#     import pandas as pd
#     import numpy as np
#     import json
#     import matplotlib.pyplot as plt
#     from sklearn.model_selection import train_test_split, cross_val_score
#     from sklearn.metrics import confusion_matrix
#     from sklearn.preprocessing import scale
#     from sklearn.ensemble import ExtraTreesClassifier
#     from db_utils_genes import get_sql_data_for_man_group_conditions
#
#     #get data
#     data = get_sql_data_for_man_group_conditions()  # retrieves data from sql
#
#     #clean data of nans and infs
#     data[data == np.inf] = np.nan  # flags infinity values as NaN
#     data.fillna(0, inplace=True)  # replaces Nan values with zero... not a great idea, but just trying it
#
#     #define columns to consider and map strings to numbers
#     labels = data[['man_group']]
#
#     labels_mapping = {'WT-like': 0, 'dev delay': 1, 'gRupt, gMix': 2, 'gMix': 3, 'gRupt': 4, '2x Crunch': 5,
#                       'rupture': 6, 'hiNeur': 7,
#                       'hiNeur, dorsBend': 8, 'dorsBend': 9, 'dors rupt': 10, 'sect rupt, dors rupt': 11,
#                       'sect rupt': 12,
#                       'sect rupt, CFS_wnt': 13, 'sect rupt, dors rupt, CFS_wnt': 14, 'CFS_lGlY': 15, 'lCFS_R': 16,
#                       'CFS_wnt': 17, 'lMEX': 18, 'eMEX, lMEX': 19, 'eMEX': 20, 'eCFS_R': 21, 'no_markers': 22}
#     labels['man_group_mapped'] = labels.man_group.map(labels_mapping)
#     print(labels.man_group_mapped.value_counts())
#     labels_mapped = labels[['man_group_mapped']]
#     considered_data = data[
#         ['maxG_GLS', 'aG_GLS', 'sY_GLS', 'fracR_GLS', 'fracG_GLS', 'CoM0Rtd_GLS_GLS', 'CoM1Rscale_GLS_GLS',
#          'MoI0GavgRes_GLS_GLS', 'MoI1Gscale_GLS_GLS', 'MoI1GstdTail_GLS_GLS', 'MoI0Rtd_GLS_GLS', 'MoI0RavgRes_GLS_GLS',
#          'MoI1Rtd_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'MoI0YavgRes_GLS_GLS', 'MoI1Yscale_GLS_GLS', 'MoI1YavgRes_GLS_GLS',
#          'sG_MS', 'maxR_MS', 'aR_MS', 'sR_MS', 'maxHead_MS', 'tailHead_MS', 'devHead_MS', 'devLength_MS',
#          'tailLength_MS', 'CoM0Gtd_MS_MS', 'CoM0Gscale_MS_MS', 'CoM1Gscale_MS_MS', 'CoM1GavgRes_MS_MS',
#          'CoM0Rscale_MS_MS', 'CoM0RavgRes_MS_MS', 'CoM1Rtd_MS_MS', 'CoM1RavgRes_MS_MS', 'CoM1RstdTail_MS_MS',
#          'MoI0Gscale_MS_MS', 'MoI1Gscale_MS_MS', 'MoI1GstdTail_MS_MS', 'MoI1Rscale_MS_MS', 'MoI1RavgRes_MS_MS',
#          'MoI1RstdTail_MS_MS']]
#
#     #split data for training
#     x_train, x_test, y_train, y_test = train_test_split(labels_mapped, considered_data, test_size=0.3)
#
#     # define classifier and fitting data
#     forest = ExtraTreesClassifier(random_state=1)
#     forest.fit(x_train, y_train)
#
#     # predict and get confusion matrix
#     y_pred = forest.predict(x_test)
#     cm = confusion_matrix(y_test, y_pred)
#     print(cm)
#
#     # Applying 10-fold cross validation
#     accuracies = cross_val_score(estimator=forest, X=x_train, y=y_train, cv=10)
#     print("accuracy (10-fold): ", np.mean(accuracies))
#
#     # Features importances
#     importances = forest.feature_importances_
#     std = np.std([tree.feature_importances_ for tree in forest.estimators_],
#                  axis=0)
#     indices = np.argsort(importances)[::-1]
#     feature_list = [labels_mapped.columns[indices[f]] for f in range(labels_mapped.shape[1])]  # names of features.
#     ff = np.array(feature_list)
#
#     # Print the feature ranking
#     print("Feature ranking:")
#
#     for f in range(labels_mapped.shape[1]):
#         print("%d. feature %d (%f) name: %s" % (f + 1, indices[f], importances[indices[f]], ff[indices[f]]))
#
#     # Plot the feature importances of the forest
#     plt.figure()
#     plt.rcParams['figure.figsize'] = [16, 6]
#     plt.title("Feature importances")
#     plt.bar(range(labels_mapped.shape[1]), importances[indices],
#             color="r", yerr=std[indices], align="center")
#     plt.xticks(range(labels_mapped.shape[1]), ff[indices], rotation=90)
#     plt.xlim([-1, labels_mapped.shape[1]])
#     plt.show()
#
#     ## The new additions to get feature importance to classes:
#
#     # To get the importance according to each class:
#
#     def class_feature_importance(X, Y, feature_importances):
#         N, M = labels_mapped.shape
#         X = scale(X)
#
#         out = {}
#         for c in set(Y):
#             out[c] = dict(
#                 zip(range(M), np.mean(X[Y == c, :], axis=0) * feature_importances)
#             )
#
#         return out
#
#     result = class_feature_importance(X, y, importances)
#     print (json.dumps(result, indent=4))

def remove_param_from_use(rnai, p_name, add):
    """
    Removes parameter from use in distance and pad calculations
    :param rnai: RNAiClass object
    :param p_name: parameter name (with '_GLS'/'_MS' indicator at the end)
    :param add: indicator whether parameter is removed (-1) or added (1)
    :return: returns modified RNAiObject
    """
    if add == 1:
        rnai.weights[p_name] = 1
        if p_name[-2:] == 'MS':
            rnai.pNamesUse[1].append(
                p_name[:-3])  # adds parameter back to pNamesUse to be considered for the next round
        else:
            rnai.pNamesUse[0].append(
                p_name[:-4])
    elif add == -1:
        rnai.weights[
            p_name] = np.nan  # assigns parameter (pN) to have a weight of nan (this effectively kicks it out)
        if p_name[-2:] == 'MS':  # checks if MS or GLS pN
            p_ntemp = p_name[0:-3]  # removes _MS from pN to allow for comparison with pNamesUse
            rnai.pNamesUse[1] = np.delete(rnai.pNamesUse[1], np.where(np.array(rnai.pNamesUse[
                                                                                   1]) == p_ntemp)).tolist()  # removes parameter from pNamesUse to accurately set n-dimensional position for RNAiClass object
            rnai.setNDPosition()  # re-sets the position in n-dimensional space for each RNAiClass object (without bad parameter)
        else:
            p_ntemp = p_name[0:-4]
            rnai.pNamesUse[0] = np.delete(rnai.pNamesUse[0],
                                          np.where(np.array(rnai.pNamesUse[0]) == p_ntemp)).tolist()
            rnai.setNDPosition()
    return rnai


def optimizeByGain(gain, rnaiGroupsIni):
    """
    checks improvements in gain function by removing one parameter at a time
    INPUT: gain value (determined in optimizeByGroups- a measure of separation between in and out groups) and rnaiGroupsIni
    (also determined in optimizeByGroups- is a list of lists of manually curated groups of RNAiClass objects (genes)
    OUTPUT: prints to console the gain information and worst parameter for each iteration 
    """

    i = 0
    p_ns = ['start']  # 'start' used so list is not empty initially
    p_n_exclude = []  # list of parameters already removed (assessed as 'bad' parameters)
    while i < 150 and len(p_ns) > 0:  # goes through 100 iterations
        gain_ini = gain(None,
                        rnaiGroupsIni)  # uses gain measurement to determine separation for in/out group genes in rnaiGroupsIni without removing any params (none)
        printLog('INITIAL GAIN:{0}'.format(gain_ini))  # prints the starting gain value
        gain_change_all = []
        p_ns = []
        for pN in PARAM_NAMES_USE_GLS:  # iterates through each parameter to determine if it improves separation
            if pN + '_GLS' not in p_n_exclude:
                #                 print('GLS checking', pN)
                sep1 = gain(pN + '_GLS', rnaiGroupsIni)
                if gain_ini < sep1:  # if removal of parameter (pN) gives a bigger separation (improves gain) when compared to initial separation
                    #             print('GLS Bad parameter {p} gain improvement={q}'.format(p=pN, q=sep1-sep0))
                    gain_change_all.append(
                        sep1 - gain_ini)  # appends the difference between considered and initial gain values
                    p_ns.append(pN + '_GLS')  # appends parameter to pNs- these are the potentially bad parameters
        # else: print('{q}=gain; GLS Good parameter {p}'.format(p=pN, q=sep1-sep0))
        for pN in PARAM_NAMES_USE_MS:  # same as above for MS
            if pN + '_MS' not in p_n_exclude:
                #                 print('MS checking', pN)
                sep1 = gain(pN + '_MS', rnaiGroupsIni)
                if gain_ini < sep1:
                    #             print('MS Bad parameter {p} gain improvement={q}'.format(p=pN, q=sep1-sep0))
                    gain_change_all.append(sep1 - gain_ini)
                    p_ns.append(pN + '_MS')
        # else: print('{q}=gain; MS Good parameter {p}'.format(p=pN, q=sep1-sep0))
        if len(p_ns) > 0:
            pN = p_ns[np.argmax(gain_change_all)]  # finds the index position for parameter (pN in pNs)which gives
            # the largest difference between initial and test case (sepAll)
            p_n_exclude.append(pN)  # appends the parameter that yields the best separation (the Worst parameter)
            printLog('RUN {r} Worst parameter {p}, improvement={g}'.format(r=i, p=pN, g=max(gain_change_all)))
            for rs in rnaiGroupsIni:
                for r in rs:
                    r.weights[pN] = np.nan  # re-assigns weighting for 'bad parameter' to nan
                    if pN[-2:] == 'MS':  # checks if MS or GLS param
                        p_ntemp = pN[0:-3]  # removes the _MS component of pN for comparison with pNamesUse list
                        r.pNamesUse[1] = np.delete(r.pNamesUse[1],
                                                   np.where(np.array(r.pNamesUse[1]) == p_ntemp)).tolist()
                        r.setNDPosition()  # re-sets the position in n-dimensional space for each RNAiClass object (without bad parameter)

                    elif pN[-3:] == 'GLS':
                        p_ntemp = pN[0:-4]
                        r.pNamesUse[0] = np.delete(r.pNamesUse[0],
                                                   np.where(np.array(r.pNamesUse[0]) == p_ntemp)).tolist()
                        r.setNDPosition()  # re-sets the position in n-dimensional space for each RNAiClass object (without bad parameter)
            i += 1
        else:
            printLog('OPTIMIZATION COMPLETE')
    return p_n_exclude, gain_ini


def get_rnai_groups_obj(rnaiList):
    """
    loads RNAiClass objects
    :param rnaiList: list of lists of rnai numbers
    :return: same shape lists of RNAiClass objects
    """

    rnai_groups_ini = []  # list of RNAi objects for all groups in rnaiList
    for rs in rnaiList:
        rns = []  # list of RNAi objects within a group
        for r in rs:
            rns.append(RNAiClass(r))  # adds RNAi object to rns list
            rns[-1].setNDPosition()  # takes the last RNAi added and sets its ND position
        rnai_groups_ini.append(
            rns)  # adds rns (individual group) to rnaiGroupsIni, which is a list of all groups (as objects)
    return rnai_groups_ini


def optimize_groups_by_min_max(rnaiList, labels, simGroups):
    """
    optimizes separation between median PAD value within group relative to median value with outside group.
    The groups are defined by rnaiList.
    INPUT:
    rnaiList: list of lists of rnai embd numbers for each group
    labels: group names
    simGroups: groups that are related to
    """

    rnai_groups_ini = get_rnai_groups_obj(rnaiList)  # list of RNAi objects for all groups in rnaiList

    def gain(pN, rnaiGroups):
        """this is a special gain function that is specific to optimizeByGroups, which weights considered parameter (pN) as nan to determine
        whether this parameter is helping or hurting with respect to separation of in vs. out group within rnaiGroups
        INPUT: pN = individual parameter, rnaiGroups= a list of lists containing manually curated groups of RNAiClass objects (genes)
        OUTPUT= gain value *may change depending on algorithm used
        """

        sep = []
        allD = []
        allDc = []
        if pN != None:  # if pN exists
            for i in range(len(rnaiGroups)):  # rs is a group of RNAiClass objects (genes) in rnaiGroups
                for j in range(len(rnaiGroups[i])):  # r is an RNAiClass object (gene) from rs
                    rnaiGroups[i][j] = remove_param_from_use(rnaiGroups[i][j], pN, -1)
        for i in range(len(rnaiGroups)):  # for each group
            dist = getPAD2EO(rnaiGroups[
                                 i])  # calculates PAD values for all possible gene combinations within ith group (1 group: a list of genes)
            distContr = []  # this is outgroup
            for m in range(len(rnaiGroups[i])):  # m is a RNAiClass object (gene) within the ith group
                for j in range(len(rnaiGroups)):  # j is the index for 1 group of RNAiClass objects (genes)
                    if int(rnaiGroups[j][0].label[-4:]) not in simGroups[
                        i]:  # if first gene within j group is not in similarity group i
                        for r in rnaiGroups[j]:  # r is an RNAiClass object from jth group in rnaiGroup
                            if i != j: distContr.append(rnaiGroups[i][m].getPAD(
                                r))  # appends PAD values for in group gene (rnaiGroups[i][m]) vs. out of group gene (r)
            distContr = np.array(distContr)
            sep.append(np.mean(dist) - np.mean(distContr))  # appends the difference of median for in vs. out group PAD
            allD.append(np.mean(dist))  # appends median in group PAD
            allDc.append(np.mean(distContr))  # appends median out group PAD
        if pN != None:  # returns the weight for the considered parameter to 1
            for i in range(len(rnaiGroups)):  # rs is a group of RNAiClass objects (genes) in rnaiGroups
                for j in range(len(rnaiGroups[i])):  # r is an RNAiClass object (gene) from rs
                    rnaiGroups[i][j] = remove_param_from_use(rnaiGroups[i][j], pN, 1)
                    #         return np.mean(sep)#gain value that indicates the separation of in and out group
                    #         return np.mean(allD)-np.mean(allDc) #gain value that indicates the difference between the average in-group and average out-group
        return min(allD) - max(allDc)
        # gain value that indicates the difference between the lowest in group and highest out group averages

    optimizeByGain(gain, rnai_groups_ini)


def optimizeByInteractome():
    '''
    optimizes F1 metric when PAD prediction is compared with interactome.
    The groups are defined by rnaiList.
    '''

    interact = dataValidation.getInteractionsIP([embdFunc.getGeneName(i) for i in range(1, 504)], interactome=True,
                                                phenome=False)
    allGenes = interact.ravel()
    allGenes = np.unique(allGenes).astype(np.int)
    rnaList = [RNAiClass(i) for i in allGenes]
    for r in rnaList: r.setNDPosition()

    def gain(pN):
        if pN != None:
            for r in rnaList:
                r.weights[pN] = 0
                r.setNDPosition()
        allPads = []
        for i in range(len(rnaList) - 1):
            for j in range(i + 1, len(rnaList)):
                allPads.append((rnaList[i].RNAi, rnaList[j].RNAi, rnaList[i].getPAD(rnaList[j])))
        allPads = np.array(allPads)
        f1, rec, pr, x = dataValidation.getF1metricsAcrossThr(allPads, interact)
        if pN != None:
            for r in rnaList:
                r.weights[pN] = 1
        return max(f1)

    optimizeByGain(gain)


def flatten_groups(list_of_lists):
    """
    Flattens list of lists
    :param list_of_lists: list of lists
    :return: list
    """
    list_flat = []
    for l in list_of_lists:
        list_flat += l
    return list_flat


def optimize_groups_by_rank(args):
    """
    optimizes ranking within each group.
    The groups are defined by rnaiList.
    INPUT:
    args: tuple of the following arguments
    rnaiList: list of lists of rnai embd numbers for each group
    labels: group names
    simGroups: groups that are related to
    groupExclude: group number to be excluded from computation
    OUTPUT:
    tuple: the excluded group number, the excluded parameters and the gain
    """
    rnai_list, labels, sim_groups, defect_label, group_exclude = args
    if group_exclude is not None:
        rnai_list.pop(group_exclude)
        labels.pop(group_exclude)
        sim_groups.pop(group_exclude)
        defect_label.pop(group_exclude)
    defect_label = np.array(defect_label)
    rnai_groups_ini = get_rnai_groups_obj(rnai_list)

    def gain(pN, rnai_groups):
        """this is a special gain function that is specific to optimizeByGroups, which weights considered parameter (pN) as nan to determine
        whether this parameter is helping or hurting with respect to separation of in vs. out group within rnaiGroups
        INPUT: pN = individual parameter, rnaiGroups= a list of lists containing manually curated groups of RNAiClass objects (genes)
        OUTPUT= gain value *may change depending on algorithm used
        """

        if pN is not None:  # if pN exists
            for i in range(len(rnai_groups)):  # rs is a group of RNAiClass objects (genes) in rnaiGroups
                for j in range(len(rnai_groups[i])):  # r is an RNAiClass object (gene) from rs
                    rnai_groups[i][j] = remove_param_from_use(rnai_groups[i][j], pN, -1)
        all_pad_dist = calc_pad_normdist(flatten_groups(rnai_groups))
        group_ranks_pad = []
        for r_in_group in rnai_list:  # for each group (use rnai_list to get int vals)
            group_ranks_pad.append(get_in_group_rank(r_in_group, all_pad_dist, 2))
        if pN is not None:  # returns the weight for the considered parameter to 1
            for i in range(len(rnai_groups)):  # rs is a group of RNAiClass objects (genes) in rnaiGroups
                for j in range(len(rnai_groups[i])):  # r is an RNAiClass object (gene) from rs
                    rnai_groups[i][j] = remove_param_from_use(rnai_groups[i][j], pN, 1)
        group_ranks_pad = np.array(group_ranks_pad)
        ind_cf = np.where(defect_label == 'cf')
        ind_m = np.where(defect_label == 'm')

        group_ranks_pad[ind_cf] /= 1. * group_ranks_pad[ind_cf].size
        group_ranks_pad[ind_m] /= 1. * group_ranks_pad[ind_m].size
        return -np.sum(group_ranks_pad) / 2.  # minus sign due to high rank being a disadvantage
        # return -np.mean(group_ranks_pad)  # minus sign due to high rank being a disadvantage

    p_exclude, final_gain = optimizeByGain(gain, rnai_groups_ini)
    return group_exclude, p_exclude, final_gain


def calc_pad_normdist(rna_list_obj):
    """
    calculate PAD and normalized distance pairwise for list of RNAs (rnaList is list of ints)
    :param rna_list_obj: list of RNAiClass objects
    :return: list of lists structured as label 1 (int), label 2 (int), PAD (float), normalized distance (float)
    """
    data = []
    for i in range(len(rna_list_obj) - 1):
        r1 = rna_list_obj[i]
        for j in range(i + 1, len(rna_list_obj)):
            r2 = rna_list_obj[j]
            pad = r1.getPAD(r2)
            ndist = r1.get_norm_dist(r2)
            cos_dist = r1.get_cosine_dist(r2)
            # ndist = 0  # FIXME: removed distance calculation for speed. Change back to r1.get_norm_dist(r2)
            if (r1.RNAi, r2.RNAi, pad, ndist, cos_dist) not in data:
                data.append((r1.RNAi, r2.RNAi, pad, ndist, cos_dist))
            if (r2.RNAi, r1.RNAi, pad, ndist, cos_dist) not in data:
                data.append((r2.RNAi, r1.RNAi, pad, ndist, cos_dist))
    data.sort(key=lambda x: (x[0], x[1]))
    return np.array(data)


def get_top_hits(rnai, all_pad_dist, index):
    """
    Gets top hits using index position
    :param rnai: reference rnai (int)
    :param all_pad_dist: list of lists with rnai pairs, pads and dist
    :param index: which position to sort all_pad_dist by
    :return: returns list of lists with top_n hits for reference rnai
    """
    sub_rnai_list = all_pad_dist[all_pad_dist[:, 0] == rnai]
    inds = np.argsort(sub_rnai_list[:, index])
    return sub_rnai_list[inds[::-1]]


def get_in_group_rank(rnai_list, all_pad_dist, index):
    """
    Calculates average rank for a list of rnai
    :param rnai_list: list of rnais
    :param all_pad_dist: list of lists with rnai pairs, pads and dist
    :param index: index position in all_pad_dist to use for ranking
    :return: average rank divided by the smallest possible sum of ranks for the group (arithmetic progression:n*(n-1)/2,
    where n is total number of rnais in the group)
    Note: number of ranks summed for each seed is n-1 and ranks start from 1, so arithmetic sum is n*(n-1)/2,
    another factor of 1/n comes from averaging.
    """
    rank = 0.
    for i in range(len(rnai_list)):
        top_hits = get_top_hits(rnai_list[i], all_pad_dist, index)
        for j in range(len(rnai_list)):
            if j != i:
                rank += np.where(top_hits[:, 1] == rnai_list[j])[0][0] + 1
    return 2. * rank / len(rnai_list) ** 2 / (len(rnai_list) - 1)


def evaluate_ranking(rlist_test, labels, defect_label):
    from random import shuffle
    # rnai_list, labels, sim_groups, defect_label = getRNAiGroupsPheno()
    # test with a random list
    # flat_rnai = [item for sublist in rnai_list for item in sublist]
    # shuffle(flat_rnai)
    # rand_list = []
    # i = 0
    # for sub in rnai_list:
    #     rand_list.append(flat_rnai[i:i+len(sub)])
    #     i += len(sub)
    # rnai_list = rand_list

    rnai_list_obj = get_rnai_groups_obj(rlist_test)
    rnai_list_flat = flatten_groups(rnai_list_obj)
    all_pad_dist = calc_pad_normdist(rnai_list_flat)
    group_ranks_pad = []
    group_ranks_dist = []
    group_ranks_cos_dist = []
    for r_in_group in rlist_test:
        group_ranks_pad.append(get_in_group_rank(r_in_group, all_pad_dist, 2))
        group_ranks_dist.append(get_in_group_rank(r_in_group, all_pad_dist, 3))
        group_ranks_cos_dist.append(get_in_group_rank(r_in_group, all_pad_dist, 4))

    for i in range(len(labels)):
        print '{l}, {p}, {d}, {c}'.format(l=labels[i], p=group_ranks_pad[i], d=group_ranks_dist[i], c=group_ranks_cos_dist[i])
    group_ranks_pad = np.array(group_ranks_pad)
    group_ranks_dist = np.array(group_ranks_dist)
    group_ranks_cos_dist = np.array(group_ranks_cos_dist)
    defect_label = np.array(defect_label)
    ind_cf = np.where(defect_label == 'cf')
    ind_m = np.where(defect_label == 'm')

    group_ranks_pad[ind_cf] /= 1. * group_ranks_pad[ind_cf].size
    group_ranks_pad[ind_m] /= 1. * group_ranks_pad[ind_m].size
    group_ranks_dist[ind_cf] /= 1. * group_ranks_dist[ind_cf].size
    group_ranks_dist[ind_m] /= 1. * group_ranks_dist[ind_m].size
    group_ranks_cos_dist[ind_cf] /= 1. * group_ranks_cos_dist[ind_cf].size
    group_ranks_cos_dist[ind_m] /= 1. * group_ranks_cos_dist[ind_m].size

    print 'AVG, {p}, {d}, {c}'.format(p=np.sum(group_ranks_pad) / 2., d=np.sum(group_ranks_dist) / 2., c=np.sum(group_ranks_cos_dist) / 2.)


def cross_corr_params(th):
    # '''this was from the 3/7/18 group opt by rank run, using training set (70% training, with at least 3 members to train) for RNAi list. re-set the origin prior to running this'''

    opt_res = [(None, ['CoM0RstdTail_GLS_GLS', 'mG_MS', 'tMove_MS', 'scaleHead_MS', 'CoM0Yscale_GLS_GLS', 'scaleLength_MS', 'sG_GLS', 'CoM1YstdTail_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM1RstdTail_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'MoI1GavgRes_MS_MS', 'mRmG_MS', 'MoI0Gscale_GLS_GLS', 'MoI0RstdTail_GLS_GLS', 'mYmG_GLS', 'tScale_MS', 'mR_MS', 'MoI1GavgRes_GLS_GLS', 'Rtail_GLS', 'MoI0Rscale_MS_MS', 'MoI0Rtd_MS_MS', 'CoM1Gtd_MS_MS', 'CoM0YavgRes_GLS_GLS', 'tScale_GLS', 'Gtail_GLS', 'CoM0RavgRes_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'mRmG_GLS', 'MoI1RstdTail_GLS_GLS', 'CoM1Yscale_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI0Gtd_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'aY_GLS', 'CoM1GstdTail_MS_MS', 'CoM1Rscale_MS_MS', 'MoI1Ytd_GLS_GLS', 'Ytail_GLS', 'MoI0Rtd_GLS_GLS', 'MoI1Rtd_MS_MS', 'CoM0Gscale_MS_MS', 'CoM0GavgRes_MS_MS', 'CoM1RavgRes_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'CoM0RstdTail_MS_MS', 'MoI0GstdTail_GLS_GLS', 'aR_GLS', 'aG_MS', 'MoI0Rscale_GLS_GLS', 'MoI0RstdTail_MS_MS', 'MoI1Gtd_MS_MS', 'MoI0YstdTail_GLS_GLS', 'MoI0GstdTail_MS_MS', 'maxR_GLS'], -2.5013766061980349), (0, ['CoM0RstdTail_GLS_GLS', 'mG_MS', 'tMove_MS', 'scaleHead_MS', 'MoI0GavgRes_MS_MS', 'scaleLength_MS', 'CoM1YstdTail_GLS_GLS', 'sG_GLS', 'CoM0Yscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'MoI1GavgRes_MS_MS', 'tScale_MS', 'aG_MS', 'CoM1Rtd_GLS_GLS', 'tScale_GLS', 'MoI1YstdTail_GLS_GLS', 'Rtail_GLS', 'mRmG_GLS', 'Ytail_GLS', 'aY_GLS', 'mR_MS', 'CoM0Rscale_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'MoI0RstdTail_GLS_GLS', 'aR_GLS', 'MoI0Rtd_MS_MS', 'MoI1Gtd_MS_MS', 'Gtail_GLS', 'MoI1RstdTail_GLS_GLS', 'tMove_GLS', 'CoM0Ytd_GLS_GLS', 'mYmG_GLS', 'MoI0GstdTail_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'mRmG_MS', 'CoM0RavgRes_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'fracY_GLS', 'CoM1Rtd_MS_MS', 'MoI1RavgRes_GLS_GLS', 'CoM0GstdTail_MS_MS', 'sR_GLS', 'CoM0YstdTail_GLS_GLS', 'CoM1GstdTail_MS_MS', 'MoI0YstdTail_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI0Gtd_MS_MS', 'CoM1Yscale_GLS_GLS', 'MoI0RavgRes_MS_MS'], -2.5670534139581758), (1, ['CoM0RstdTail_GLS_GLS', 'mG_MS', 'tMove_MS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'scaleLength_MS', 'CoM0Yscale_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'sG_GLS', 'mRmG_MS', 'CoM1RstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'mRmG_GLS', 'MoI0Rscale_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'Gtail_GLS', 'CoM0RavgRes_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'aY_GLS', 'mR_MS', 'tScale_GLS', 'fracY_GLS', 'maxR_GLS', 'CoM1Gtd_MS_MS', 'CoM0YavgRes_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI1YstdTail_GLS_GLS', 'MoI0RstdTail_MS_MS', 'MoI0Rtd_MS_MS', 'tScale_MS', 'Ytail_GLS', 'maxG_MS', 'MoI1RstdTail_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI0GavgRes_MS_MS', 'CoM1Rtd_GLS_GLS', 'MoI1RavgRes_GLS_GLS', 'aR_GLS', 'CoM0GavgRes_MS_MS', 'Rtail_GLS', 'MoI1GavgRes_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'CoM1GstdTail_MS_MS', 'CoM0YstdTail_GLS_GLS', 'mYmG_GLS', 'MoI1Gtd_MS_MS', 'CoM0RstdTail_MS_MS', 'MoI0Gtd_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'sR_GLS', 'aG_MS', 'MoI0GstdTail_MS_MS', 'MoI0Rscale_MS_MS', 'CoM0Ytd_GLS_GLS', 'MoI0RstdTail_GLS_GLS', 'CoM1RstdTail_MS_MS', 'MoI0GstdTail_GLS_GLS', 'CoM1Yscale_GLS_GLS', 'CoM1GavgRes_MS_MS', 'devHead_MS'], -2.1917346938775513), (2, ['CoM0RstdTail_GLS_GLS', 'mG_MS', 'tMove_MS', 'scaleHead_MS', 'CoM0Yscale_GLS_GLS', 'scaleLength_MS', 'sG_GLS', 'CoM1YstdTail_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM1RstdTail_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'mRmG_MS', 'MoI1Rscale_GLS_GLS', 'MoI1GavgRes_MS_MS', 'Gtail_GLS', 'MoI1YstdTail_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'Rtail_GLS', 'CoM0RavgRes_GLS_GLS', 'tScale_MS', 'CoM0YstdTail_GLS_GLS', 'tScale_GLS', 'mRmG_GLS', 'MoI1RstdTail_GLS_GLS', 'aY_GLS', 'sR_GLS', 'Ytail_GLS', 'CoM1Rscale_MS_MS', 'fracY_GLS', 'MoI0Ytd_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'aR_GLS', 'MoI0Rtd_MS_MS', 'maxR_GLS', 'MoI0Gtd_MS_MS', 'CoM0GavgRes_MS_MS', 'MoI0RstdTail_MS_MS', 'mYmG_GLS', 'maxG_MS', 'CoM1Rtd_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'MoI0Rscale_MS_MS', 'CoM1RstdTail_MS_MS', 'MoI0GstdTail_MS_MS', 'MoI0Rscale_GLS_GLS', 'MoI1RavgRes_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'CoM1Gtd_MS_MS', 'MoI1Gtd_MS_MS', 'tailLength_MS', 'tMove_GLS', 'MoI0RstdTail_GLS_GLS', 'CoM0RstdTail_MS_MS', 'mR_MS', 'CoM1GstdTail_MS_MS', 'CoM0GstdTail_MS_MS'], -2.3921598639455786), (3, ['mG_MS', 'tMove_MS', 'CoM0RstdTail_GLS_GLS', 'scaleHead_MS', 'MoI0GavgRes_MS_MS', 'mRmG_MS', 'CoM0Yscale_GLS_GLS', 'MoI1GavgRes_MS_MS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'sG_GLS', 'sR_GLS', 'mR_MS', 'tScale_MS', 'CoM0YstdTail_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'Rtail_GLS', 'Ytail_GLS', 'CoM1Ytd_GLS_GLS', 'tScale_GLS', 'scaleLength_MS', 'CoM0Rscale_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'MoI0RstdTail_GLS_GLS', 'MoI1RavgRes_GLS_GLS', 'mYmG_GLS', 'aY_GLS', 'MoI0Gscale_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'CoM1Gtd_MS_MS', 'Gtail_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI0Rtd_MS_MS', 'CoM0Rtd_MS_MS', 'CoM1YavgRes_GLS_GLS', 'MoI0RstdTail_MS_MS', 'CoM1Rscale_MS_MS', 'CoM0GavgRes_MS_MS', 'maxG_MS', 'MoI0Gtd_MS_MS', 'CoM0GstdTail_MS_MS', 'mRmG_GLS', 'MoI0Gtd_GLS_GLS', 'MoI1GstdTail_MS_MS', 'CoM0Ytd_GLS_GLS', 'maxG_GLS', 'MoI0RavgRes_MS_MS', 'MoI0GstdTail_MS_MS', 'CoM0RstdTail_MS_MS', 'MoI1YavgRes_GLS_GLS', 'CoM0Rscale_MS_MS', 'MoI1Yscale_GLS_GLS', 'MoI1Ytd_GLS_GLS'], -2.2232501889644745), (4, ['CoM0RstdTail_GLS_GLS', 'mG_MS', 'tMove_MS', 'CoM0Yscale_GLS_GLS', 'scaleHead_MS', 'scaleLength_MS', 'sG_GLS', 'mRmG_MS', 'MoI1GavgRes_MS_MS', 'CoM1RstdTail_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'mR_MS', 'Gtail_GLS', 'MoI1YstdTail_GLS_GLS', 'tScale_MS', 'CoM1Rtd_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'MoI0RstdTail_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM0RavgRes_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'Rtail_GLS', 'CoM0YstdTail_GLS_GLS', 'maxR_GLS', 'MoI0Rscale_MS_MS', 'MoI1GavgRes_GLS_GLS', 'mRmG_GLS', 'tScale_GLS', 'CoM1RavgRes_GLS_GLS', 'aY_GLS', 'aR_GLS', 'mYmG_GLS', 'Ytail_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'CoM1Rscale_MS_MS', 'CoM1Gtd_MS_MS', 'MoI0GstdTail_GLS_GLS', 'sR_GLS', 'CoM1Yscale_GLS_GLS', 'MoI1Gtd_MS_MS', 'MoI0Gtd_GLS_GLS', 'aG_MS', 'CoM1RstdTail_MS_MS', 'MoI1Rtd_MS_MS', 'MoI0YstdTail_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'CoM1GstdTail_MS_MS', 'CoM0Rtd_MS_MS'], -2.5757533383723863), (5, ['CoM0RstdTail_GLS_GLS', 'mG_MS', 'tMove_MS', 'scaleHead_MS', 'MoI0GavgRes_MS_MS', 'CoM0Yscale_GLS_GLS', 'scaleLength_MS', 'sG_GLS', 'MoI1Rscale_GLS_GLS', 'MoI1GavgRes_MS_MS', 'MoI0Yscale_GLS_GLS', 'mRmG_MS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'Gtail_GLS', 'tScale_GLS', 'mR_MS', 'CoM1Ytd_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'mRmG_GLS', 'Rtail_GLS', 'MoI1GavgRes_GLS_GLS', 'Ytail_GLS', 'MoI0Rscale_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI1RavgRes_GLS_GLS', 'maxG_MS', 'CoM1Gtd_MS_MS', 'MoI1RstdTail_GLS_GLS', 'maxR_GLS', 'aY_GLS', 'CoM0RavgRes_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'mYmG_GLS', 'MoI1Ytd_GLS_GLS', 'CoM0Rtd_MS_MS', 'MoI1Gtd_GLS_GLS', 'MoI0Ytd_GLS_GLS'], -2.663158304988662), (6, ['CoM0RstdTail_GLS_GLS', 'scaleHead_MS', 'CoM0Yscale_GLS_GLS', 'MoI1GavgRes_MS_MS', 'tMove_MS', 'scaleLength_MS', 'CoM1YstdTail_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'sG_GLS', 'MoI1YstdTail_GLS_GLS', 'mRmG_MS', 'MoI0Yscale_GLS_GLS', 'mR_MS', 'MoI1Rscale_GLS_GLS', 'Gtail_GLS', 'tScale_MS', 'Ytail_GLS', 'mRmG_GLS', 'maxR_GLS', 'CoM1Ytd_GLS_GLS', 'tScale_GLS', 'CoM0YstdTail_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'CoM1Rscale_MS_MS', 'CoM1Yscale_GLS_GLS', 'Rtail_GLS', 'MoI0GavgRes_MS_MS', 'CoM0Rscale_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI0Rtd_MS_MS', 'mG_MS', 'CoM0RavgRes_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'CoM0GavgRes_MS_MS', 'MoI1GavgRes_GLS_GLS', 'sR_GLS', 'MoI1GstdTail_MS_MS', 'MoI0Gtd_GLS_GLS', 'MoI0Gtd_MS_MS', 'mYmG_GLS', 'MoI0RstdTail_GLS_GLS', 'aY_GLS', 'CoM0Ytd_GLS_GLS', 'CoM1Gtd_MS_MS', 'MoI0GstdTail_MS_MS', 'aG_MS', 'MoI0Rscale_MS_MS', 'MoI0Rscale_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'maxG_MS', 'CoM1RavgRes_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'aR_GLS', 'CoM0RstdTail_MS_MS', 'CoM0GstdTail_MS_MS'], -2.4233475056689344), (7, ['CoM0RstdTail_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'mG_MS', 'tMove_MS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'sG_GLS', 'mRmG_MS', 'CoM1YstdTail_GLS_GLS', 'scaleLength_MS', 'CoM1RstdTail_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'mR_MS', 'Gtail_GLS', 'MoI1YstdTail_GLS_GLS', 'tScale_MS', 'CoM1Rtd_GLS_GLS', 'Rtail_GLS', 'CoM0Rscale_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'MoI0GavgRes_MS_MS', 'tScale_GLS', 'mRmG_GLS', 'MoI1GavgRes_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'aY_GLS', 'CoM1Rscale_MS_MS', 'Ytail_GLS', 'MoI0GavgRes_GLS_GLS', 'CoM0GstdTail_MS_MS', 'MoI1Gtd_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'CoM0GavgRes_MS_MS', 'mYmG_GLS', 'CoM1Ytd_GLS_GLS', 'CoM1GstdTail_MS_MS', 'MoI0Gtd_GLS_GLS', 'MoI1Rtd_MS_MS', 'MoI0Rtd_MS_MS', 'MoI0Gtd_MS_MS', 'MoI1Ytd_GLS_GLS', 'MoI1Gtd_MS_MS', 'MoI0YstdTail_GLS_GLS'], -2.3298037131519274), (8, ['mG_MS', 'tMove_MS', 'CoM0RstdTail_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'mRmG_MS', 'CoM1RstdTail_GLS_GLS', 'sG_GLS', 'MoI0Yscale_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'mR_MS', 'CoM1YstdTail_GLS_GLS', 'Gtail_GLS', 'CoM0Rscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'MoI1RstdTail_GLS_GLS', 'mYmG_GLS', 'tScale_MS', 'aY_GLS', 'mRmG_GLS', 'tScale_GLS', 'CoM1Gtd_MS_MS', 'sR_GLS', 'MoI0RavgRes_GLS_GLS', 'fracY_GLS', 'CoM0Ytd_GLS_GLS', 'Ytail_GLS', 'CoM1Ytd_GLS_GLS', 'Rtail_GLS', 'MoI0Gtd_GLS_GLS', 'MoI1GavgRes_MS_MS', 'maxG_MS', 'MoI0RstdTail_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI0Rtd_MS_MS', 'MoI0RavgRes_MS_MS', 'MoI0Rscale_GLS_GLS', 'scaleHead_MS', 'CoM1Rscale_MS_MS', 'CoM0YavgRes_GLS_GLS', 'scaleLength_MS', 'CoM0YstdTail_GLS_GLS', 'MoI1Rtd_MS_MS', 'MoI1Ytd_GLS_GLS', 'CoM1GstdTail_MS_MS'], -2.3484126984126985), (9, ['CoM0RstdTail_GLS_GLS', 'mG_MS', 'tMove_MS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'scaleLength_MS', 'MoI0Yscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'sG_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'mRmG_MS', 'MoI1Rscale_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'mR_MS', 'CoM1YavgRes_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'Rtail_GLS', 'CoM0Rscale_GLS_GLS', 'tScale_GLS', 'CoM1Ytd_GLS_GLS', 'Ytail_GLS', 'Gtail_GLS', 'MoI0GavgRes_MS_MS', 'MoI0Rtd_MS_MS', 'tScale_MS', 'CoM0YavgRes_GLS_GLS', 'mRmG_GLS', 'CoM1RavgRes_GLS_GLS', 'aY_GLS', 'tMove_GLS', 'MoI1GavgRes_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI1RstdTail_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'CoM1Gtd_MS_MS', 'CoM1Rscale_MS_MS', 'mYmG_GLS', 'CoM0RstdTail_MS_MS', 'MoI0Gtd_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'aR_GLS', 'CoM0GstdTail_MS_MS'], -2.5211018990929706), (10, ['CoM0RstdTail_GLS_GLS', 'mG_MS', 'tMove_MS', 'MoI0GavgRes_MS_MS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'CoM0Yscale_GLS_GLS', 'mRmG_MS', 'sG_GLS', 'scaleLength_MS', 'MoI1Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'Gtail_GLS', 'CoM1Rtd_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'mR_MS', 'MoI1GavgRes_GLS_GLS', 'tScale_MS', 'MoI0RstdTail_GLS_GLS', 'Ytail_GLS', 'CoM1Rtd_MS_MS', 'MoI0Rscale_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'tScale_GLS', 'mRmG_GLS', 'CoM1Gtd_MS_MS', 'aY_GLS', 'MoI0Gtd_MS_MS', 'CoM0YstdTail_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI1YstdTail_GLS_GLS', 'MoI0RstdTail_MS_MS', 'fracY_GLS', 'MoI1RstdTail_GLS_GLS', 'maxG_GLS', 'CoM0Rtd_MS_MS', 'MoI0Gtd_GLS_GLS', 'CoM0Gscale_MS_MS', 'maxR_GLS', 'CoM0RstdTail_MS_MS', 'CoM0YavgRes_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'MoI0RavgRes_MS_MS', 'Rtail_GLS', 'CoM0GstdTail_MS_MS', 'MoI1Ytd_GLS_GLS', 'MoI1Rtd_MS_MS', 'MoI1Yscale_GLS_GLS', 'CoM0GavgRes_MS_MS', 'MoI1Rscale_MS_MS', 'MoI0GstdTail_MS_MS', 'sR_GLS', 'maxG_MS', 'tailLength_MS', 'MoI0GavgRes_GLS_GLS'], -2.4053125), (11, ['CoM0RstdTail_GLS_GLS', 'mG_MS', 'scaleHead_MS', 'CoM0Yscale_GLS_GLS', 'MoI1GavgRes_MS_MS', 'tMove_MS', 'MoI1Rscale_GLS_GLS', 'scaleLength_MS', 'sG_GLS', 'MoI0GavgRes_MS_MS', 'MoI0Yscale_GLS_GLS', 'mRmG_MS', 'MoI0Gscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'tScale_GLS', 'mR_MS', 'Gtail_GLS', 'CoM0Rscale_GLS_GLS', 'maxR_GLS', 'Rtail_GLS', 'MoI1YstdTail_GLS_GLS', 'mRmG_GLS', 'aY_GLS', 'Ytail_GLS', 'CoM1Ytd_GLS_GLS', 'MoI0Rscale_MS_MS', 'MoI0RstdTail_GLS_GLS', 'CoM1Gtd_MS_MS', 'MoI0Rtd_MS_MS', 'MoI0Gtd_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'tScale_MS', 'CoM0YstdTail_GLS_GLS', 'tMove_GLS', 'CoM1Rscale_MS_MS', 'CoM1Rtd_GLS_GLS', 'mYmG_GLS', 'CoM0Ytd_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'maxG_MS', 'MoI1Rscale_MS_MS', 'MoI1Gtd_GLS_GLS', 'MoI1Gtd_MS_MS', 'sR_GLS', 'CoM0YavgRes_GLS_GLS', 'CoM1GstdTail_MS_MS', 'MoI0GstdTail_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'MoI1Yscale_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI0GstdTail_MS_MS', 'maxHead_MS', 'sG_MS', 'fracY_GLS', 'tailLength_MS', 'MoI1GavgRes_GLS_GLS', 'aR_GLS', 'MoI1Ytd_GLS_GLS', 'MoI0RavgRes_GLS_GLS', 'MoI1Rtd_MS_MS'], -2.1233177437641721), (12, ['CoM0RstdTail_GLS_GLS', 'mG_MS', 'tMove_MS', 'scaleHead_MS', 'CoM0Yscale_GLS_GLS', 'scaleLength_MS', 'MoI1GavgRes_MS_MS', 'MoI1Rscale_GLS_GLS', 'tScale_MS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'Rtail_GLS', 'maxG_MS', 'MoI0GavgRes_MS_MS', 'CoM0Ytd_GLS_GLS', 'mR_MS', 'Ytail_GLS', 'tScale_GLS', 'CoM1Ytd_GLS_GLS', 'CoM1Gtd_MS_MS', 'MoI0Gtd_MS_MS', 'MoI1RstdTail_GLS_GLS', 'aY_GLS', 'mRmG_MS', 'sG_GLS', 'mRmG_GLS', 'CoM1RavgRes_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'Gtail_GLS', 'maxR_GLS', 'CoM0RavgRes_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI0YstdTail_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI0Rtd_GLS_GLS', 'MoI0GavgRes_GLS_GLS', 'MoI0Rtd_MS_MS', 'MoI1Rtd_MS_MS'], -2.7836309523809524), (13, ['mG_MS', 'tMove_MS', 'CoM0RstdTail_GLS_GLS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'scaleLength_MS', 'MoI1Rscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'tScale_MS', 'MoI1GavgRes_GLS_GLS', 'mR_MS', 'CoM0Rscale_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM1Gtd_MS_MS', 'sG_GLS', 'mRmG_MS', 'mRmG_GLS', 'aY_GLS', 'tScale_GLS', 'MoI1RstdTail_GLS_GLS', 'CoM1Rscale_MS_MS', 'CoM0Yscale_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'MoI1Rscale_MS_MS', 'CoM1Ytd_GLS_GLS', 'fracY_GLS', 'CoM0RavgRes_GLS_GLS', 'MoI0Rtd_MS_MS', 'MoI0RstdTail_GLS_GLS', 'tMove_GLS', 'MoI0Gtd_MS_MS', 'CoM0YstdTail_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI0GstdTail_MS_MS', 'MoI0Yscale_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'CoM0Rtd_GLS_GLS', 'devHead_MS', 'MoI0Rscale_MS_MS', 'MoI0GstdTail_GLS_GLS', 'Ytail_GLS', 'CoM0GstdTail_MS_MS', 'CoM0Rtd_MS_MS', 'MoI0Gtd_GLS_GLS'], -2.7552579365079364), (14, ['CoM0RstdTail_GLS_GLS', 'tMove_MS', 'scaleHead_MS', 'CoM0Yscale_GLS_GLS', 'mG_MS', 'MoI1GavgRes_MS_MS', 'sG_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'tScale_MS', 'scaleLength_MS', 'CoM0Rscale_GLS_GLS', 'mRmG_GLS', 'Rtail_GLS', 'Gtail_GLS', 'MoI1YstdTail_GLS_GLS', 'mR_MS', 'tScale_GLS', 'Ytail_GLS', 'aR_GLS', 'MoI1GavgRes_GLS_GLS', 'aY_GLS', 'MoI0GavgRes_MS_MS', 'CoM1RavgRes_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'CoM0RstdTail_MS_MS', 'CoM0YavgRes_GLS_GLS', 'CoM0GavgRes_MS_MS', 'maxG_MS', 'MoI0Rtd_MS_MS', 'fracG_GLS', 'MoI1RavgRes_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'CoM1Gtd_MS_MS', 'mYmG_GLS', 'MoI0GstdTail_GLS_GLS', 'MoI0RavgRes_MS_MS', 'MoI0Gtd_GLS_GLS', 'MoI0Gtd_MS_MS', 'devHead_MS', 'MoI0GstdTail_MS_MS', 'CoM0YstdTail_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'mRmG_MS', 'CoM0Gscale_MS_MS', 'MoI1Ytd_GLS_GLS', 'MoI0RstdTail_MS_MS', 'tailHead_MS'], -2.4840501385739477), (15, ['mG_MS', 'scaleHead_MS', 'CoM0RstdTail_GLS_GLS', 'scaleLength_MS', 'CoM0Yscale_GLS_GLS', 'sG_GLS', 'tMove_MS', 'mRmG_MS', 'MoI0GavgRes_MS_MS', 'MoI1GavgRes_MS_MS', 'MoI1Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'Rtail_GLS', 'CoM1Rtd_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'MoI0RstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'mR_MS', 'MoI1YstdTail_GLS_GLS', 'maxR_GLS', 'tScale_GLS', 'CoM1Gtd_MS_MS', 'Gtail_GLS', 'mYmG_GLS', 'mRmG_GLS', 'MoI1GavgRes_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI0YstdTail_GLS_GLS', 'Ytail_GLS', 'MoI0GstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'maxG_MS', 'MoI1Gtd_MS_MS', 'CoM0YstdTail_GLS_GLS', 'aY_GLS', 'CoM0Ytd_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'MoI1Gscale_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'CoM0RstdTail_MS_MS', 'MoI0Gtd_MS_MS', 'tScale_MS', 'CoM0RavgRes_MS_MS', 'CoM0RavgRes_GLS_GLS', 'sR_GLS', 'MoI0Rtd_MS_MS', 'CoM1GstdTail_MS_MS', 'CoM1Rtd_MS_MS', 'CoM1YavgRes_GLS_GLS', 'MoI0Rscale_MS_MS', 'aR_GLS', 'CoM1Rscale_MS_MS', 'aG_MS', 'MoI0RavgRes_GLS_GLS'], -2.1241843033509697), (16, ['CoM0RstdTail_GLS_GLS', 'mG_MS', 'tMove_MS', 'scaleHead_MS', 'MoI0GavgRes_MS_MS', 'CoM0Yscale_GLS_GLS', 'mRmG_MS', 'MoI1GavgRes_MS_MS', 'sG_GLS', 'MoI1Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'scaleLength_MS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'tScale_MS', 'mR_MS', 'CoM1Rtd_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'Rtail_GLS', 'CoM1Ytd_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'tScale_GLS', 'Ytail_GLS', 'mRmG_GLS', 'maxR_GLS', 'MoI0RstdTail_GLS_GLS', 'aY_GLS', 'Gtail_GLS', 'CoM0RavgRes_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'maxG_MS', 'MoI0Rtd_MS_MS', 'CoM1Gtd_MS_MS', 'tMove_GLS', 'fracY_GLS', 'MoI0YstdTail_GLS_GLS', 'devHead_MS', 'CoM0YavgRes_GLS_GLS', 'MoI0Rscale_MS_MS', 'CoM0YstdTail_GLS_GLS', 'mYmG_GLS', 'CoM1Rscale_MS_MS', 'MoI0Gtd_GLS_GLS', 'CoM0GavgRes_MS_MS', 'MoI0Gtd_MS_MS', 'MoI0RavgRes_MS_MS', 'CoM0RstdTail_MS_MS', 'CoM0Rtd_MS_MS', 'MoI1Ytd_GLS_GLS', 'CoM1Yscale_GLS_GLS', 'MoI1Gtd_MS_MS', 'CoM1GstdTail_MS_MS', 'MoI0Rscale_GLS_GLS'], -2.6346082136558326)]


    # '''this was from the 3/4/18 group opt by rank run, using training set (70% training, with at least 3 members to train) for RNAi list'''
    # opt_res = [(None, ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'mRmG_MS', 'MoI1Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'tMove_MS', 'CoM1Gtd_MS_MS', 'CoM0Yscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'mR_MS', 'sG_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'devHead_MS', 'CoM0Ytd_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM0Rtd_MS_MS', 'tMove_GLS', 'MoI1YstdTail_GLS_GLS', 'mRmG_GLS', 'aY_GLS', 'MoI1GavgRes_GLS_GLS', 'CoM0GavgRes_MS_MS', 'tScale_GLS', 'tScale_MS', 'CoM1Ytd_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'scaleLength_MS', 'CoM0YavgRes_GLS_GLS', 'Gtail_GLS', 'MoI0Gtd_MS_MS', 'CoM1Rscale_MS_MS', 'MoI0RstdTail_GLS_GLS', 'MoI0Rscale_MS_MS', 'mYmG_GLS', 'Rtail_GLS', 'tailHead_MS', 'maxR_GLS', 'MoI0Rtd_MS_MS', 'CoM0RavgRes_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'Ytail_GLS', 'aR_GLS', 'MoI0Rscale_GLS_GLS', 'MoI1Gtd_MS_MS', 'sR_GLS', 'MoI0RstdTail_MS_MS', 'CoM1Yscale_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'CoM1Gscale_MS_MS', 'MoI0GstdTail_GLS_GLS', 'MoI1RavgRes_GLS_GLS', 'maxG_MS', 'CoM0RstdTail_MS_MS', 'MoI0GstdTail_MS_MS', 'MoI0Gtd_GLS_GLS', 'CoM1GstdTail_MS_MS', 'MoI0RavgRes_MS_MS', 'CoM0GstdTail_MS_MS', 'CoM1Rtd_MS_MS', 'MoI1Ytd_GLS_GLS'], -2.4274514203829685), (0, ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'CoM1Gtd_MS_MS', 'tMove_MS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'mRmG_MS', 'MoI1Rscale_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM1RstdTail_GLS_GLS', 'tailHead_MS', 'CoM1Rtd_MS_MS', 'MoI0Yscale_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'tScale_MS', 'sG_GLS', 'tScale_GLS', 'scaleLength_MS', 'MoI0RavgRes_MS_MS', 'CoM1Rtd_GLS_GLS', 'mRmG_GLS', 'MoI1GavgRes_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'tMove_GLS', 'MoI0Rtd_MS_MS', 'CoM1Rscale_MS_MS', 'aY_GLS', 'Rtail_GLS', 'CoM0YavgRes_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'MoI0Gtd_MS_MS', 'mR_MS', 'fracY_GLS', 'MoI1RstdTail_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'Gtail_GLS', 'MoI0Rscale_GLS_GLS', 'Ytail_GLS', 'aR_GLS', 'MoI0Rscale_MS_MS', 'devHead_MS', 'CoM0YstdTail_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'mYmG_GLS', 'aG_MS', 'MoI1Rtd_MS_MS', 'MoI0RstdTail_MS_MS', 'sR_GLS', 'MoI1Gtd_MS_MS', 'MoI0RstdTail_GLS_GLS'], -2.5463932980599644), (1, ['mG_MS', 'tailHead_MS', 'CoM0RstdTail_GLS_GLS', 'tMove_MS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'mRmG_MS', 'MoI0Yscale_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'CoM1Gtd_MS_MS', 'MoI1Rscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM0Rscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM0GavgRes_MS_MS', 'aY_GLS', 'sG_GLS', 'CoM0Ytd_GLS_GLS', 'mRmG_GLS', 'devHead_MS', 'scaleLength_MS', 'tScale_MS', 'MoI0Rtd_MS_MS', 'sR_GLS', 'MoI1YstdTail_GLS_GLS', 'mYmG_GLS', 'tScale_GLS', 'mR_MS', 'MoI0Gscale_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI0Gtd_MS_MS', 'CoM1Ytd_GLS_GLS', 'MoI0Rscale_MS_MS', 'MoI0Gtd_GLS_GLS', 'MoI1GstdTail_GLS_GLS', 'Rtail_GLS', 'Gtail_GLS', 'Ytail_GLS', 'MoI0Rscale_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'MoI1Gtd_MS_MS', 'CoM0Rtd_MS_MS', 'CoM0RavgRes_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'fracY_GLS', 'aR_GLS', 'MoI1RavgRes_GLS_GLS', 'CoM0GstdTail_MS_MS', 'MoI1Rtd_MS_MS', 'CoM1GstdTail_MS_MS', 'MoI0Ytd_GLS_GLS'], -2.3955769715293522), (2, ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'CoM1Gtd_MS_MS', 'mRmG_MS', 'tMove_MS', 'MoI1Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM0Rscale_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'devHead_MS', 'sG_GLS', 'tScale_GLS', 'CoM1YstdTail_GLS_GLS', 'mRmG_GLS', 'scaleLength_MS', 'aY_GLS', 'MoI0Gscale_GLS_GLS', 'mYmG_GLS', 'CoM1Rtd_MS_MS', 'tScale_MS', 'CoM0Rtd_MS_MS', 'CoM1Ytd_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'MoI0Rtd_MS_MS', 'Gtail_GLS', 'aR_GLS', 'CoM1Rscale_MS_MS', 'CoM1Rtd_GLS_GLS', 'Rtail_GLS', 'MoI1Rtd_MS_MS', 'tMove_GLS', 'CoM0RavgRes_GLS_GLS', 'Ytail_GLS', 'sR_GLS', 'CoM1GstdTail_MS_MS', 'mR_MS', 'MoI0Ytd_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI1RstdTail_GLS_GLS', 'maxG_MS', 'CoM1YavgRes_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'MoI1RavgRes_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'MoI0RavgRes_MS_MS', 'MoI1Ytd_GLS_GLS', 'fracY_GLS'], -2.6041087175610986), (3, ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'tMove_MS', 'mRmG_MS', 'MoI0Yscale_GLS_GLS', 'CoM1Gtd_MS_MS', 'CoM1RstdTail_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'mR_MS', 'CoM1YstdTail_GLS_GLS', 'tScale_MS', 'CoM1Rtd_MS_MS', 'MoI1GavgRes_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'CoM0GavgRes_MS_MS', 'tScale_GLS', 'sG_GLS', 'scaleLength_MS', 'MoI0RavgRes_MS_MS', 'CoM1Ytd_GLS_GLS', 'CoM0Rtd_MS_MS', 'sR_GLS', 'Ytail_GLS', 'MoI0RavgRes_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI0Ytd_GLS_GLS', 'tMove_GLS', 'CoM0Rscale_MS_MS', 'CoM0YstdTail_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'CoM1Gscale_MS_MS', 'MoI0Gscale_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'mYmG_GLS', 'MoI1Gtd_MS_MS', 'MoI1Yscale_GLS_GLS', 'MoI1Rscale_MS_MS', 'CoM0Gscale_MS_MS', 'CoM1Rtd_GLS_GLS', 'devHead_MS', 'MoI0RstdTail_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'maxR_GLS', 'CoM0RavgRes_MS_MS', 'MoI1Ytd_GLS_GLS', 'CoM0GstdTail_MS_MS', 'MoI0Rscale_MS_MS', 'MoI0GavgRes_GLS_GLS', 'Rtail_GLS', 'sG_MS', 'tailLength_MS', 'maxG_MS', 'MoI0GstdTail_MS_MS', 'MoI1Rtd_MS_MS', 'MoI0Rtd_MS_MS'], -2.2000296044343663), (4, ['mG_MS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'CoM0RstdTail_GLS_GLS', 'tailHead_MS', 'mRmG_MS', 'tMove_MS', 'MoI0Yscale_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'CoM1Gtd_MS_MS', 'MoI0GavgRes_MS_MS', 'CoM1RstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'sG_GLS', 'tScale_GLS', 'devHead_MS', 'mR_MS', 'MoI1YstdTail_GLS_GLS', 'mRmG_GLS', 'aY_GLS', 'mYmG_GLS', 'aR_GLS', 'MoI1RstdTail_GLS_GLS', 'CoM1Yscale_GLS_GLS', 'tScale_MS', 'Gtail_GLS', 'CoM1Rtd_GLS_GLS', 'CoM0Rtd_MS_MS', 'CoM1RavgRes_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'MoI0YstdTail_GLS_GLS', 'MoI0Rtd_MS_MS'], -3.0630637440161252), (5, ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'mRmG_MS', 'MoI1Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'tMove_MS', 'CoM0Yscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'tailHead_MS', 'CoM1Gtd_MS_MS', 'CoM1YstdTail_GLS_GLS', 'sG_GLS', 'CoM0Rscale_GLS_GLS', 'scaleLength_MS', 'tScale_MS', 'mRmG_GLS', 'CoM1Rtd_GLS_GLS', 'CoM0Rtd_MS_MS', 'tScale_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'Gtail_GLS', 'CoM1Gscale_MS_MS', 'MoI0Rtd_MS_MS', 'mR_MS', 'devHead_MS', 'CoM0RavgRes_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI0Rscale_GLS_GLS', 'Ytail_GLS', 'MoI1YstdTail_GLS_GLS', 'mYmG_GLS', 'MoI0Gtd_MS_MS', 'MoI0Ytd_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'Rtail_GLS', 'MoI0YstdTail_GLS_GLS', 'aY_GLS', 'CoM1RavgRes_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'maxG_MS', 'CoM0YstdTail_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'sR_GLS', 'MoI1Ytd_GLS_GLS', 'MoI1Gtd_MS_MS', 'MoI1Gtd_GLS_GLS', 'CoM1GstdTail_MS_MS', 'CoM0GavgRes_MS_MS', 'MoI1Rtd_MS_MS'], -2.4300885770975054), (6, ['CoM0RstdTail_GLS_GLS', 'CoM1Gtd_MS_MS', 'mRmG_MS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'tMove_MS', 'MoI0Yscale_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'mG_MS', 'scaleLength_MS', 'CoM0Rscale_GLS_GLS', 'tScale_GLS', 'CoM1Rtd_MS_MS', 'CoM1RstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'devHead_MS', 'MoI0Rtd_MS_MS', 'mR_MS', 'CoM0YstdTail_GLS_GLS', 'MoI0GavgRes_MS_MS', 'sG_GLS', 'MoI1Rscale_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'mRmG_GLS', 'aY_GLS', 'Ytail_GLS', 'CoM1Ytd_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI0Gtd_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'CoM0Rtd_MS_MS', 'MoI0YstdTail_GLS_GLS', 'Rtail_GLS', 'tScale_MS', 'CoM1RavgRes_GLS_GLS', 'Gtail_GLS', 'mYmG_GLS', 'aR_GLS', 'CoM0RavgRes_GLS_GLS', 'maxR_GLS', 'maxG_MS', 'tailHead_MS', 'MoI0GstdTail_GLS_GLS', 'MoI1Gtd_MS_MS', 'MoI0GstdTail_MS_MS'], -2.6499964569160994), (7, ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'tMove_MS', 'mRmG_MS', 'tailHead_MS', 'CoM0Yscale_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'CoM1Gtd_MS_MS', 'CoM1RstdTail_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'sG_GLS', 'mR_MS', 'MoI0Yscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM0Rscale_GLS_GLS', 'CoM1Rtd_MS_MS', 'tScale_GLS', 'mRmG_GLS', 'MoI0Rtd_MS_MS', 'MoI1Gtd_MS_MS', 'MoI0Gtd_GLS_GLS', 'scaleLength_MS', 'tScale_MS', 'MoI0Rscale_MS_MS', 'CoM1YstdTail_GLS_GLS', 'CoM0Rtd_MS_MS', 'aY_GLS', 'mYmG_GLS', 'MoI1YstdTail_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'Gtail_GLS', 'Rtail_GLS', 'MoI1Ytd_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'Ytail_GLS', 'fracY_GLS', 'CoM1GstdTail_MS_MS'], -2.4203599773242628), (8, ['mG_MS', 'CoM1Gtd_MS_MS', 'tMove_MS', 'MoI1Rscale_GLS_GLS', 'mRmG_MS', 'MoI0Yscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'tScale_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM0RstdTail_GLS_GLS', 'MoI0Rtd_MS_MS', 'mR_MS', 'sG_GLS', 'tailHead_MS', 'MoI0GavgRes_MS_MS', 'CoM0Rtd_MS_MS', 'MoI0Rscale_MS_MS', 'tScale_MS', 'mYmG_GLS', 'Rtail_GLS', 'Ytail_GLS', 'mRmG_GLS', 'aY_GLS', 'scaleHead_MS', 'tailLength_MS', 'MoI1Rscale_MS_MS', 'CoM1Rtd_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'MoI0Gtd_MS_MS', 'CoM0Ytd_GLS_GLS', 'CoM1Yscale_GLS_GLS', 'fracY_GLS', 'CoM1Rscale_MS_MS', 'CoM0RavgRes_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'Gtail_GLS', 'devHead_MS', 'MoI1Gtd_MS_MS', 'MoI0Rscale_GLS_GLS', 'sR_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI0Rtd_GLS_GLS', 'aR_GLS', 'MoI0RavgRes_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'MoI0YstdTail_GLS_GLS', 'MoI1Rtd_MS_MS'], -2.3896818310657597), (9, ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'tMove_MS', 'mRmG_MS', 'CoM1Gtd_MS_MS', 'MoI1Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM0Yscale_GLS_GLS', 'mR_MS', 'sG_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI1Rscale_MS_MS', 'tScale_GLS', 'MoI0Rtd_MS_MS', 'mRmG_GLS', 'CoM1YstdTail_GLS_GLS', 'tScale_MS', 'CoM1YavgRes_GLS_GLS', 'CoM0Rtd_MS_MS', 'MoI1YstdTail_GLS_GLS', 'MoI0Gtd_MS_MS', 'CoM0GstdTail_MS_MS', 'CoM0Gtd_MS_MS', 'CoM1Ytd_GLS_GLS', 'scaleLength_MS', 'MoI0RstdTail_GLS_GLS', 'MoI0RavgRes_MS_MS', 'CoM0YavgRes_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'aY_GLS', 'fracY_GLS', 'MoI0RavgRes_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI0Gscale_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'tailHead_MS', 'Gtail_GLS', 'tMove_GLS', 'CoM1Rtd_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'devHead_MS', 'CoM1RavgRes_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'maxR_GLS', 'MoI0YstdTail_GLS_GLS', 'MoI1Rtd_MS_MS', 'Rtail_GLS'], -2.5140086451247168), (10, ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'mRmG_MS', 'MoI1Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'tMove_MS', 'CoM0Yscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM1Gtd_MS_MS', 'mR_MS', 'sG_GLS', 'MoI1Rscale_MS_MS', 'CoM0Rtd_MS_MS', 'CoM1RstdTail_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'MoI0Rtd_MS_MS', 'CoM1YstdTail_GLS_GLS', 'MoI0RstdTail_GLS_GLS', 'MoI1Gtd_MS_MS', 'tScale_GLS', 'tScale_MS', 'CoM1YavgRes_GLS_GLS', 'mRmG_GLS', 'devHead_MS', 'aY_GLS', 'CoM1Ytd_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI1GavgRes_GLS_GLS', 'scaleLength_MS', 'MoI0Rscale_GLS_GLS', 'tMove_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI1YstdTail_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'Ytail_GLS', 'tailHead_MS', 'CoM0YstdTail_GLS_GLS', 'maxR_GLS', 'CoM1RavgRes_GLS_GLS', 'Rtail_GLS', 'mYmG_GLS', 'MoI1RstdTail_GLS_GLS', 'sR_GLS', 'aG_MS', 'Gtail_GLS', 'MoI0RstdTail_MS_MS', 'MoI0YstdTail_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'MoI0Rscale_MS_MS', 'CoM1Rscale_GLS_GLS', 'CoM0GstdTail_MS_MS', 'MoI1Yscale_GLS_GLS', 'fracY_GLS', 'MoI0GstdTail_MS_MS', 'MoI1RavgRes_GLS_GLS', 'MoI0RavgRes_GLS_GLS'], -2.3932291666666665), (11, ['mG_MS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'tMove_MS', 'CoM0RstdTail_GLS_GLS', 'mRmG_MS', 'CoM1Gtd_MS_MS', 'MoI1Rscale_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM1RstdTail_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'mR_MS', 'tMove_GLS', 'MoI0Rtd_MS_MS', 'sG_GLS', 'CoM1YstdTail_GLS_GLS', 'tScale_MS', 'CoM0Rscale_GLS_GLS', 'tScale_GLS', 'Ytail_GLS', 'mRmG_GLS', 'MoI0Gtd_MS_MS', 'MoI0RstdTail_GLS_GLS', 'tailHead_MS', 'aY_GLS', 'scaleLength_MS', 'MoI1Yscale_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI1Gtd_MS_MS', 'MoI0Rscale_MS_MS', 'mYmG_GLS', 'MoI0Gtd_GLS_GLS', 'maxR_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'Rtail_GLS', 'Gtail_GLS', 'CoM1Ytd_GLS_GLS', 'tailLength_MS', 'sR_GLS', 'MoI1Ytd_GLS_GLS', 'MoI0RavgRes_MS_MS', 'MoI0YstdTail_GLS_GLS', 'MoI1Rscale_MS_MS', 'MoI0GstdTail_GLS_GLS', 'devHead_MS', 'fracY_GLS', 'CoM0Ytd_GLS_GLS'], -2.3559233276643989), (12, ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'mRmG_MS', 'MoI1Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'tMove_MS', 'CoM1Gtd_MS_MS', 'MoI0GavgRes_MS_MS', 'CoM0Yscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'mR_MS', 'CoM1YstdTail_GLS_GLS', 'tScale_MS', 'CoM1Rtd_MS_MS', 'MoI1GavgRes_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'sG_GLS', 'tScale_GLS', 'scaleLength_MS', 'CoM1Rtd_GLS_GLS', 'devHead_MS', 'MoI0RavgRes_MS_MS', 'CoM0Gscale_MS_MS', 'MoI0Gtd_GLS_GLS', 'CoM0Rtd_MS_MS', 'Rtail_GLS', 'mRmG_GLS', 'aY_GLS', 'tMove_GLS', 'MoI1RstdTail_GLS_GLS', 'sR_GLS', 'CoM1Ytd_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'Ytail_GLS', 'MoI0RstdTail_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'MoI0Rtd_MS_MS', 'CoM0YavgRes_GLS_GLS', 'maxG_MS', 'MoI1RavgRes_GLS_GLS', 'MoI0Gtd_MS_MS'], -2.8691836734693874), (13, ['mG_MS', 'tailHead_MS', 'mRmG_MS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'MoI1Rscale_GLS_GLS', 'CoM0RstdTail_GLS_GLS', 'tMove_MS', 'CoM1Gtd_MS_MS', 'CoM0Yscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'MoI0Yscale_GLS_GLS', 'CoM0Rtd_MS_MS', 'mR_MS', 'maxR_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'sG_GLS', 'MoI1YstdTail_GLS_GLS', 'mRmG_GLS', 'aY_GLS', 'MoI1RstdTail_GLS_GLS', 'fracY_GLS', 'CoM1YstdTail_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'tScale_MS', 'tScale_GLS', 'MoI0Gscale_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'scaleLength_MS', 'CoM0YavgRes_GLS_GLS', 'devHead_MS', 'MoI1GstdTail_GLS_GLS', 'MoI0Rtd_MS_MS', 'MoI0Gtd_MS_MS', 'MoI0Gtd_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'MoI0RstdTail_GLS_GLS', 'CoM1Rscale_MS_MS', 'CoM1Yscale_GLS_GLS', 'CoM1GstdTail_MS_MS', 'Gtail_GLS', 'CoM0RstdTail_MS_MS', 'mYmG_GLS', 'Ytail_GLS', 'CoM1Ytd_GLS_GLS', 'MoI1Rtd_MS_MS', 'maxG_MS', 'CoM1Rtd_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'MoI1Gtd_MS_MS', 'sR_GLS', 'CoM0YstdTail_GLS_GLS', 'aR_GLS', 'tMove_GLS', 'MoI1RavgRes_GLS_GLS', 'MoI0RstdTail_MS_MS', 'MoI0YstdTail_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'MoI0Rscale_MS_MS', 'MoI1Rscale_MS_MS', 'MoI0Rtd_GLS_GLS', 'MoI0GstdTail_MS_MS', 'MoI0RavgRes_MS_MS', 'MoI1Yscale_GLS_GLS', 'tailLength_MS', 'sG_MS', 'MoI0RavgRes_GLS_GLS', 'CoM1GavgRes_MS_MS'], -2.353094529478458), (14, ['tailHead_MS', 'CoM0RstdTail_GLS_GLS', 'mG_MS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'tMove_MS', 'tScale_MS', 'MoI1Rscale_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'CoM1Gtd_MS_MS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'devHead_MS', 'CoM1Ytd_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'mR_MS', 'sG_GLS', 'MoI0Rtd_MS_MS', 'tScale_GLS', 'mRmG_GLS', 'aY_GLS', 'fracY_GLS', 'CoM1Rtd_GLS_GLS', 'scaleLength_MS', 'MoI0GavgRes_MS_MS', 'Ytail_GLS', 'CoM1Rscale_MS_MS', 'MoI0Gtd_MS_MS', 'CoM1Rtd_MS_MS', 'Rtail_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'Gtail_GLS', 'CoM0GavgRes_MS_MS', 'MoI0GstdTail_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'mRmG_MS', 'CoM0RavgRes_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'mYmG_GLS', 'MoI1Gtd_MS_MS', 'CoM0GstdTail_MS_MS', 'CoM1YavgRes_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'MoI0RstdTail_GLS_GLS', 'MoI1Rscale_MS_MS', 'aG_MS', 'MoI1Ytd_GLS_GLS', 'MoI0RstdTail_MS_MS', 'CoM0RstdTail_MS_MS', 'MoI0Gtd_GLS_GLS'], -2.4675453514739232), (15, ['mG_MS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'mRmG_MS', 'tailHead_MS', 'MoI1Rscale_GLS_GLS', 'CoM0RstdTail_GLS_GLS', 'MoI0GavgRes_MS_MS', 'MoI0Yscale_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'CoM1Gtd_MS_MS', 'MoI1GavgRes_GLS_GLS', 'devHead_MS', 'sG_GLS', 'tMove_MS', 'mR_MS', 'CoM0Ytd_GLS_GLS', 'tScale_GLS', 'CoM1Ytd_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'scaleLength_MS', 'MoI0Gscale_GLS_GLS', 'MoI0RavgRes_MS_MS', 'MoI1Gtd_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI1Rtd_MS_MS', 'MoI0Rtd_MS_MS', 'MoI0YstdTail_GLS_GLS', 'mYmG_GLS', 'MoI0Gtd_GLS_GLS', 'mRmG_GLS', 'CoM0Rtd_MS_MS', 'Gtail_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'maxR_GLS', 'CoM0Rtd_GLS_GLS', 'CoM1GstdTail_MS_MS', 'MoI0GstdTail_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'Rtail_GLS', 'aY_GLS', 'MoI1Gscale_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'Ytail_GLS', 'CoM1RavgRes_GLS_GLS', 'tScale_MS', 'MoI1Gtd_MS_MS', 'CoM1Rscale_MS_MS', 'MoI0Gtd_MS_MS', 'CoM0YavgRes_GLS_GLS', 'maxG_MS', 'MoI1Ytd_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'MoI0RstdTail_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'MoI0GstdTail_MS_MS', 'aR_GLS', 'CoM0GstdTail_MS_MS', 'MoI1RstdTail_MS_MS', 'sG_MS', 'MoI1Yscale_GLS_GLS', 'tMove_GLS'], -2.0837515747039559), (16, ['mG_MS', 'tailHead_MS', 'mRmG_MS', 'CoM0RstdTail_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI1GavgRes_MS_MS', 'scaleHead_MS', 'CoM1Gtd_MS_MS', 'CoM0Yscale_GLS_GLS', 'tMove_MS', 'MoI0GavgRes_MS_MS', 'MoI0Yscale_GLS_GLS', 'devHead_MS', 'CoM1YstdTail_GLS_GLS', 'tScale_GLS', 'CoM0Rscale_GLS_GLS', 'mR_MS', 'sG_GLS', 'MoI1YstdTail_GLS_GLS', 'scaleLength_MS', 'tScale_MS', 'MoI1GavgRes_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'tMove_GLS', 'Gtail_GLS', 'CoM1Gscale_MS_MS', 'MoI0Gtd_GLS_GLS', 'Rtail_GLS', 'MoI0RstdTail_GLS_GLS', 'MoI0Rtd_MS_MS', 'CoM1Rscale_MS_MS', 'mYmG_GLS', 'aY_GLS', 'mRmG_GLS', 'fracY_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI0Gtd_MS_MS', 'aG_MS', 'MoI1Rtd_MS_MS', 'CoM1RavgRes_GLS_GLS', 'CoM0Rtd_MS_MS', 'MoI1Gtd_MS_MS', 'MoI0GstdTail_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'CoM1GstdTail_MS_MS', 'Ytail_GLS', 'CoM0YstdTail_GLS_GLS', 'MoI1YavgRes_GLS_GLS', 'MoI1Rscale_MS_MS', 'aR_GLS', 'MoI1Ytd_GLS_GLS', 'MoI0RavgRes_MS_MS', 'MoI1Yscale_GLS_GLS', 'MoI0Rscale_MS_MS', 'CoM0GstdTail_MS_MS'], -2.6013510959939534)]


    # '''this was from the 3/2/18 group opt by rank run, using training set (50%) for RNAi list'''
    # opt_res = [(0, ['mG_MS', 'mRmG_GLS', 'aY_GLS', 'tScale_MS', 'CoM0RstdTail_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'CoM1Gtd_MS_MS', 'Ytail_GLS', 'fracY_GLS', 'mYmG_GLS', 'mRmG_MS', 'CoM1Rscale_MS_MS', 'MoI1Yscale_GLS_GLS', 'tScale_GLS', 'MoI0GavgRes_MS_MS', 'fracR_GLS', 'MoI0Rtd_MS_MS', 'CoM1YavgRes_GLS_GLS', 'maxG_GLS', 'MoI1Rscale_MS_MS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'MoI0Gtd_MS_MS', 'Rtail_GLS', 'devLength_MS', 'CoM0YstdTail_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'CoM1RstdTail_MS_MS', 'CoM0Rscale_GLS_GLS', 'aR_GLS', 'MoI1GavgRes_MS_MS', 'scaleLength_MS', 'MoI1GstdTail_GLS_GLS', 'CoM0RstdTail_MS_MS', 'MoI1Gtd_MS_MS', 'CoM0Gtd_MS_MS', 'MoI0Rtd_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'tailLength_MS', 'MoI1Rtd_MS_MS', 'CoM1Ytd_GLS_GLS', 'MoI0YavgRes_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'Gtail_GLS', 'maxR_GLS', 'sG_GLS', 'MoI1GavgRes_GLS_GLS', 'MoI1YavgRes_GLS_GLS', 'aG_GLS', 'CoM1Rtd_MS_MS', 'mR_MS', 'MoI1Rtd_GLS_GLS', 'sR_GLS', 'MoI1Gtd_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'CoM0Rtd_GLS_GLS', 'tMove_MS', 'MoI1Ytd_GLS_GLS', 'maxHead_MS', 'scaleHead_MS', 'devHead_MS', 'MoI0YstdTail_GLS_GLS', 'CoM1GstdTail_MS_MS'], -1.9516534391534393), (1, ['mG_MS', 'mRmG_GLS', 'aY_GLS', 'CoM1RavgRes_GLS_GLS', 'CoM1Gtd_MS_MS', 'CoM0RstdTail_GLS_GLS', 'mYmG_GLS', 'fracY_GLS', 'CoM1YstdTail_GLS_GLS', 'Ytail_GLS', 'Gtail_GLS', 'tScale_GLS', 'mR_MS', 'mRmG_MS', 'MoI1Rscale_GLS_GLS', 'MoI1GavgRes_MS_MS', 'CoM0Rscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI0GavgRes_MS_MS', 'Rtail_GLS', 'CoM1Rscale_MS_MS', 'maxR_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'aR_GLS', 'MoI1Rscale_MS_MS', 'CoM0RstdTail_MS_MS', 'CoM0Ytd_GLS_GLS', 'sG_GLS', 'maxG_GLS', 'tMove_GLS', 'MoI0Gtd_MS_MS', 'MoI0Rtd_MS_MS', 'CoM1Gscale_MS_MS', 'tMove_MS', 'MoI0Ytd_GLS_GLS', 'fracR_GLS', 'MoI1Gtd_GLS_GLS', 'tailLength_MS', 'MoI1YavgRes_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'aG_GLS', 'MoI0Gscale_MS_MS', 'MoI1Gtd_MS_MS', 'MoI1Gscale_MS_MS', 'MoI1Rtd_GLS_GLS', 'MoI1RavgRes_GLS_GLS'], -2.1523236331569664), (2, ['mG_MS', 'mRmG_GLS', 'aY_GLS', 'tScale_MS', 'CoM0RstdTail_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'CoM1Gtd_MS_MS', 'Ytail_GLS', 'MoI0GavgRes_MS_MS', 'scaleLength_MS', 'Gtail_GLS', 'fracY_GLS', 'mYmG_GLS', 'tScale_GLS', 'aR_GLS', 'mR_MS', 'maxG_GLS', 'fracR_GLS', 'CoM0YstdTail_GLS_GLS', 'Rtail_GLS', 'MoI1YavgRes_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'scaleHead_MS', 'CoM1YavgRes_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'maxR_GLS', 'CoM0GstdTail_MS_MS', 'MoI1GavgRes_MS_MS', 'MoI0GstdTail_MS_MS', 'CoM1GavgRes_MS_MS', 'MoI1RstdTail_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI0Rtd_MS_MS'], -3.0361287477954146), (3, ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'mRmG_MS', 'MoI0GavgRes_MS_MS', 'MoI1GavgRes_MS_MS', 'CoM1RavgRes_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'CoM1RstdTail_MS_MS', 'MoI1Rscale_MS_MS', 'MoI1Gscale_GLS_GLS', 'maxR_GLS', 'tMove_GLS', 'MoI0Gtd_MS_MS', 'aY_GLS', 'Ytail_GLS', 'Rtail_GLS', 'CoM0Ytd_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'MoI0Rtd_MS_MS', 'MoI0GstdTail_GLS_GLS', 'CoM1Rtd_MS_MS', 'CoM1YavgRes_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'sG_GLS', 'aR_GLS', 'CoM1Rtd_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'tScale_GLS', 'MoI0Rtd_GLS_GLS', 'MoI1GstdTail_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'mYmG_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI1Rtd_GLS_GLS', 'MoI0GstdTail_MS_MS', 'CoM0Rscale_MS_MS', 'MoI0Rscale_GLS_GLS'], -2.5620855379188709), (4, ['mG_MS', 'mRmG_GLS', 'aY_GLS', 'CoM1RavgRes_GLS_GLS', 'CoM1Gtd_MS_MS', 'CoM0RstdTail_GLS_GLS', 'mYmG_GLS', 'fracY_GLS', 'Gtail_GLS', 'mR_MS', 'tScale_GLS', 'mRmG_MS', 'MoI1Yscale_GLS_GLS', 'Ytail_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'maxG_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'Rtail_GLS', 'maxR_GLS', 'CoM0Ytd_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'MoI1Rscale_MS_MS', 'aR_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI1YavgRes_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'MoI1GavgRes_MS_MS', 'CoM0RstdTail_MS_MS', 'MoI0RstdTail_MS_MS', 'MoI0Rtd_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'CoM1Rscale_MS_MS', 'sG_GLS', 'MoI1GstdTail_GLS_GLS', 'tailHead_MS', 'MoI1GavgRes_GLS_GLS', 'CoM1Rtd_MS_MS', 'CoM1Yscale_GLS_GLS', 'CoM1GstdTail_MS_MS', 'CoM0GstdTail_MS_MS', 'MoI0RavgRes_GLS_GLS', 'MoI1Rtd_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'MoI0Gtd_MS_MS'], -2.2552336860670192), (5, ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'mRmG_MS', 'CoM1RavgRes_GLS_GLS', 'tScale_MS', 'mRmG_GLS', 'aY_GLS', 'fracY_GLS', 'Ytail_GLS', 'mYmG_GLS', 'mR_MS', 'MoI1Rscale_GLS_GLS', 'Gtail_GLS', 'CoM1Gtd_MS_MS', 'tScale_GLS', 'MoI1Rscale_MS_MS', 'CoM1RstdTail_GLS_GLS', 'aR_GLS', 'maxG_GLS', 'Rtail_GLS', 'CoM1YavgRes_GLS_GLS', 'scaleLength_MS', 'MoI0GstdTail_GLS_GLS', 'maxR_GLS', 'CoM1Ytd_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI0Rtd_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI1GavgRes_MS_MS', 'MoI1YavgRes_GLS_GLS', 'tMove_MS', 'CoM0Rtd_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'sG_GLS', 'MoI0Gtd_MS_MS', 'MoI0Rtd_MS_MS', 'MoI1Gtd_MS_MS'], -2.4614583333333333), (6, ['CoM0RstdTail_GLS_GLS', 'devHead_MS', 'CoM1YstdTail_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'mRmG_GLS', 'aY_GLS', 'tScale_GLS', 'fracY_GLS', 'mYmG_GLS', 'MoI1GavgRes_MS_MS', 'CoM1Gtd_MS_MS', 'Ytail_GLS', 'tMove_MS', 'MoI1YstdTail_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'Rtail_GLS', 'CoM1YavgRes_GLS_GLS', 'MoI0GavgRes_MS_MS', 'mR_MS', 'CoM0YstdTail_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'sG_GLS', 'CoM0RstdTail_MS_MS', 'maxG_GLS', 'MoI0GstdTail_GLS_GLS', 'CoM1RstdTail_MS_MS', 'MoI0Rtd_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'MoI0GstdTail_MS_MS', 'MoI0Ytd_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI0Rtd_MS_MS', 'MoI0YavgRes_GLS_GLS'], -2.9000694444444441), (7, ['mG_MS', 'mRmG_GLS', 'aY_GLS', 'tScale_MS', 'CoM0RstdTail_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'CoM1Gtd_MS_MS', 'Ytail_GLS', 'fracY_GLS', 'mRmG_MS', 'mYmG_GLS', 'CoM1Rscale_MS_MS', 'MoI1Yscale_GLS_GLS', 'MoI1GavgRes_MS_MS', 'CoM1YstdTail_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'MoI1RavgRes_MS_MS', 'MoI1GstdTail_GLS_GLS', 'CoM1RstdTail_MS_MS', 'MoI0Gtd_MS_MS', 'tScale_GLS', 'Gtail_GLS', 'Rtail_GLS', 'CoM0RstdTail_MS_MS', 'CoM0Ytd_GLS_GLS', 'aR_GLS', 'sG_GLS', 'maxR_GLS', 'MoI1Rscale_MS_MS', 'MoI0GavgRes_MS_MS', 'CoM0YstdTail_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'aG_GLS', 'CoM1GavgRes_MS_MS', 'MoI0Rtd_GLS_GLS', 'devLength_MS', 'MoI0Rtd_MS_MS', 'MoI1Rtd_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'scaleLength_MS', 'CoM1Ytd_GLS_GLS', 'MoI1YavgRes_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'tailLength_MS', 'CoM1GstdTail_MS_MS'], -2.0178819444444445), (8, ['mG_MS', 'mRmG_MS', 'aY_GLS', 'mRmG_GLS', 'CoM1RavgRes_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1Gtd_MS_MS', 'mYmG_GLS', 'fracY_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI1Rscale_MS_MS', 'Ytail_GLS', 'mR_MS', 'CoM1Rscale_MS_MS', 'MoI0GavgRes_MS_MS', 'tScale_GLS', 'maxR_GLS', 'MoI1Rscale_GLS_GLS', 'CoM0RstdTail_GLS_GLS', 'Rtail_GLS', 'CoM0Rscale_GLS_GLS', 'CoM0RstdTail_MS_MS', 'MoI1GstdTail_GLS_GLS', 'aR_GLS', 'Gtail_GLS', 'CoM0Ytd_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'sG_GLS', 'MoI1Gtd_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'scaleHead_MS', 'MoI1RstdTail_GLS_GLS', 'aG_MS', 'MoI0Gtd_MS_MS', 'MoI0Rtd_GLS_GLS', 'tailLength_MS', 'MoI1YavgRes_GLS_GLS', 'CoM1GstdTail_MS_MS', 'MoI0Gtd_GLS_GLS', 'MoI0Rtd_MS_MS', 'MoI1Gtd_MS_MS', 'CoM0GstdTail_MS_MS', 'MoI0RstdTail_GLS_GLS', 'MoI1Ytd_GLS_GLS'], -2.0725347222222221), (9, ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'devHead_MS', 'CoM1RavgRes_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'mRmG_GLS', 'aY_GLS', 'mYmG_GLS', 'mRmG_MS', 'MoI1GavgRes_MS_MS', 'fracY_GLS', 'MoI0GavgRes_MS_MS', 'MoI1Rscale_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'Ytail_GLS', 'CoM1Gtd_MS_MS', 'CoM1Rscale_MS_MS', 'CoM1RstdTail_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'tScale_GLS', 'Rtail_GLS', 'mR_MS', 'CoM1RstdTail_MS_MS', 'maxR_GLS', 'MoI1GstdTail_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'sG_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI0Rtd_GLS_GLS', 'tMove_GLS', 'MoI0Gtd_GLS_GLS', 'MoI0GstdTail_MS_MS', 'MoI1Rscale_MS_MS', 'aR_GLS', 'maxG_GLS', 'CoM0Ytd_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'MoI1YavgRes_GLS_GLS', 'CoM1Gscale_MS_MS', 'aG_GLS', 'CoM0YstdTail_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI0Rtd_MS_MS', 'CoM1Ytd_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'scaleLength_MS', 'CoM0RstdTail_MS_MS', 'MoI1Gtd_MS_MS', 'tScale_MS', 'scaleHead_MS', 'MoI1Gtd_GLS_GLS', 'CoM0Rtd_GLS_GLS', 'MoI0RavgRes_GLS_GLS', 'tMove_MS', 'CoM0Gtd_MS_MS', 'MoI1Rtd_MS_MS'], -1.9523958333333336), (10, ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'mRmG_MS', 'CoM1RavgRes_GLS_GLS', 'aY_GLS', 'mRmG_GLS', 'fracY_GLS', 'Ytail_GLS', 'mYmG_GLS', 'CoM1YavgRes_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'devHead_MS', 'MoI1Rscale_GLS_GLS', 'MoI1GavgRes_MS_MS', 'MoI1Gtd_GLS_GLS', 'CoM1Gtd_MS_MS', 'MoI1Yscale_GLS_GLS', 'mR_MS', 'CoM1Rscale_MS_MS', 'MoI0GavgRes_MS_MS', 'CoM1RstdTail_MS_MS', 'tScale_GLS', 'CoM0RstdTail_MS_MS', 'sR_GLS', 'MoI1Gtd_MS_MS', 'MoI1Ytd_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI0Rtd_MS_MS', 'MoI0Rtd_GLS_GLS'], -3.2013888888888888), (11, ['mG_MS', 'mRmG_MS', 'aY_GLS', 'MoI0GavgRes_MS_MS', 'CoM0RstdTail_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'mRmG_GLS', 'fracY_GLS', 'CoM1Gtd_MS_MS', 'mYmG_GLS', 'Gtail_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI1Rscale_MS_MS', 'Ytail_GLS', 'MoI1Gtd_GLS_GLS', 'mR_MS', 'tScale_GLS', 'CoM1Rscale_MS_MS', 'aR_GLS', 'MoI1Rscale_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'Rtail_GLS', 'maxR_GLS', 'MoI0GstdTail_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI1YavgRes_GLS_GLS', 'sG_GLS', 'MoI1GavgRes_GLS_GLS', 'tailLength_MS', 'CoM0Ytd_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'aG_GLS', 'CoM1Rtd_GLS_GLS', 'MoI1Rtd_GLS_GLS', 'CoM0GstdTail_MS_MS', 'MoI1RstdTail_GLS_GLS', 'sY_GLS', 'CoM0YstdTail_GLS_GLS', 'scaleHead_MS', 'MoI0Rtd_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI1Gtd_MS_MS', 'MoI0Rtd_MS_MS', 'MoI1Rtd_MS_MS'], -2.0948958333333332), (12, ['mG_MS', 'mRmG_MS', 'CoM0RstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'MoI1GavgRes_MS_MS', 'MoI0GavgRes_MS_MS', 'aY_GLS', 'mRmG_GLS', 'fracY_GLS', 'CoM1Gtd_MS_MS', 'tScale_GLS', 'mYmG_GLS', 'MoI1Rscale_GLS_GLS', 'mR_MS', 'aR_GLS', 'CoM0Rscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'Ytail_GLS', 'CoM1Rscale_MS_MS', 'CoM1YavgRes_GLS_GLS', 'MoI1Yscale_GLS_GLS', 'MoI1Rscale_MS_MS', 'MoI0RstdTail_MS_MS', 'Gtail_GLS', 'CoM0YstdTail_GLS_GLS', 'MoI0RavgRes_MS_MS', 'scaleLength_MS', 'CoM0RstdTail_MS_MS', 'CoM1Rtd_MS_MS', 'tailLength_MS', 'tScale_MS', 'tMove_GLS', 'CoM0GstdTail_MS_MS', 'CoM1Rtd_GLS_GLS', 'MoI1Gtd_MS_MS', 'MoI1RstdTail_GLS_GLS', 'Rtail_GLS', 'MoI0Gscale_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI0Rtd_GLS_GLS', 'maxG_GLS', 'MoI1GavgRes_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI1Rtd_GLS_GLS'], -1.7673263888888888), (13, ['mG_MS', 'mRmG_MS', 'CoM0RstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'MoI1GavgRes_MS_MS', 'CoM0Rscale_GLS_GLS', 'aY_GLS', 'mRmG_GLS', 'fracY_GLS', 'Ytail_GLS', 'mR_MS', 'CoM1Gtd_MS_MS', 'mYmG_GLS', 'CoM1Rscale_MS_MS', 'MoI1Rscale_GLS_GLS', 'tScale_GLS', 'aR_GLS', 'CoM1YavgRes_GLS_GLS', 'aG_GLS', 'maxR_GLS', 'Gtail_GLS', 'MoI0Rtd_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'Rtail_GLS', 'maxG_GLS', 'CoM0Ytd_GLS_GLS', 'MoI1YavgRes_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'CoM1Gscale_MS_MS', 'MoI0Rtd_MS_MS', 'MoI1Gtd_GLS_GLS', 'sG_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI1Rscale_MS_MS', 'CoM1RstdTail_MS_MS', 'MoI0GavgRes_MS_MS', 'scaleHead_MS', 'tailLength_MS', 'aG_MS', 'MoI1Gtd_MS_MS', 'MoI0Gtd_MS_MS', 'MoI0Ytd_GLS_GLS', 'sY_GLS', 'CoM1GstdTail_MS_MS', 'MoI0Gtd_GLS_GLS', 'MoI0RavgRes_GLS_GLS'], -2.2499652777777777), (14, ['mG_MS', 'mRmG_MS', 'tailHead_MS', 'CoM0RstdTail_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'mRmG_GLS', 'aY_GLS', 'fracY_GLS', 'mYmG_GLS', 'MoI1Rscale_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'Ytail_GLS', 'CoM1Gtd_MS_MS', 'fracG_GLS', 'CoM0Rscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'tScale_GLS', 'aG_GLS', 'mR_MS', 'MoI1GavgRes_MS_MS', 'aR_GLS', 'Rtail_GLS', 'MoI0GavgRes_MS_MS', 'CoM1RstdTail_MS_MS', 'MoI1Gtd_GLS_GLS', 'maxR_GLS', 'MoI1Rscale_MS_MS', 'MoI1YstdTail_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'maxG_GLS', 'MoI0GstdTail_GLS_GLS', 'sG_GLS', 'scaleHead_MS', 'MoI1GstdTail_GLS_GLS', 'MoI1Gtd_MS_MS', 'MoI1RstdTail_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI0Rtd_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'MoI1YavgRes_GLS_GLS', 'tMove_MS', 'CoM1Rtd_GLS_GLS', 'sR_GLS', 'MoI0Gtd_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI0Rtd_MS_MS', 'CoM0GstdTail_MS_MS', 'MoI0RstdTail_GLS_GLS', 'aG_MS'], -2.2749867724867725), (15, ['mG_MS', 'mRmG_MS', 'aY_GLS', 'CoM1RavgRes_GLS_GLS', 'mRmG_GLS', 'fracY_GLS', 'mYmG_GLS', 'CoM0RstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'Ytail_GLS', 'tScale_GLS', 'CoM1Rscale_MS_MS', 'aR_GLS', 'tailHead_MS', 'mR_MS', 'MoI1Rscale_GLS_GLS', 'Gtail_GLS', 'MoI1GavgRes_MS_MS', 'maxG_GLS', 'MoI0Gscale_GLS_GLS', 'MoI0Rtd_GLS_GLS', 'MoI1Rtd_MS_MS', 'MoI0Gtd_MS_MS', 'Rtail_GLS', 'fracG_GLS', 'MoI0Rscale_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'MoI0GavgRes_MS_MS', 'scaleLength_MS', 'tScale_MS', 'CoM0YstdTail_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'CoM1GavgRes_MS_MS', 'MoI1Gtd_MS_MS', 'maxR_GLS', 'CoM0Ytd_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'sG_GLS', 'CoM1YstdTail_GLS_GLS', 'MoI1GstdTail_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'scaleHead_MS', 'MoI1YstdTail_GLS_GLS', 'CoM0Rtd_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI0RstdTail_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'MoI1YavgRes_GLS_GLS', 'MoI0Rtd_MS_MS', 'MoI0GstdTail_MS_MS', 'CoM1Rtd_MS_MS', 'tailLength_MS', 'aG_MS', 'tMove_MS', 'CoM1GstdTail_MS_MS', 'MoI1Gtd_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'sY_GLS', 'CoM1Gscale_MS_MS', 'CoM1Rtd_GLS_GLS', 'MoI1Rtd_GLS_GLS'], -1.8357936507936508), (16, ['mG_MS', 'mRmG_GLS', 'CoM0RstdTail_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'MoI1GavgRes_MS_MS', 'aY_GLS', 'fracY_GLS', 'CoM1RavgRes_GLS_GLS', 'CoM1Gtd_MS_MS', 'mRmG_MS', 'CoM0Rscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'mR_MS', 'Ytail_GLS', 'mYmG_GLS', 'CoM1Rscale_MS_MS', 'tScale_GLS', 'MoI1Rscale_GLS_GLS', 'aR_GLS', 'MoI0GavgRes_MS_MS', 'MoI1Yscale_GLS_GLS', 'MoI1Rscale_MS_MS', 'Rtail_GLS', 'CoM0RstdTail_MS_MS', 'maxR_GLS', 'MoI1Gtd_GLS_GLS', 'MoI1GstdTail_GLS_GLS', 'MoI0Gtd_MS_MS', 'MoI0Rtd_MS_MS', 'MoI0YavgRes_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI1Gtd_MS_MS', 'MoI0Rtd_GLS_GLS', 'MoI1Ytd_GLS_GLS'], -2.8524206349206347)]


    # '''this was from the 1/7/18 group opt by rank run- norm by # GS and MS phenotypes, included 3 additional groups'''
    # opt_res = [(0, ['mRmG_GLS', 'mRmG_MS', 'CoM0Gtd_MS_MS', 'tMove_GLS', 'CoM0Ytd_GLS_GLS', 'MoI0Rscale_MS_MS', 'mG_MS', 'tScale_MS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM0YavgRes_GLS_GLS', 'MoI1Rscale_MS_MS', 'MoI0Yscale_GLS_GLS', 'MoI0Rtd_MS_MS', 'scaleLength_MS', 'CoM1Rscale_MS_MS', 'MoI1Gtd_MS_MS', 'sR_GLS', 'CoM0GavgRes_MS_MS', 'tMove_MS', 'CoM0RstdTail_GLS_GLS', 'CoM1Rscale_GLS_GLS', 'CoM1Rtd_MS_MS', 'MoI1Gscale_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'MoI1Rtd_MS_MS', 'MoI0Ytd_GLS_GLS', 'CoM1GavgRes_MS_MS', 'MoI1Rscale_GLS_GLS', 'devLength_MS', 'MoI1GstdTail_MS_MS', 'CoM0RavgRes_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'CoM0RavgRes_MS_MS', 'MoI0RstdTail_MS_MS', 'MoI0YstdTail_GLS_GLS', 'fracY_GLS', 'MoI1Rtd_GLS_GLS', 'MoI0GstdTail_GLS_GLS'], -2.5616077869096152),(1, ['mRmG_MS', 'mRmG_GLS', 'tScale_MS', 'mG_MS', 'MoI0Rscale_MS_MS', 'CoM0Ytd_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM0RstdTail_GLS_GLS', 'tMove_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI1Ytd_GLS_GLS', 'CoM0GavgRes_MS_MS', 'tMove_MS', 'MoI1Rscale_MS_MS', 'MoI0Yscale_GLS_GLS', 'scaleLength_MS', 'MoI0GstdTail_MS_MS', 'CoM1Rscale_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI0Rtd_MS_MS', 'CoM0Gtd_MS_MS', 'CoM1Rtd_MS_MS', 'MoI1Gscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM0Rscale_GLS_GLS', 'MoI1Rtd_MS_MS', 'MoI0Rtd_GLS_GLS', 'MoI1Gtd_MS_MS', 'MoI0RavgRes_MS_MS', 'sG_GLS', 'MoI0GstdTail_GLS_GLS', 'MoI0YstdTail_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'mR_MS', 'MoI1Rscale_GLS_GLS', 'tScale_GLS', 'CoM0RavgRes_GLS_GLS', 'devLength_MS', 'CoM0RavgRes_MS_MS', 'MoI1YstdTail_GLS_GLS', 'Ytail_GLS'], -2.4001894057966688), (2, ['mRmG_GLS', 'mRmG_MS', 'tScale_MS', 'mG_MS', 'MoI0Rscale_MS_MS', 'CoM0Ytd_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI0GstdTail_MS_MS', 'CoM1Rscale_MS_MS', 'MoI1RstdTail_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'MoI1Rscale_MS_MS', 'MoI0GavgRes_MS_MS', 'scaleLength_MS', 'MoI1Gtd_MS_MS', 'CoM0GavgRes_MS_MS', 'CoM0YavgRes_GLS_GLS', 'tMove_MS', 'sR_GLS', 'CoM0RstdTail_GLS_GLS', 'MoI0RavgRes_MS_MS', 'MoI0Rtd_MS_MS', 'tMove_GLS', 'CoM0Rscale_GLS_GLS', 'CoM1GstdTail_MS_MS', 'CoM1Rscale_GLS_GLS', 'CoM0Gtd_MS_MS', 'MoI0GstdTail_GLS_GLS', 'MoI1Rtd_MS_MS', 'CoM1Rtd_MS_MS', 'MoI1Gscale_GLS_GLS', 'sG_MS', 'CoM1Gscale_MS_MS', 'mR_MS', 'MoI1Rtd_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'MoI1GstdTail_MS_MS', 'devLength_MS', 'CoM0RavgRes_GLS_GLS', 'CoM1Yscale_GLS_GLS', 'tScale_GLS', 'CoM0RavgRes_MS_MS', 'MoI0RstdTail_MS_MS', 'aG_MS'], -2.4904831113651462), (3, ['mRmG_MS', 'tScale_MS', 'mG_MS', 'MoI0Rscale_MS_MS', 'CoM0Ytd_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'tMove_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM0Gtd_MS_MS', 'MoI0GstdTail_MS_MS', 'CoM0YavgRes_GLS_GLS', 'CoM0RstdTail_GLS_GLS', 'MoI1Rscale_MS_MS', 'MoI0Yscale_GLS_GLS', 'tMove_MS', 'CoM0GavgRes_MS_MS', 'MoI1Gtd_MS_MS', 'MoI1Ytd_GLS_GLS', 'scaleLength_MS', 'CoM1Rscale_MS_MS', 'MoI0Rtd_MS_MS', 'CoM0Rscale_GLS_GLS', 'sR_GLS', 'CoM1Rscale_GLS_GLS', 'CoM0RavgRes_MS_MS', 'MoI1RstdTail_GLS_GLS', 'CoM1Rtd_MS_MS', 'MoI1Gscale_GLS_GLS', 'MoI1Rtd_MS_MS', 'sG_MS', 'mR_MS', 'tScale_GLS', 'MoI0RstdTail_MS_MS', 'CoM0RavgRes_GLS_GLS', 'MoI1Rtd_GLS_GLS', 'CoM1GavgRes_MS_MS', 'MoI1Rscale_GLS_GLS', 'fracY_GLS', 'sG_GLS', 'mYmG_GLS', 'MoI1Gtd_GLS_GLS', 'MoI1GstdTail_MS_MS', 'MoI0YstdTail_GLS_GLS', 'aY_GLS', 'devLength_MS', 'CoM1Rtd_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'MoI0RstdTail_GLS_GLS', 'tailHead_MS', 'MoI0Ytd_GLS_GLS', 'MoI1Yscale_GLS_GLS', 'tailLength_MS', 'scaleHead_MS', 'MoI1RavgRes_GLS_GLS', 'CoM1Gscale_MS_MS', 'MoI0RavgRes_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'CoM0GstdTail_MS_MS', 'MoI1YstdTail_GLS_GLS'], -2.2629684494108497), (4, ['mRmG_GLS', 'mRmG_MS', 'CoM0Ytd_GLS_GLS', 'MoI0Rscale_MS_MS', 'mG_MS', 'CoM1RstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1Rscale_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'CoM0RstdTail_GLS_GLS', 'MoI1Rscale_MS_MS', 'MoI0Yscale_GLS_GLS', 'CoM0GavgRes_MS_MS', 'tMove_MS', 'MoI0Rtd_MS_MS', 'CoM0YavgRes_GLS_GLS', 'CoM1Rscale_MS_MS', 'tMove_GLS', 'scaleLength_MS', 'tScale_MS', 'MoI0GavgRes_MS_MS', 'MoI1Gtd_MS_MS', 'CoM0Rscale_GLS_GLS', 'devLength_MS', 'MoI1Rtd_GLS_GLS', 'MoI1Rtd_MS_MS', 'CoM1Rtd_MS_MS', 'MoI1Gscale_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'CoM1GavgRes_MS_MS', 'sG_GLS', 'CoM1RavgRes_GLS_GLS', 'MoI1GstdTail_MS_MS', 'CoM0RavgRes_GLS_GLS', 'CoM0Gtd_MS_MS', 'CoM0Rtd_MS_MS', 'MoI1Rscale_GLS_GLS', 'sR_GLS', 'MoI1Gtd_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'MoI0GstdTail_MS_MS', 'fracY_GLS', 'MoI1GstdTail_GLS_GLS', 'CoM0GstdTail_MS_MS'], -2.535771078590944), (5, ['mRmG_GLS', 'mRmG_MS', 'tScale_MS', 'CoM1YstdTail_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'mG_MS', 'MoI0Rscale_MS_MS', 'CoM1RstdTail_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI0Rtd_MS_MS', 'MoI0GavgRes_MS_MS', 'scaleLength_MS', 'CoM0GavgRes_MS_MS', 'MoI1Rscale_MS_MS', 'MoI1Ytd_GLS_GLS', 'CoM1Rscale_GLS_GLS', 'MoI0GstdTail_MS_MS', 'CoM0RavgRes_GLS_GLS', 'tScale_GLS', 'CoM0RstdTail_GLS_GLS', 'tMove_MS', 'MoI1RstdTail_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'MoI1Gtd_MS_MS', 'tMove_GLS', 'CoM1Rtd_MS_MS', 'MoI1Rtd_MS_MS', 'sR_GLS', 'CoM1YavgRes_GLS_GLS', 'MoI1GstdTail_MS_MS', 'CoM0Gtd_MS_MS', 'sY_GLS'], -2.5293501870253765), (6, ['mRmG_GLS', 'mRmG_MS', 'CoM0Gtd_MS_MS', 'sG_GLS', 'CoM0Ytd_GLS_GLS', 'MoI0Rscale_MS_MS', 'mG_MS', 'CoM0RstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'MoI1Rscale_MS_MS', 'CoM1RstdTail_GLS_GLS', 'tMove_MS', 'CoM0YavgRes_GLS_GLS', 'tScale_MS', 'MoI1RstdTail_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'MoI0GstdTail_MS_MS', 'MoI0Yscale_GLS_GLS', 'CoM0GavgRes_MS_MS', 'scaleLength_MS', 'MoI1Ytd_GLS_GLS', 'MoI1RavgRes_GLS_GLS', 'MoI1Rtd_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI0Rtd_MS_MS', 'CoM1Rscale_GLS_GLS', 'CoM1Rtd_MS_MS', 'MoI1Gscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'tMove_GLS', 'CoM0Rscale_GLS_GLS', 'MoI1Gtd_MS_MS', 'MoI1Rtd_MS_MS', 'MoI0RavgRes_MS_MS', 'CoM1Ytd_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'MoI1GstdTail_MS_MS', 'devLength_MS', 'CoM0RavgRes_MS_MS', 'CoM1GavgRes_MS_MS', 'aY_GLS', 'CoM1RavgRes_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'CoM1Gscale_MS_MS', 'CoM0Yscale_GLS_GLS', 'sR_GLS', 'MoI1Gtd_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'MoI0RstdTail_MS_MS', 'CoM1RavgRes_MS_MS', 'tScale_GLS', 'CoM1GstdTail_MS_MS', 'CoM0GstdTail_MS_MS', 'mYmG_GLS', 'MoI0YstdTail_GLS_GLS', 'sG_MS', 'MoI0Rscale_GLS_GLS', 'CoM0Gscale_MS_MS', 'tailHead_MS', 'mR_MS', 'MoI1Yscale_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'MoI0Rtd_GLS_GLS', 'scaleHead_MS', 'MoI0RavgRes_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'MoI1GavgRes_MS_MS'], -2.2863815393485982), (7, ['mRmG_MS', 'mRmG_GLS', 'tScale_MS', 'mG_MS', 'MoI0Rscale_MS_MS', 'CoM0Ytd_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM0RstdTail_GLS_GLS', 'CoM0GavgRes_MS_MS', 'tMove_GLS', 'tMove_MS', 'MoI1Rscale_MS_MS', 'CoM0Gtd_MS_MS', 'MoI1Gtd_MS_MS', 'CoM0RavgRes_GLS_GLS', 'MoI1GstdTail_MS_MS', 'MoI0Rtd_MS_MS', 'MoI0Yscale_GLS_GLS', 'scaleLength_MS', 'MoI0GavgRes_MS_MS', 'MoI1Ytd_GLS_GLS', 'CoM1Rscale_GLS_GLS', 'CoM1Rtd_MS_MS', 'CoM1Rscale_MS_MS', 'CoM0RavgRes_MS_MS', 'CoM0Rscale_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'CoM0GstdTail_MS_MS', 'MoI0YstdTail_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'tScale_GLS', 'MoI1Rtd_MS_MS', 'sR_GLS', 'MoI0RstdTail_MS_MS', 'MoI1Gscale_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'CoM1GavgRes_MS_MS', 'devLength_MS', 'tailHead_MS', 'MoI1Rtd_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'fracY_GLS'], -2.2384217378315729), (8, ['mRmG_GLS', 'mRmG_MS', 'tScale_MS', 'mG_MS', 'MoI0Rscale_MS_MS', 'CoM1Rscale_MS_MS', 'CoM0Ytd_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI1Gtd_MS_MS', 'CoM1RstdTail_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM1YstdTail_GLS_GLS', 'devLength_MS', 'MoI1Ytd_GLS_GLS', 'MoI0RavgRes_MS_MS', 'tMove_MS', 'MoI0Rtd_MS_MS', 'MoI1RstdTail_GLS_GLS', 'CoM0Rtd_GLS_GLS', 'MoI1Rtd_GLS_GLS', 'fracY_GLS'], -2.7030631398373552), (9, ['mRmG_GLS', 'mRmG_MS', 'tScale_MS', 'mG_MS', 'MoI0Rscale_MS_MS', 'CoM0Ytd_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM0Gtd_MS_MS', 'MoI1Ytd_GLS_GLS', 'tMove_MS', 'MoI1Rscale_MS_MS', 'CoM0YavgRes_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI0GstdTail_MS_MS', 'MoI0Yscale_GLS_GLS', 'aY_GLS', 'MoI1YstdTail_GLS_GLS', 'CoM0GstdTail_MS_MS'], -3.0562440621622846), (10, ['mRmG_GLS', 'mRmG_MS', 'tScale_MS', 'mG_MS', 'MoI0Rscale_MS_MS', 'CoM0Ytd_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'tMove_GLS', 'MoI1Gtd_MS_MS', 'MoI1Ytd_GLS_GLS', 'tMove_MS', 'CoM1RstdTail_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI1RstdTail_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'MoI1Rscale_MS_MS', 'MoI0GavgRes_MS_MS', 'scaleLength_MS', 'CoM0RstdTail_GLS_GLS', 'MoI0Rtd_MS_MS', 'CoM1Rscale_GLS_GLS', 'CoM0GavgRes_MS_MS', 'CoM0Rscale_GLS_GLS', 'MoI0GstdTail_MS_MS', 'MoI1GstdTail_MS_MS', 'MoI0GstdTail_GLS_GLS', 'CoM1Rtd_MS_MS', 'MoI1Gscale_GLS_GLS', 'MoI0RstdTail_MS_MS', 'CoM0YavgRes_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'CoM0GstdTail_MS_MS'], -2.7080130582325523), (11, ['mRmG_MS', 'mRmG_GLS', 'CoM0Gtd_MS_MS', 'CoM1RstdTail_GLS_GLS', 'tScale_MS', 'MoI0Rscale_MS_MS', 'mG_MS', 'CoM0Ytd_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM1YstdTail_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'MoI1Rscale_MS_MS', 'MoI0Yscale_GLS_GLS', 'scaleLength_MS', 'tScale_GLS', 'CoM1Rtd_GLS_GLS', 'MoI1RavgRes_GLS_GLS', 'tMove_MS', 'CoM0YavgRes_GLS_GLS', 'MoI1Gtd_MS_MS', 'CoM1Rscale_MS_MS', 'MoI1Rtd_GLS_GLS', 'MoI0GstdTail_MS_MS', 'MoI1Rscale_GLS_GLS', 'MoI1Gscale_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'CoM0GavgRes_MS_MS', 'MoI0Rtd_MS_MS', 'CoM1Rscale_GLS_GLS', 'CoM1Rtd_MS_MS', 'CoM0Rscale_GLS_GLS', 'tMove_GLS', 'CoM1Gscale_MS_MS', 'sR_GLS', 'MoI1Rtd_MS_MS', 'devLength_MS', 'CoM0RavgRes_GLS_GLS', 'MoI0YstdTail_GLS_GLS', 'MoI0RstdTail_MS_MS', 'MoI0RstdTail_GLS_GLS', 'MoI1GstdTail_MS_MS', 'CoM0RstdTail_GLS_GLS', 'CoM0Rtd_GLS_GLS', 'mR_MS', 'fracY_GLS', 'MoI1RstdTail_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'CoM0RavgRes_MS_MS', 'MoI1YstdTail_GLS_GLS', 'CoM0GstdTail_MS_MS'], -2.3006929235862543), (12, ['mRmG_GLS', 'mRmG_MS', 'tScale_MS', 'mG_MS', 'MoI0Rscale_MS_MS', 'CoM0Ytd_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM0Gtd_MS_MS', 'MoI1Ytd_GLS_GLS', 'tMove_MS', 'CoM0RstdTail_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI1Rscale_MS_MS', 'MoI0Yscale_GLS_GLS', 'CoM0GavgRes_MS_MS', 'MoI0GstdTail_MS_MS', 'tMove_GLS', 'scaleLength_MS', 'MoI1Gtd_MS_MS', 'CoM1Rscale_MS_MS', 'MoI0Rtd_MS_MS', 'CoM1Rscale_GLS_GLS', 'CoM1Rtd_MS_MS', 'MoI1Gscale_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'sR_GLS', 'MoI0RstdTail_MS_MS', 'MoI0Ytd_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'MoI1Rtd_MS_MS', 'sG_MS', 'CoM1GavgRes_MS_MS', 'tScale_GLS', 'MoI1Rtd_GLS_GLS', 'CoM1Gscale_MS_MS', 'MoI1YstdTail_GLS_GLS', 'MoI1RavgRes_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'devLength_MS', 'MoI0GstdTail_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'fracG_GLS'], -2.6123976216309144), (13, ['mRmG_GLS', 'mRmG_MS', 'tScale_MS', 'mG_MS', 'MoI0Rscale_MS_MS', 'CoM0Ytd_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM0Gtd_MS_MS', 'MoI1Ytd_GLS_GLS', 'tMove_MS', 'CoM0RstdTail_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI1Rscale_MS_MS', 'MoI0Yscale_GLS_GLS', 'scaleLength_MS', 'MoI1Gtd_MS_MS', 'CoM1Rscale_MS_MS', 'MoI0Rtd_MS_MS', 'sR_GLS', 'CoM0Rscale_GLS_GLS', 'MoI1Gscale_GLS_GLS', 'CoM1Rtd_MS_MS', 'tMove_GLS', 'CoM1Rscale_GLS_GLS', 'MoI0GstdTail_MS_MS', 'MoI1Rscale_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'MoI1Rtd_MS_MS', 'CoM0GavgRes_MS_MS', 'CoM1GavgRes_MS_MS', 'sG_MS', 'CoM0YstdTail_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'CoM1Gscale_MS_MS', 'tScale_GLS', 'CoM0RavgRes_GLS_GLS', 'devLength_MS', 'MoI1Rtd_GLS_GLS', 'fracG_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI1GstdTail_MS_MS', 'CoM1GstdTail_MS_MS'], -2.5926034014386286), (14, ['mRmG_GLS', 'tScale_MS', 'CoM1YstdTail_GLS_GLS', 'CoM1Rscale_MS_MS', 'CoM0Ytd_GLS_GLS', 'mG_MS', 'MoI0GstdTail_MS_MS', 'CoM1Ytd_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM0Gtd_MS_MS', 'MoI1RstdTail_GLS_GLS', 'tScale_GLS', 'CoM0RstdTail_GLS_GLS', 'CoM0GavgRes_MS_MS', 'CoM0YavgRes_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'mRmG_MS', 'tMove_MS', 'MoI1Rscale_MS_MS', 'MoI0GavgRes_MS_MS', 'CoM1Rscale_GLS_GLS', 'MoI0Rtd_MS_MS', 'MoI0Rscale_MS_MS', 'CoM1Rtd_MS_MS', 'MoI1Gscale_GLS_GLS', 'MoI1Gtd_MS_MS', 'CoM0Rscale_GLS_GLS', 'tMove_GLS', 'scaleLength_MS', 'MoI1Rtd_MS_MS', 'MoI1Rscale_GLS_GLS', 'devLength_MS', 'MoI0Ytd_GLS_GLS', 'CoM0RavgRes_GLS_GLS', 'CoM0RavgRes_MS_MS', 'sG_GLS', 'MoI1GstdTail_MS_MS', 'MoI0GstdTail_GLS_GLS', 'MoI0YstdTail_GLS_GLS'], -2.5611298581858346), (15, ['mRmG_MS', 'mRmG_GLS', 'CoM0Ytd_GLS_GLS', 'mG_MS', 'MoI0Rscale_MS_MS', 'CoM1YstdTail_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'tMove_MS', 'scaleLength_MS', 'MoI0Yscale_GLS_GLS', 'MoI1GavgRes_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI0GavgRes_MS_MS', 'tScale_MS', 'MoI1Gscale_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'MoI1Rtd_MS_MS', 'sR_GLS', 'MoI1RavgRes_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'MoI0RavgRes_MS_MS', 'MoI1Gtd_MS_MS', 'MoI1Rscale_GLS_GLS', 'CoM0Rtd_GLS_GLS', 'CoM0Gtd_MS_MS', 'MoI1Rscale_MS_MS', 'tMove_GLS', 'CoM0RstdTail_GLS_GLS', 'CoM0RavgRes_MS_MS', 'CoM1Ytd_GLS_GLS', 'maxR_GLS', 'CoM0GstdTail_MS_MS', 'CoM1GstdTail_MS_MS', 'CoM1RstdTail_GLS_GLS', 'MoI1GstdTail_MS_MS'], -2.1628162220998322), (16, ['mRmG_MS', 'mRmG_GLS', 'tScale_MS', 'mG_MS', 'MoI0Rscale_MS_MS', 'CoM0Ytd_GLS_GLS', 'CoM1YstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'CoM0Gtd_MS_MS', 'tMove_MS', 'CoM1Rscale_MS_MS', 'MoI1Rscale_MS_MS', 'MoI0GavgRes_MS_MS', 'MoI0Yscale_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'tMove_GLS', 'CoM0GavgRes_MS_MS', 'CoM0RstdTail_GLS_GLS', 'MoI0GstdTail_MS_MS', 'CoM0Rscale_GLS_GLS', 'scaleLength_MS', 'MoI1Gtd_MS_MS', 'MoI0Rtd_MS_MS', 'MoI1Gscale_GLS_GLS', 'CoM1Rtd_MS_MS', 'CoM1Rscale_GLS_GLS', 'MoI1Rtd_MS_MS', 'MoI1Ytd_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'tScale_GLS', 'MoI1Rtd_GLS_GLS', 'devLength_MS', 'sG_GLS', 'CoM0RavgRes_MS_MS', 'CoM1GavgRes_MS_MS', 'CoM0Gscale_MS_MS', 'MoI1Gtd_GLS_GLS'], -2.6415984523421994)]

    p_stats = []
    p_names = []
    for ind, pars, res in opt_res:
        for par in pars:
            if par not in p_names:
                p_names.append(par)
                p_stats.append(1)
            else:
                p_stats[np.where(np.array(p_names) == par)[0][0]] += 1
    p_names = np.array(p_names)
    p_stats = np.array(p_stats)
    sort_inds = np.argsort(p_stats)
    j = 0
    for i in sort_inds[::-1]:
        if p_stats[i] >= th:
            # print p_stats[i], p_names[i]
            j += 1
    print 'total number of parameters to exclude {}'.format(j)
    print(repr(p_names[sort_inds[::-1]]))  # repr prints with comma separation. ::-1 gives the reverse order
    print(repr(p_stats[sort_inds[::-1]]))
    fig = myFigure()
    fig.hist(p_stats)
    fig.show()


def get_PAD_matrix():
    """
    Calculates PAD values for all conditions listed in rnaiList
    :return: pad_array: a matrix of all PAD values to be used for tSNE
    """
    rnaiList, labels, simGroups = get_rnai_groups_pheno_distinct()
    rna_objs = get_rnai_groups_obj(rnaiList)
    group_length = []
    rna_objs_list = []
    for group in rna_objs:
        group_length.append(len(group))
        for obj in group:
            rna_objs_list.append(obj)
    pad_array = get_all_PAD2EO(rna_objs_list)
    pad_array = np.array(pad_array)
    print('finished retrieving PAD')
    return pad_array.reshape((int(np.sqrt(pad_array.size)), int(np.sqrt(pad_array.size)))), group_length, labels


def get_dist_matrix():
    """
    Calculates dist values for all conditions listed in rnaiList
    :return: dist_array: a matrix of all dist values to be used for tSNE
    """
    rnaiList, labels, simGroups = get_rnai_groups_pheno_distinct()
    rna_objs = get_rnai_groups_obj(rnaiList)
    group_length = []
    rna_objs_list = []
    for group in rna_objs:
        group_length.append(len(group))
        for obj in group:
            rna_objs_list.append(obj)
    dist_array = get_all_dist2EO(rna_objs_list)
    dist_array = np.array(dist_array)
    print('finished retrieving dist')
    return dist_array.reshape((int(np.sqrt(dist_array.size)), int(np.sqrt(dist_array.size)))), group_length, labels


def get_tSNE(array, seed):
    from sklearn.manifold import TSNE
    X_embedded = TSNE(n_components=2, random_state=seed, method='exact').fit_transform(array)
    return X_embedded


def show_tsne(points, group_length, labels):
    """
    Plots t_sne for man groups
    :param points: 2D point positions
    :param group_length: length of each group
    :param labels: group labels
    :return:
    """
    fig = myFigure()
    l0 = 0
    for i in range(len(group_length)):
        for point in points[l0:l0 + group_length[i] - 1]:
            fig.scatter(point[0], point[1], color=fig.colors[i % len(fig.colors)])
        fig.scatter(points[l0 + group_length[i] - 1][0], points[l0 + group_length[i] - 1][1],
                    color=fig.colors[i % len(fig.colors)],
                    label=labels[i])
        l0 += group_length[i]
    fig.legend(1)
    # fig.set_axes_equal()
    return fig


def verify_tsne(org, tsne):
    tsne_dist = np.zeros_like(org)
    for i in range(tsne.shape[0]):
        for j in range(tsne.shape[0]):
            tsne_dist[i, j] = np.linalg.norm(tsne[i] - tsne[j])
    print('eCFS', (tsne_dist - org)[:4, :4])
    print('wt', (tsne_dist - org)[-7:, -7:])
    print('tsne', tsne_dist[:4, :4], tsne_dist[-7:, -7:])
    fig = myFigure()
    # fig.hist(tsne_dist.ravel())
    fig.hist(org.ravel())
    fig.show()

def show_clustermap():
    import pandas as pd
    import seaborn as sns
    from string import ascii_letters
    import matplotlib.pyplot as plt
    # df = pd.read_csv('Z:/phenoDist_rank_thr4_test.csv')
    # rnaiList = [21, 364, 408, 130, 264, 386, 31, 32, 255, 388, 422, 118, 359,417, 64, 115, 281, 77, 63, 19, 501,181, 182, 357, 184, 185, 154, 363, 117, 45, 447, 108,9, 398, 435, 45, 52, 38, 95, 4, 57, 5, 277,26, 453, 235, 327, 489, 379, 420, 225, 289, 261, 350,10, 15, 217, 18, 28, 439, 177, 291, 209, 110, 142, 98, 101, 91, 186,383, 495, 498, 414, 375, 396, 321]

    rnaiList = [21, 364, 408, 130, 264, 386, 31, 32, 255, 388, 422, 118, 359, 417, 64, 115, 281, 77, 63, 19, 501, 77,181, 182, 357, 184, 185, 154, 363, 117, 45, 447, 108, 34, 16, 31, 264, 9, 398, 435, 45, 52, 38, 95, 4,57, 5, 277, 67, 90, 235, 403, 503, 261, 404, 25, 26, 453, 235, 327, 489, 379, 420, 225, 289, 261,350, 10, 15, 217, 18, 28, 439, 177, 291, 209, 110, 142, 98, 101, 383, 495, 498, 414, 375, 396,321, 387, 385, 422, 31, 320, 58, 125, 288, 426, 197, 1, 76]  # all

    r_list = np.unique(rnaiList)
    print(len(r_list))
    # r_list = [i for i in range(100)]
    rna_obj_list = []
    rna_obj_labels = []
    for i in r_list:
        print('generating', i)
        r = RNAiClass(i)
        label = r.label
        rna_obj_list.append(r)
        rna_obj_labels.append(r.label)
    datalist = get_all_PAD2EO(rna_obj_list)
    # datalist = get_all_dist2EO(rna_obj_list)

    data = np.reshape(datalist,(len(r_list),len(r_list)))
    columns = rna_obj_labels

    # df = pd.read_csv('Z:/phenoDistCalc/test.csv')
    # print(df.head(5))
    # new_df = df[["RNA1","RNA2","PAD"]]
    # mylambda = lambda x: x.split("D")[1]
    # new_df['RNA1_num'] = df.RNA1.apply(mylambda)
    # new_df['RNA2_num'] = df.RNA2.apply(mylambda)
    # print(new_df.head(5))
    #
    # data = np.reshape(new_df["PAD"],(503,503))
    # columns = np.unique(new_df['RNA2_num'])
    # index = np.unique(new_df['RNA1_num'])
    #
    # print(data)
    # print(columns)
    # print(index)

    sns.set(style="white")
    # d = pd.DataFrame(data=(new_df["PAD"]), columns=(new_df['RNA2_num']))
    d= pd.DataFrame(data=data, columns=columns)
    print(d)
    print("calc corr")
    # d = pd.DataFrame(data = ([[1,0.5,0.2],[0.1, 1, 0.7], [0.3,0.5, 1]]), columns=['a','b','c'], index=['aa','bb','cc']) #this works
    # rs = np.random.RandomState(33)
    # d = pd.DataFrame(data = rs.normal(size=(26,26)), columns=list(ascii_letters[26:]))#this works

    corr = d.corr()
    print("moving on")
    # mask = np.zeros_like(corr,dtype=np.bool)
    # mask[np.triu_indices_from(mask)] = True
    f, ax = plt.subplots(figsize=(15,12))
    cmap = sns.diverging_palette(220,10,as_cmap=True)
    # sns.heatmap(corr, mask=False, cmap=cmap, vmax=.3, center=0,
    #         square=True, linewidths=.5, cbar_kws={"shrink": .5})
    sns.set(font_scale=0.4)
    sns.clustermap(d, metric='euclidean', method='average')
    plt.show()


def show_heatmap():
    '''
    plots a heatmap readout of the defects scored in the 50 gene pilot data set
    :return:
    '''
    import pandas as pd
    import seaborn as sns
    from string import ascii_letters
    import matplotlib.pyplot as plt

    all = True # True makes map for all data, False makes map for CFS/P data only

    if all:  # for All 50 genes heatmap
        data = pd.read_csv('Z:/pilot_for_heatmapAllSens.csv')
        print(data.head(5))
        labels = data[data.columns[2]].astype(str)  # retrieves gene names
        data = data[data.columns[9:]].astype(float)
        f, ax = plt.subplots(figsize=(12, 12))
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        sns.set(style="white")
        sns.heatmap(data, mask=False, cmap="RdBu", robust=False, vmax=3, center=0,
                    square=False, linewidths=.5, cbar_kws={"shrink": .5}, xticklabels='auto', yticklabels=labels)  # cmap="PRGn" (Green), cmap="RdBu" (Blue), cmap=cmap (Red)
        # sns.set(font_scale=0.4)
        plt.show()

    else: # for CFS/P heatmap
        data = pd.read_csv('Z:/cfsp_for_heatmap.csv')
        print(data.head(5))
        labels = data[data.columns[0]].astype(str)  # for cfsp heatmap
        data = data[data.columns[1:4]].astype(float)
        f, ax = plt.subplots(figsize=(4, 6))
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        sns.set(style="white")
        sns.heatmap(data, mask=False, cmap="bwr", robust=False, vmax=50, vmin=-50, center=0,
                    square=False, linewidths=.5, cbar_kws={"shrink": .5}, xticklabels='auto', yticklabels=labels)
        # sns.set(font_scale=0.4)
        plt.show()

def check_dist2zero():
    """
    Check distance to origin for each group and prints group names in order
    :return:
    """
    rlist_test, rlist_train, simGroups_test, simGroups_train, labels_train, labels_test, defect_label_train, \
    defect_label_test = getRNAiGroupsPheno()
    dists = []
    for rs in rlist_train:
        r_dist = []
        for r in rs:
            r_obj = RNAiClass(r)
            r_dist.append(r_obj.getDist2Zero())
        dists.append(np.mean(r_dist))
    inds = np.argsort(dists)
    for i in inds:
        print('{} distance = {}'.format(labels_train[i], dists[i]))

def test():
    print("testing")

def find_closest_groups(embd, group_ave_vects, group_names, makePlot=False, makeDataframe=False, populate_mySQL=True):  # IN PROGRESS- Problem with eCFS_R group!!
    '''
    for a queried RNAi condition, checks closest curated group and returns a continuous readout of PADs between average
    vector for each group and queried gene of interest. This program does several things. If makePlot, it generates
    a plot of all manual groups and the PAD values to the center of those groups for the seed gene (embd). If
    populate_mySQL it writes a description of the closest manual group(s) to the sql. If makeDataframe, it makes a
    nice row for saving to csv.
    :param embd: rnai condition
    :param group_ave_vects: list of vectors representing each of the curated groups
    :param group_names: group names
    :return: depends on what is enabled. May return a plot, may return a df of PAD values and labels for the best
    matched groups or may return an SQL update command.
    '''

    import pandas as pd
    import matplotlib.pyplot as plt
    from embdFunc import getGeneCommonName, getNDDistance

    r = RNAiClass(embd)  # retrieve RNAi object
    r.setNDPosition()
    rows_list = []
    print_description = True


    for i in range(len(group_ave_vects)): # find pad, dist, cos to the ave for each group, save in dict.
        d = {}
        vect = group_ave_vects[i]  # average vector for a manual group
        name = group_names[i] # name of manual group (passed from argument)
        dist = getNDDistance(r.posND, vect, get_params_use())  # distance between gene of interest and ave vector for group
        if np.isinf(dist):
            dist = 0.
        pad = 1. - dist / (r.getDist2Zero() + getNDDistance(vect, 'control', get_params_use()))
        if pad > 0:
            print pad
        else:
            pad = np.nan
        d= {'rna':embd, 'group':name, 'PAD':pad, 'Dist':dist}  # generates a dictionary for each group name
        rows_list.append(d) # adds a dictionary to rows_list for each RNAi-group comparison (this is easily converted to DF)
    df = pd.DataFrame(rows_list, columns=['rna','group','PAD','Dist'])  # makes a dataframe from dictionaries

    # group, group_desc = get_closest_group_description(df)
    group, group_desc = get_closest_group_description(df, print_description=True)


    if makePlot:  # show bar graph of RNAi proximity to manual groups
        df = df.sort_values('PAD', ascending=False)  # sorts values by PAD for display
        print(df.head(3)) # print top 3 categories
        df.plot(kind='bar',x='group', y='PAD', color='red', rot=45) # make a bar chart of closest matched groups by PAD
        plt.title('EMBD{0}, {1}(RNAi)'.format(embd, getGeneCommonName(embd)), fontsize=12, y=1.0, color='k')  # makes title
        plt.xlabel('PAD', fontsize=10)
        plt.show()

    if makeDataframe: # needed to save to csv
        piv_df = df.pivot(index='rna', columns='group', values='PAD') # pivots table to display all pads in one row per RNAi condition
        # print(piv_df)
        piv_df2 = piv_df.assign(description = [group_desc])

        return(piv_df2)

    if populate_mySQL: # needed to update column in mySQL
        from db_utils_embryos import initialize_db
        conn, curs = initialize_db()
        sql = "SELECT id, embd_id, movement_GLS, movement_MS, auto_description FROM genes WHERE embd_id={}".format(
            embd)
        curs.execute(sql)
        id, embd_id, movementG, movementM, auto_description = curs.fetchall()[0]
        if movementG == 0 or movementM == 0:
            sql2 = "UPDATE genes SET auto_description=\"{p}\" WHERE id={id};".format(p=group_desc, id=id)
            curs.execute(sql2)
            conn.commit()
        elif group[0] == 'dev delay':
            sql2 = "UPDATE genes SET auto_description=\"{p}\" WHERE id={id};".format(p=group_desc, id=id)
            curs.execute(sql2)
            conn.commit()
        elif movementG == 1 and movementM == 1 and group [0] is not 'dev delay':
            sql2 = "UPDATE genes SET auto_description=\"{p}\" WHERE id={id};".format(p='Wild-type. Embryos develop'
                                                                                       ' normally and sustain movement '
                                                                                       'within eggshell (or hatch).', id=id)
            curs.execute(sql2)
            conn.commit()


def get_closest_group_description(df, print_description=False):
    '''
    finds closest manual group and pulls description for that group. To be plotted or added to csv or sql
    :param df: dataframe df = has 'rna','group','PAD','Dist' from find_closest_groups
    :return: group and group description. Both are lists
    '''
      # retrieves a brief description based on closest manual group
    df = df.sort_values('PAD', ascending=False)  # sorts by PAD values to return highest PAD first
    df= df.reset_index(drop=True)  # resets index position based on sorted order
    ind = df.PAD.idxmax()  # finds the index of the highest PAD
    max = df.PAD.max()  # finds the PAD value of the highest PAD
    inds = df.index[df['PAD'] > max - max * 0.05].tolist()  # returns a list of index positions that are at or near the top PAD value
    group = [df.group[ind]]  # top matched group

    poor_match_flag = False
    moving_flag = False
    if df.PAD[ind] < 0.5:  # checks to see if top matched group is above 0.5 PAD, if not returns 'poor match' indicator
        if df.group[ind]== 'dev delay' or df.group[ind]=='wt' and df.PAD[ind] > 0.35: # WT and dev delay are more spread out groups- if they are top groups, permit a lower threshold bc they are more spread out (moving)
            group = group
            moving_flag = True
        else:
            group = ['poor match']
            poor_match_flag = True
    if len(inds)>1 and not poor_match_flag: # checks to see if there are more than one group with high PAD values, if so, uses that list as groups
        group = df.group[[i for i in inds]].tolist()
    if df.PAD[ind] < 0.55 and not poor_match_flag and not moving_flag:  # only retrieves top group if that group PAD is between 0.5-0.55. Want to avoid getting
        # second best group with a PAD of <0.5
        group = [df.group[ind]]  # top matched group
    group_desc = get_group_description(group)  # returns a list of descriptions based on group name(s)
    if print_description:
        print((group, group_desc))
    return(group, group_desc)


def get_group_description(group):
    '''
    provides a description of the closest manual group
    :param group: manual group distinction- a list of strings, usually one, i.e. ['2xCrunch'], i.e. ['2xCrunch', rupture].
    List consists of top group and any group whose PAD score is within 5% of the top group to account for positions in between groups
    :return: a list of descriptions (string)
    '''

    desc_list = []
    for i in group:
        if i == '2xCrunch':
            desc = 'Embryos arrest during late elongation, usually at two-fold, often with irregular epidermal structure'
        elif i == 'CFS_lGlY':
            desc = 'Embryos exhibit cell fate defects, usually with reduced endoderm, and arrest at comma to 1.5 fold' \
                   ' (germ layer reporter strain) or 2-3 fold (morphogenesis strain). '
        elif i == 'CFS_wnt':
            desc = 'Embryos exhibit severe cell fate defects and arrest early with sectored germ layers '
        elif i == 'dev delay':
            desc = 'Embryos develop at a slow pace.'
        elif i == 'dorsRupt':
            desc = 'Embryos exhibit cell fate defects, with absent/reduced epidermis, and arrest at comma to 1.5 fold stage. '
        elif i == 'dorsBend':
            desc = 'Embryos initiate elongation, but usually bend backwards (or twist) and arrest near 1.5 fold stage, often' \
                   ' with partial rupture.'
        elif i == 'eCFS_R':
            desc = 'Embryos exhibit severe cell fate defects and arrest early with excess ectoderm.'
        elif i == 'eMEX':
            desc = 'Embryos exhibit severe cell fate defects and arrest early with excess mesoderm.'
        elif i == 'gMix':
            desc = 'Embryos exhibit a variable point of arrest in the germ layer reporter strain. Embryos are mostly normal ' \
                   'in morphogenesis strain.'
        elif i == 'gRupt':
            desc = 'Embryos expressing germ layer markers exhibit variable degrees of rupture- sometimes complete, sometimes' \
                   ' arrest at two-fold with minor anterior rupture. Embryos expressing morphogenesis markers exhibit no defects.'
        elif i == 'hiNeur':
            desc = 'Embryos expressing germ layer markers exhibit elongation defects (usually 1.5 fold, but can vary from' \
                   ' comma to two-fold stage) with mild apparent cell fate/positioning defects. Embryos expressing' \
                   ' morphogenesis markers appear to have excessive neurons and often have epidermal organization defects,' \
                   ' and exhibit elongation defects.'
        elif i == 'lCFS_R':
            desc = 'Embryos exhibit cell fate defects and arrest near comma stage with excess ectoderm.'
        elif i == 'lMEX':
            desc = 'Embryos exhibit severe cell fate defects and arrest near bean/comma stage with excess mesoderm.'
        elif i == 'no_markers':
            desc = 'Embryos arrest early, prior to fluorescent marker expression.'
        elif i == 'rupture':
            desc = 'Embryos exhibit strong rupture phenotype in both reporter strains.'
        elif i == 'sect rupt':
            desc = 'Embryos exhibit severe cell fate defects and usually arrest with sectored germ layer, often attempt and fail enclosure.'
        elif i == 'wt':
            desc = 'Wild-type. Embryos develop normally and sustain movement within eggshell (or hatch).'
        elif i == 'poor match':
            desc = 'Embryo phenotype does not closely match any group. Likely a weak or highly variable phenotype.'
        else:
            desc = 'no group returned'

        desc_list.append(desc)
        description = " | ".join(desc_list)
    return(description)


def get_all_closest_groups(to_csv=False, to_SQL=False, to_plot=False):
    '''
    For each RNAi condition, measures PAD to the ave position for each manual group, generates a df, may write to csv
    make a plot and or update sql 'auto_description column', depending on which kwargs settings are enabled.
    :return: dataframe to csv, plot, or SQL update
    '''
    import pandas as pd
    import time
    from fileLookup import SKIPPED

    group_names, group_ave_vects = calculate_avg_vectors_man_groups()
    all_dfs =[]

    path = 'Z:/closest_groups{}.csv'.format(time.strftime('%Y%m%d'))
    skipped = SKIPPED

    for i in range(37,504):
        if 'EMBD{0}'.format(i) in skipped:
            print("skipped {0}".format(i))
            pass
        else:
            print(i)
            df = find_closest_groups(i, group_ave_vects, group_names,makePlot=to_plot, makeDataframe=to_csv,
                                     populate_mySQL=to_SQL)
            all_dfs.append(df)
    if to_csv:
        all = pd.concat(all_dfs)
        print(all.head(5))
        all.to_csv(path)

def find_closest_group_byEMBD(embd):
    '''
    plots and returns nearest manual group information for a given embd
    :param embd: int value representing RNAi condition, i.e. 1085
    :return: plots best match and returns dataframe containing PAD values to each group centroid
    '''
    labels_list, avg_vect_list= calculate_avg_vectors_man_groups()
    find_closest_groups(embd, avg_vect_list, labels_list, makePlot=True, populate_mySQL=False)

def plot_dist2WT_dist2avg (rnaiList, labels, simGroups):
    '''
    :param rnaiList:
    :param labels:
    :param simGroups:
    :return: plot of dist to wt vs dist to avg within groups
    '''
    fig, dist2wt, dist2wtErr = plotDist2WTvsArrest(rnaiList, labels)
    fig, dist2avg, dist2avgErr = plotDistInGroupVsArrest(rnaiList, labels)
    fig = myFigure()
    fig.errorbar(dist2wt, dist2avg, xerr=dist2wtErr, yerr=dist2avgErr, join=False)
    fig.xlabel('dist 2 wt')
    fig.ylabel('dist 2 avg')
    fig = plotPADInGroup(rnaiList, labels, simGroups)
    # fig.ylim((0.985,1))
    # fig = plotCSIInGroupVsArrest(rnaiList, labels)
    fig.show()

def make_tSNE_all_data():
    """
    Calculates PAD values for all conditions listed in rnaiList
    :return: pad_array: a matrix of all PAD values to be used for tSNE
    IN PROGRESS-- this runs as is...but still need to add in color overlay for query gene and colors, labels for man groups, remove axis

    if int(labels) in X_list:
        np.where(points[0][0]== points[int(labels[0])-1][0]) #checks
    """
    from sklearn.manifold import TSNE
    import pandas as pd
    import time
    from embdFunc import getGeneCommonName

    rnaiList = [i for i in range(1,504)] # a list of EMBD numbers to generate RNAi objects from
    labels = [str(x) for x in rnaiList] # a list of EMBD numbers in string format, as labels-- use to cross ref group numbers to indicies later

    # rna_objs = get_rnai_groups_obj(rnaiList)
    rna_objs = []
    for r in rnaiList:
        rna_objs.append(RNAiClass(r))  # adds RNAi object to rns list
        rna_objs[-1].setNDPosition()  # takes the last RNAi added and sets its ND position

    pad_array = get_all_PAD2EO(rna_objs) #getsPADmatrix
    pad_array = np.array(pad_array)
    print('finished retrieving PAD')
    array = pad_array.reshape((int(np.sqrt(pad_array.size)), int(np.sqrt(pad_array.size))))

    seed = 71076
    X_embedded = TSNE(n_components=2, random_state=seed, method='exact').fit_transform(array)  #getTSNE
    pd.DataFrame(X_embedded).to_csv('Z:\\Automated_analysis\\tSNE\\X_embedded_500genes_{0}_{d}.csv'.format(seed, d=time.strftime('%m%d%Y')), header=None, index=None)
    print('DF Saved')

    points = X_embedded # array of points for plotting tSNE

    fig = myFigure() #show_tsne
    fig.noAxis()
    for point in points: #2D point positions
        fig.scatter(point[0], point[1], color="grey", marker=".", alpha=0.5) #label=labels[i] # plots all points in X_embedded


    'get overlay positions to plot'
    query = [[408]]
    groupLists = [[21, 364, 408, 130], [417, 64, 115, 281, 77, 63, 19, 501], [264, 386, 31, 32, 255, 388, 422, 118, 359],
                [181, 182, 357, 184, 185, 154, 363, 117, 45, 447, 108], [9, 398, 435, 45, 52, 38], [95, 4, 57, 5, 277],
                [26, 453, 235, 327, 489, 379, 420, 225, 289, 261, 350],
                [10, 15, 217, 18, 28, 439, 177, 291, 209], [110, 142, 98, 101, 186],
                [383, 495, 498, 414, 375, 396, 321]]
    # groupLists = [[21], [31, 32], [64,77, 63, 19],
    #             [45], [9, 45, 52, 38], [95, 4, 57, 5],
    #             [26],
    #             [10, 15, 18, 28], [98, 91],
    #             [50]]#use this set to test out without having to wait for full dataset
    # groupLists = [[1,2], [3,4], [5],
    #             [6], [7,8,9], [10],
    #             [11],
    #             [12, 13, 14], [15],
    #             [16]]  #use this set to test out without having to wait for full dataset
    groupLists += query  # adds the query onto the last position in the list
    groupLabels = ['eCFS_R',  'MEX','sect rupt', 'hiNeur', 'dorsBend', 'rupture', 'gMix', '2xCrunch', 'dev delay', 'wt', 'query']
    groupColors = ['sienna',  'orangered','firebrick', 'gold', 'y', 'green', 'turquoise', 'teal', 'dodgerblue', 'b', 'm']
    plot_labels = ['Early Arrest, \n High Ectoderm', 'Early Arrest, \n High Mesoderm',
               'Early Arrest Sectored  \n Germ-Layers', 'Elongation Defect,  \n High Neurons',
              'Enclosure/Early Elong.  \n Dorsal Defect', 'Enclosure/Early Elong. \n Ventral Rupture',
              'Variable \n Defects', 'Late Elongation \n Arrest', 'Dev. Delay', 'Wild-type', 'Query']

    for g in range(len(groupLists)):
        man_group_list = groupLists[g] #individual manual group
        lab = plot_labels[g]
        inds1= [] # this will be populated with the index position for each member of the man group within the X_embedded
        for p in man_group_list:
            inds1.append(labels.index(str(p)))

        for ind in inds1:
            fig.scatter(points[ind][0], points[ind][1], color=groupColors[g]) #plots the individual members of a group in a specific color
        if lab == 'Query':
            r = RNAiClass(query[0][0])
            query_name = getGeneCommonName(query[0][0])
            lab = '{}(RNAi)\n {}'.format(query_name, r.label)
            fig.scatter(np.mean(points[inds1], axis=0)[0], np.mean(points[inds1], axis=0)[1], label=lab, marker="*", s=fig.markerSize*20,           alpha=0.5)
            fig.text(np.mean(points[inds1], axis=0)[0], np.mean(points[inds1], axis=0)[1], lab,
                 color=groupColors[g], size=14, style='italic')
        elif (g % 2)==0: # if g is even, shift label up and over
            fig.scatter(np.mean(points[inds1], axis=0)[0], np.mean(points[inds1], axis=0)[1], label=lab, marker=".",
                        alpha=0.5)
            fig.text(np.mean(points[inds1], axis=0)[0] + 4, np.mean(points[inds1], axis=0)[1] + 4, lab,
                     color=groupColors[g], size=11)
        else:
            fig.scatter(np.mean(points[inds1], axis=0)[0], np.mean(points[inds1], axis=0)[1], label=lab, marker=".",
                        alpha=0.5)
            fig.text(np.mean(points[inds1], axis=0)[0] - 4, np.mean(points[inds1], axis=0)[1] - 4, lab,
                     color=groupColors[g], size=11)
        # if g == 0:
    #             fig.text(np.mean(points[inds1], axis=0)[0]-1, np.mean(points[inds1], axis=0)[1]-1, lab, color=groupColors[g], size=10)
    #         elif g == 2:
    #             fig.text(np.mean(points[inds1], axis=0)[0] - 1, np.mean(points[inds1], axis=0)[1] + 1, lab,
    #                      color=groupColors[g], size=10)
    #         elif (g % 2) == 0:  # if g is even, shift label up and over
    #             fig.text(np.mean(points[inds1], axis=0)[0]+1, np.mean(points[inds1], axis=0)[1]+1, lab, color=groupColors[g], size=10)
    #         else: # if i is odd, shift text down and over
    #             fig.text(np.mean(points[inds1], axis=0)[0]-1, np.mean(points[inds1], axis=0)[1]-1, lab, color=groupColors[g], size=10)

    fig.noClip()
    # fig.legend(1)
    fig.show()
    return fig ######

def make_tSNE_all_data():
    """
    Calculates PAD values for all conditions listed in rnaiList
    :return: pad_array: a matrix of all PAD values to be used for tSNE
    IN PROGRESS-- this runs as is...but still need to add in color overlay for query gene and colors, labels for man groups, remove axis

    if int(labels) in X_list:
        np.where(points[0][0]== points[int(labels[0])-1][0]) #checks
    """
    from sklearn.manifold import TSNE
    import pandas as pd
    import time
    from embdFunc import getGeneCommonName

    rnaiList = [i for i in range(1,504)] # a list of EMBD numbers to generate RNAi objects from
    labels = [str(x) for x in rnaiList] # a list of EMBD numbers in string format, as labels-- use to cross ref group numbers to indicies later

    # rna_objs = get_rnai_groups_obj(rnaiList)
    rna_objs = []
    for r in rnaiList:
        rna_objs.append(RNAiClass(r))  # adds RNAi object to rns list
        rna_objs[-1].setNDPosition()  # takes the last RNAi added and sets its ND position

    pad_array = get_all_PAD2EO(rna_objs) #getsPADmatrix
    pad_array = np.array(pad_array)
    print('finished retrieving PAD')
    array = pad_array.reshape((int(np.sqrt(pad_array.size)), int(np.sqrt(pad_array.size))))

    seed = 71076
    X_embedded = TSNE(n_components=2, random_state=seed, method='exact').fit_transform(array)  #getTSNE
    pd.DataFrame(X_embedded).to_csv('Z:\\Automated_analysis\\tSNE\\X_embedded_500genes_{0}_{d}.csv'.format(seed, d=time.strftime('%m%d%Y')), header=None, index=None)
    print('DF Saved')

    points = X_embedded # array of points for plotting tSNE

    fig = myFigure() #show_tsne
    fig.noAxis()
    for point in points: #2D point positions
        fig.scatter(point[0], point[1], color="grey", marker=".", alpha=0.5) #label=labels[i] # plots all points in X_embedded


    'get overlay positions to plot'
    query = [[408]]
    groupLists = [[21, 364, 408, 130], [417, 64, 115, 281, 77, 63, 19, 501], [264, 386, 31, 32, 255, 388, 422, 118, 359],
                [181, 182, 357, 184, 185, 154, 363, 117, 45, 447, 108], [9, 398, 435, 45, 52, 38], [95, 4, 57, 5, 277],
                [26, 453, 235, 327, 489, 379, 420, 225, 289, 261, 350],
                [10, 15, 217, 18, 28, 439, 177, 291, 209], [110, 142, 98, 101, 186],
                [383, 495, 498, 414, 375, 396, 321]]
    # groupLists = [[21], [31, 32], [64,77, 63, 19],
    #             [45], [9, 45, 52, 38], [95, 4, 57, 5],
    #             [26],
    #             [10, 15, 18, 28], [98, 91],
    #             [50]]#use this set to test out without having to wait for full dataset
    # groupLists = [[1,2], [3,4], [5],
    #             [6], [7,8,9], [10],
    #             [11],
    #             [12, 13, 14], [15],
    #             [16]]  #use this set to test out without having to wait for full dataset
    groupLists += query  # adds the query onto the last position in the list
    groupLabels = ['eCFS_R',  'MEX','sect rupt', 'hiNeur', 'dorsBend', 'rupture', 'gMix', '2xCrunch', 'dev delay', 'wt', 'query']
    groupColors = ['sienna',  'orangered','firebrick', 'gold', 'y', 'green', 'turquoise', 'teal', 'dodgerblue', 'b', 'm']
    plot_labels = ['Early Arrest, \n High Ectoderm', 'Early Arrest, \n High Mesoderm',
               'Early Arrest Sectored  \n Germ-Layers', 'Elongation Defect,  \n High Neurons',
              'Enclosure/Early Elong.  \n Dorsal Defect', 'Enclosure/Early Elong. \n Ventral Rupture',
              'Variable \n Defects', 'Late Elongation \n Arrest', 'Dev. Delay', 'Wild-type', 'Query']

    for g in range(len(groupLists)):
        man_group_list = groupLists[g] #individual manual group
        lab = plot_labels[g]
        inds1= [] # this will be populated with the index position for each member of the man group within the X_embedded
        for p in man_group_list:
            inds1.append(labels.index(str(p)))

        for ind in inds1:
            fig.scatter(points[ind][0], points[ind][1], color=groupColors[g]) #plots the individual members of a group in a specific color
        if lab == 'Query':
            r = RNAiClass(query[0][0])
            query_name = getGeneCommonName(query[0][0])
            lab = '{}(RNAi)\n {}'.format(query_name, r.label)
            fig.scatter(np.mean(points[inds1], axis=0)[0], np.mean(points[inds1], axis=0)[1], label=lab, marker="*", s=fig.markerSize*20,           alpha=0.5)
            fig.text(np.mean(points[inds1], axis=0)[0], np.mean(points[inds1], axis=0)[1], lab,
                 color=groupColors[g], size=14, style='italic')
        elif (g % 2)==0: # if g is even, shift label up and over
            fig.scatter(np.mean(points[inds1], axis=0)[0], np.mean(points[inds1], axis=0)[1], label=lab, marker=".",
                        alpha=0.5)
            fig.text(np.mean(points[inds1], axis=0)[0] + 4, np.mean(points[inds1], axis=0)[1] + 4, lab,
                     color=groupColors[g], size=11)
        else:
            fig.scatter(np.mean(points[inds1], axis=0)[0], np.mean(points[inds1], axis=0)[1], label=lab, marker=".",
                        alpha=0.5)
            fig.text(np.mean(points[inds1], axis=0)[0] - 4, np.mean(points[inds1], axis=0)[1] - 4, lab,
                     color=groupColors[g], size=11)
        # if g == 0:
    #             fig.text(np.mean(points[inds1], axis=0)[0]-1, np.mean(points[inds1], axis=0)[1]-1, lab, color=groupColors[g], size=10)
    #         elif g == 2:
    #             fig.text(np.mean(points[inds1], axis=0)[0] - 1, np.mean(points[inds1], axis=0)[1] + 1, lab,
    #                      color=groupColors[g], size=10)
    #         elif (g % 2) == 0:  # if g is even, shift label up and over
    #             fig.text(np.mean(points[inds1], axis=0)[0]+1, np.mean(points[inds1], axis=0)[1]+1, lab, color=groupColors[g], size=10)
    #         else: # if i is odd, shift text down and over
    #             fig.text(np.mean(points[inds1], axis=0)[0]-1, np.mean(points[inds1], axis=0)[1]-1, lab, color=groupColors[g], size=10)

    fig.noClip()
    # fig.legend(1)
    fig.show()
    return fig ######

def make_UMAP_content_files():
    '''
    generates PAD matrix and distance matrix files (plus labels and colors) to be used for making UMAPs- needs python 3.6 to make the UMAP itself.
    :return: pad_array_reshaped, dist_array_reshaped, groupLists, groupLabels, groupColors, plot_labels
    '''

    from sklearn.manifold import TSNE
    import pandas as pd
    import time
    from embdFunc import getGeneCommonName

    rnaiList = [i for i in range(1, 504)]  # a list of EMBD numbers to generate RNAi objects from
    labels = [str(x) for x in
              rnaiList]  # a list of EMBD numbers in string format, as labels-- use to cross ref group numbers to indicies later

    rna_objs = []
    for r in rnaiList:
        rna_objs.append(RNAiClass(r))  # adds RNAi object to rns list
        rna_objs[-1].setNDPosition()  # takes the last RNAi added and sets its ND position

    pad_array = get_all_PAD2EO(rna_objs)  # getsPADmatrix
    pad_array = np.array(pad_array)
    print('finished retrieving PAD')
    dist_array = get_all_dist2EO(rna_objs)
    dist_array = np.array(dist_array)
    print('finished retrieving dist')

    pad_array_reshaped = pad_array.reshape((int(np.sqrt(pad_array.size)), int(np.sqrt(pad_array.size))))
    dist_array_reshaped = dist_array.reshape((int(np.sqrt(dist_array.size)), int(np.sqrt(dist_array.size))))
    np.save('Z:/PAD_Matrix_msParams.npy', pad_array_reshaped)
    np.save('Z:/Dist_Matrix_msParams.npy', dist_array_reshaped)
    np.save('Z:/Matrix_labels_msParams.npy', labels)


    'overlay position information for plot'
    groupLists = [[21, 364, 408, 130], [417, 64, 115, 281, 77], [77, 63, 19, 501],
                  [264, 386, 31, 32, 255, 388, 422, 118, 359],
                  [181, 182, 357, 184, 185, 154, 363, 117, 45, 447, 108], [9, 398, 435, 45, 52, 38],
                  [95, 4, 57, 5, 277],
                  [26, 453, 235, 327, 489, 379, 420, 225, 289, 261, 350],
                  [10, 15, 217, 18, 28, 439, 177, 291, 209], [110, 142, 98, 101, 186],
                  [383, 495, 498, 414, 375, 396, 321]]
    # groupLists = [[21], [31, 32], [64,77, 63, 19],
    #             [45], [9, 45, 52, 38], [95, 4, 57, 5],
    #             [26],
    #             [10, 15, 18, 28], [98, 91],
    #             [50]]#use this set to test out without having to wait for full dataset
    # groupLists = [[1,2], [3,4], [18], [5],
    #             [6], [7,8,9], [10],
    #             [11],
    #             [12, 13, 14], [15],
    #             [16]]  #use this set to test out without having to wait for full dataset
    groupLabels = ['eCFS_R', 'eMEX', 'lMEX', 'sect rupt', 'hiNeur', 'dorsBend', 'rupture', 'gMix', '2xCrunch',
                   'dev delay', 'wt-like']
    groupColors = ['sienna', 'orangered', 'firebrick', 'darkorange', 'gold', 'y', 'green', 'turquoise', 'teal',
                   'dodgerblue', 'b']
    plot_labels = ['Early Arrest, \n High Ectoderm', 'Early Arrest, \n High Mesoderm',
                   'Mid Arrest, \n High Mesoderm',
                   'Early Arrest Sectored  \n Germ-Layers', 'Elongation Defect,  \n High Neurons',
                   'Enclosure/Early Elong.  \n Dorsal Defect', 'Enclosure/Early Elong. \n Ventral Rupture',
                   'Variable \n Defects', 'Late Elongation \n Arrest', 'Dev. Delay', 'Wild-type']
    np.save('Z:/group_Lists_msParams.npy', groupLists)
    np.save('Z:/group_Labels_msParams.npy', groupLabels)
    np.save('Z:/group_colors_msParams.npy', groupColors)
    np.save('Z:/group_Labels_toPlot_msParams.npy', plot_labels)

    return(pad_array_reshaped, dist_array_reshaped, groupLists, groupLabels, groupColors, plot_labels)


def initialize():
    # array, group_length, labels = get_PAD_matrix()
    # array, group_length, labels = get_dist_matrix()

    # tsne = get_tSNE(array, 71076)
    # fig1 = show_tsne(tsne, group_length, labels)
    # tsne = get_tSNE(array, 12483)
    # fig2 = show_tsne(tsne, group_length, labels)
    # tsne = get_tSNE(array, 20715)
    # fig3 = show_tsne(tsne, group_length, labels)
    # for seed in range(15):  # note that tsne doesnt plot correctly when run on this workstation-- run from BG computer or RK computer
    #     tsne = get_tSNE(array, seed)
    #     fig1 = show_tsne(tsne, group_length, labels)
    #     fig1.title('seed = {}'.format(seed))
    # fig1.show()
    # print('FINISHED')

    make_tSNE_all_data()


    '''evaluate optimization'''
    # pg, pm = get_params_use()
    # print('number of parameters used GLS = {}, MS = {}'.format(len(pg), len(pm)))

    # rlist_test, rlist_train, simGroups_test, simGroups_train, labels_train, labels_test, defect_label_train, defect_label_test, rnaiList, labels, simGroups = getRNAiGroupsPheno()

    # print('test list = {0}'.format(rlist_test))
    # print('labels = {0}'.format(labels_test))
    # print('test df label = {0}'.format(defect_label_test))
    # print('train list = {0}'.format(rlist_train))

    'this overwrites the test set and removes dev delay group, which had a bad member (91- this has half dev delay and half early arrest) which was skewing ranks'
    rlist_test = [[31, 32, 255], [77, 281], [184, 447, 185], [435, 9], [5, 57], [261, 25, 235], [489, 235, 26, 350],
            [291, 439, 209], [414, 375, 383]]  # removed [186, 91] dev delay group because 91 is not a good member
    defect_label_test = ['cf', 'cf', 'm', 'm', 'm', 'm', 'm', 'm', 'm']
    labels_test = ['sect rupt', 'eMEX', 'hiNeur', 'dorsBend', 'rupture', 'gRupt', 'gMix', '2xCrunch', 'wt']

    # evaluate_ranking(rlist_test, labels_test, defect_label_test)
    # evaluate_ranking(rlist_train, labels_train, defect_label_train)


    # '''run optimization on groups- by rank'''
    # rlist_test, rlist_train, simGroups_test, simGroups_train, labels_train, labels_test, defect_label_train, defect_label_test = getRNAiGroupsPheno()
    # pool = mp.Pool(processes=6)
    # results_all = pool.map(optimize_groups_by_rank,
    #                        zip(repeat(rlist_train), repeat(labels_train), repeat(simGroups_train), repeat(defect_label_train),
    #                            [None]+range(len(rlist_train))))
    # printLog(results_all)


#     # optimize_groups_by_min_max(rnaiList, labels, simGroups)
#     plot_dist2WT_dist2avg(rnaiList, labels, simGroups)
#
    # fig, dist2wt, dist2wtErr = plotDist2WTvsArrest(rnaiList, labels)
#     fig, dist2avg, dist2avgErr = plotDistInGroupVsArrest(rnaiList, labels)
#     fig = myFigure()
#     fig.errorbar(dist2wt, dist2avg, xerr=dist2wtErr, yerr=dist2avgErr, join=False)
#     fig.xlabel('dist 2 wt')
#     fig.ylabel('dist 2 avg')

    # fig = plotPADInGroup(rnaiList, labels, simGroups)

    # fig.ylim((0.985,1))
    # fig = plotCSIInGroupVsArrest(rnaiList, labels)

    # fig.show()



if __name__ == '__main__':
    # initialize()
    print('starting')
    # check_dist2zero()
    # cross_corr_params(0)

    # rnas = [320, 398, 54, 95, 59, 397, 288, 177, 336, 10, 52, 209, 217, 363, 277, 447, 184, 123, 117, 408, 125,
    #         130, 364, 13, 359, 21, 77]
    # rnai_list, labels, sim_groups, defect_label = getRNAiGroupsPheno()
    # r_flat = [item for sub in rnai_list for item in sub]
    # for r in rnas:
    #     if r not in r_flat:
    #         print('RNAi # {} is not in training'.format(r))

    # test()
    # show_clustermap()
    # show_heatmap()

    '''find the closest matching manual group'''
    labels_list, avg_vect_list= calculate_avg_vectors_man_groups()
    find_closest_groups(6, avg_vect_list, labels_list, makePlot=True, populate_mySQL=False)

    # get_all_closest_groups(to_csv=False, to_SQL=True, to_plot=False)

    # conditions = [i for i in range(1,504)]
    # conditions.append(1138)
    # # conditions.append(3025)
    # rnaList = []
    # from RNAiClass import RNAiClass
    # for gene in conditions:
    #     r = RNAiClass(gene)
    #     rnaList.append(r)
    # print("step1 done")
    # arraylist = getPAD2EO(rnaList)
    # arraylist=np.array(arraylist)
    # select = arraylist[np.where(arraylist[:,1] == "EMBD1138")[0]]
    # print(select)
    # select2 = arraylist[np.where(arraylist[:,0] == "EMBD1138")[0]]
    # print(select2)
    # print("done")

    # find_closest_group_byEMBD(419)
    # testing_importances()
    # labels, data = retrieve_and_map_man_group_autodata()
    # test_two_sample_ks_test(labels, data)
    # log_regression()
    # x_train, x_test, y_train, y_test = run_man_group_classifier()
    # search_cv(x_train, x_test, y_train, y_test)

    # find_closest_group_byEMBD(1903)
    # make_UMAP_content_files()