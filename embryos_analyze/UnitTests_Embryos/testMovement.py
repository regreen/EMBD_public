'''
Created on Jun 24, 2016

@author: Admin
'''

from varLookup import FOLDER_IN, printLog
import csv
from emb_handler import loadEmbs, refreshEmbsMulti
from myFunc import getStrVal
import numpy as np


def loadMovement(fileName):
    csvFile = csv.reader(open(fileName, 'r'), delimiter=',')
    result = []
    for row in csvFile:
        if row[1] != 'Target':
            result.append([])
            for val in row[:5]:
                if len(result[-1]) == 4 and result[-1][3] == 'yes':
                    result[-1].append(float(val))
                else:
                    result[-1].append(val)
    return result


def initializeGLS():
    fileName = FOLDER_IN+'Movement_Validation/EMBD_Movement_Validation_EMBD0000_GLS_corrected.csv'
    # fileName = FOLDER_IN + 'Movement_Validation/EMBD_Movement_Validation_RNAi_GLS.csv'
    moveList = loadMovement(fileName)
    correct, fPos, fNeg, noM = 0, 0, 0, 0
    diff = []
    tMoveAll = []
    ls = []
    for embInf in moveList:
        if fileName.split('_')[-3] == 'EMBD0000':
            folder = FOLDER_IN + 'cropped/{0}/GLS/{1}/{2}/'.format(embInf[1], embInf[0], embInf[2])
        else:
            folder = FOLDER_IN + 'cropped/{0}/GLS/{1}/'.format(embInf[1], embInf[2])
        embs = loadEmbs(folder)
        embs = refreshEmbsMulti(embs)
        if len(embs) > 0:
            emb = embs[0]
            #             tMove, success = moveDetect.getMovementGS(emb)
            #             if success and not emb.movement: print('!!!!!!!!!!!!!!problem emb', emb.label, emb.RNAi, emb.date)
            tMove = emb.tMove
            #             if emb.movement: print('resids Late=',emb.params['residsLateY'])
            success = emb.movement
            if success and embInf[3] == 'yes':
                printLog('Movement detected for {0}/{1}/{2} difference {3}'.format(embInf[1], embInf[0], embInf[2],
                                                                                   embInf[4] - tMove))
                diff.append(embInf[4] - tMove)
                ls.append(emb.label)
                correct += 1
                tMoveAll.append(emb.params['tMove'])
            if not success and embInf[3] == 'yes':
                printLog('MISSED movement for {0}/{1}/{2} at {3}'.format(embInf[1], embInf[0], embInf[2], embInf[4]))
                print('tail-b=', emb.params['Ytail'] - emb.params['bY'], 'resids=', emb.params['residsLateY'],
                      emb.params['YtailErr'])
                fNeg += 1
            if success and embInf[3] != 'yes':
                printLog('FALSE movement for {0}/{1}/{2} at {3}'.format(embInf[1], embInf[0], embInf[2], tMove))
                print('tail-b=', emb.params['Ytail'] - emb.params['bY'], 'resids=', emb.params['residsLateY'],
                      emb.params['YtailErr'])
                fPos += 1
            if not success and embInf[3] != 'yes':
                printLog('No movement for {0}/{1}/{2}'.format(embInf[1], embInf[0], embInf[2]))
                correct += 1
                noM += 1
    print('movement find precision {0}\n{1} correct\n{2} missed\n{3} false positive\n{4} not moving'.format(
        getStrVal(np.mean(diff), np.std(diff)), correct, fNeg, fPos, noM))
    print('tMove from mG', np.mean(tMoveAll), np.std(tMoveAll))


def initializeMS():
    fileName = FOLDER_IN + 'Movement_Validation/EMBD_Movement_Validation_EMBD0000_MS.csv'
    # fileName = FOLDER_IN+'Movement_Validation/EMBD_Movement_Validation_RNAi_MS.csv'
    moveList = loadMovement(fileName)
    correct, fPos, fNeg, noM = 0, 0, 0, 0
    diff = []
    tMoveAll = []
    ls = []
    for embInf in moveList[:]:
        if fileName.split('_')[-2] == 'EMBD0000':
            folder = FOLDER_IN + 'cropped/{0}/MS/{1}/{2}/'.format(embInf[1], embInf[0], embInf[2])
        else:
            folder = FOLDER_IN + 'cropped/{0}/MS/{1}/'.format(embInf[1], embInf[2])
        embs = loadEmbs(folder)
        embs = refreshEmbsMulti(embs)
        if len(embs) > 0:
            emb = embs[0]
            #             if emb.loaded: emb.refresh()
            tMove = emb.tMove
            #             if emb.movement: print('resids Late=',emb.params['residsLateY'])
            success = emb.movement
            #             if emb.movement: print('resids Late=',emb.params['residsLateY'])
            if success and embInf[3] == 'yes':
                printLog('Movement detected for {0}/{1}/{2} difference {3}'.format(embInf[1], embInf[0], embInf[2],
                                                                                   embInf[4] - tMove))
                diff.append(embInf[4] - tMove)
                ls.append(emb.label)
                correct += 1
                tMoveAll.append(emb.params['tMove'])
            elif not success and embInf[3] == 'yes':
                printLog('MISSED movement for {0}/{1}/{2} at {3}'.format(embInf[1], embInf[0], embInf[2], embInf[4]))
                print('tail-b=', emb.params['Ytail'] - emb.params['bY'], 'resids=', emb.params['residsLateY'])
                fNeg += 1
            elif success and embInf[3] != 'yes':
                printLog('FALSE movement for {0}/{1}/{2} at {3}'.format(embInf[1], embInf[0], embInf[2], tMove))
                fPos += 1
            elif not success and embInf[3] != 'yes':
                printLog('No movement for {0}/{1}/{2}'.format(embInf[1], embInf[0], embInf[2]))
                correct += 1
                noM += 1
    printLog('movement find precision {0}\n{1} correct\n{2} missed\n{3} false positive\n {4} not moving'.format(
        getStrVal(np.mean(diff), np.std(diff)), correct, fNeg, fPos, noM))
    print('tMove from mG', np.mean(tMoveAll), np.std(tMoveAll))
    print('total={0}, counted={1}'.format(len(moveList), correct + fPos + fNeg))


if __name__ == '__main__':
    # print('______________MS________________')
    # initializeMS()
    print('______________GLS________________')
    initializeGLS()
