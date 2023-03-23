"""
Created on May 19, 2016

@author: ODLab

Collection of movement detection functions
"""

import numpy as np
from myFigure import myFigure
from myFunc import find_consec_inds
from scipy import ndimage
import cv2
from varLookup import GLS_MOVE_PRECISION, MS_MOVE_PRECISION, FOLDER_IN


def getMovementGS(emb):
    """uses spots data to determine"""
    return getMovementGS1(emb)


def getMovementMS(emb):
    """uses match template approach to find movement for MS strain- finds dif between pixels in consecutive images and sums residuals.
    Checks for threshold of 0.7 after timepoint 10. Must move for >=3 frames to be called movement. Must continue moving for
    >6 frames to be scored as successful movement (if not, arrests at 3-fold)."""
    return getMovement2(emb)


def getDead(emb):
    return getDead1(emb)


def getMovement1(emb):
    colors = ['g', 'r', 'k']
    fig = myFigure()
    for ch in [1]:
        ims = []
        for i in range(31):
            im = np.float32(emb.image.getMaxProject(i, ch))
            a, b = im.shape
            ims.append(im[a / 8:7 * a / 8, b / 8:7 * b / 8])  # shaves edges to eliminate debris or neighboring embryos
        imCorr = []
        for i in range(30):
            xt, yt = cv2.phaseCorrelate(ims[i], ims[i + 1])
            imCorr.append(np.sqrt(xt ** 2 + yt ** 2))
        imCorr = ndimage.filters.median_filter(imCorr, size=5)
        # fig.plot(np.arange(1, 31) + 1, imCorr, label='{0}_{1}'.format(emb.label, emb.date))
    # fig.show()
    maxV = np.max(imCorr[10:])
    print('1 max val = ', maxV)
    try:
        tMove = np.arange(1, 31)[10:][np.where(imCorr[10:] > 25)[0][0]]
        print('{0}_{1} moves at t={2}'.format(emb.label, emb.date, tMove))
        return tMove, True
    except:
        print('{0}_{1} No movement'.format(emb.label, emb.date))
        return np.nan, False


def getMovement2(emb):
    """ Uses match template approach to find movement for MS strain- finds dif between pixels in consecutive images and sums residuals.
        Checks for threshold of 0.7 after timepoint 10. Must move for >=3 frames to be called movement. If stops for 3 consecutive frames, after onset of movement, movement is called
        false (3-fold), unless movement is observed for 3 additional timepoints after the stoppage point (movement called true). Also checks if embryo swishes out of view (often due to hatching)
        -this is detected based on drop in total intensity levels, relative to max intensity. If swished after onset of movement, movement is called true
        :arg: emb
        :returns:
        tMove-  the first timepoint of movement
        movement- true/false -- indicates whether the embryo moves and continues to move throughout the data series
        """

    show = False
    #     from totalintensityFit import fitSigm
    #     y = emb.red
    #     x = np.arange(y.size)
    #     popt, perr, residsEarly, residsLate, moveT, fig = fitSigm(x, y, y.size)
    maxG = 12225
    #     maxR = 19025
    maxR = 30000
    #     maxV = [maxG, maxR]
    colors = ['g', 'r', 'k']
    if show: fig = myFigure()
    #     im0=True

    # avg_median_cont = 463032304.425
    # std_median_cont = 81186307.9047

    for ch in [1]:
        ims = [] #max projections of the red channel for each timepoint
        for i in range(31):
            im = np.float32(emb.image.getMaxProject(i, ch))
            a, b = im.shape
            ims.append(im)
        totInt = emb.tInt[1]
        # print totInt
        # print "max = {0}".format(np.max(totInt))
        # print totInt/np.max(totInt)
        # print totInt/avg_median_cont
        # print np.median(totInt)/avg_median_cont
        max_intens = np.max(totInt) #intensity measures are used to check if embryo is swished out of view
        max_ind = np.where(totInt == np.max(totInt))[0] #index where totInt is max
        ratio_intens = totInt/max_intens #list of ratios
        thr = 0.5 #the ratio_intens is always near 1 after the max unless the emb has swished out. We use thresh of 0.5 for ratio intensity max.

        tf_ratio_intens = (ratio_intens <= thr)  # converts to true/false, where true means that the ratio of intensity/med int is below specified threshold
        ratio_inds = np.where(tf_ratio_intens == True) #gets index positions where ratio is true (below specified thresh)
        swish, swish_ind = find_consec_inds(ratio_inds[0], 3) #checks for 3 consecutive true evaluations- true/false indicates swish out of view- swish_ind is a list of the first index of a series of consecutive indices
        if len(swish_ind) > 0:
            swish_time = ratio_inds[0][swish_ind] #gives list of the first index position (of 3 consecutive) that are below specified threshold
            if len(swish_time) > 0:
                swishFlag = False #set to false as a default, since swish must occur after median
                for s in swish_time:
                    if s > max_ind: #checks if swish occurs after max index
                        swishFlag = True
                        break
                # print "swish is {0} at {1}".format(swishFlag, s)
        else:
            swishFlag = False
            # print "swish is {0}".format(swishFlag)


        imCorr = [] #a list of match template measures between each consecutive frame
        for i in range(30):
            mtCoeff = cv2.matchTemplate(ims[i], ims[i + 1], 0)
            imCorr.append(mtCoeff[0])
        imCorr = np.array(imCorr)
        norm = []
        for im in ims[-5:]:
            norm.append(np.sum(im ** 2))
        norm = np.mean(norm)
        imCorr = imCorr / norm
        #         imCorr = ndimage.filters.median_filter(imCorr, size=4)
        dy = ((np.roll(totInt, -1) - totInt) / np.median(totInt[-10:]))[:-1]
        imCorr[np.where(dy < -0.5)] = 0.
        if show: fig.plot(np.arange(1, 31) + 1, imCorr, label='{0}_{1}'.format(emb.label, emb.date))

    # try:
    if show:
        fig.ylim((0, 1.8))
        fig.show()

    if np.sum(imCorr[10:] > 0.7) < 3: #checks to see if embryo moves for at least 3 timepoints. If not, returns nan for tMove and False for movement
        # print('{0}_{1} No movement'.format(emb.label, emb.date))
        return np.nan, False
    else:
        inds = np.where(imCorr[10:] > 0.7)[0] #gets all indices after the 10th timepoint where imCorr is above the 0.7 threshold
        tMove = np.nan
        for i in range(inds.size-1): #checks to be sure that tMove is found accurately, by requiring movement to continue within 3 timepoints (avoids issue with elongation being falsely called as movement)
            if inds[i+1] - inds[i] < 3:
                tMove = np.arange(1, 31)[10:][inds[i]]
                break
        if tMove == np.nan: #if movement does not continue within three points of movement at any point, tMove is called nan and movement is returned False
            # print('{0}_{1} No movement'.format(emb.label, emb.date))
            return np.nan, False
        t_move_stop = np.arange(1, 31)[tMove:][np.where(imCorr[tMove:] <= 0.7)[0]]  #checks for indices after tMove where the embryo is stopped
        if len(t_move_stop) >= 3: #checks to see if after movement starts, if it ever stops (must have at least 3 stopped points to be considered stopped)
            movement = True
            n = 3
            consec, ind_list = find_consec_inds(t_move_stop, n) #checks to see if 3 consecutive indices are present in t_move_stop- 3 stopped timepoints in a row, after the onset of movement
            if consec:
                movement = False #3 consecutive indices indicates that the embryo has stopped moving (3-fold arrest)
                stop_ind = ind_list[0]
                stop_time = t_move_stop[ind_list[0]] #this reports the first timepoint (of n consecutive timepoints) where the embryo is considered stopped moving
                t_move_again = np.arange(1, 31)[stop_time+n:][np.where(imCorr[stop_time+n:] >= 0.7)[0]]
                if len(t_move_again) >= 3: #checks if there are 3 or more moving timepoints after the stoppage point
                    movement = True
                if swishFlag:
                    movement = True  # swished out of view


        else:
            movement = True #moves and continues to move
        # print('{0}_{1} moves at t={2} and movement is {3}'.format(emb.label, emb.date, tMove, movement))
        return tMove + MS_MOVE_PRECISION, movement

        # tMove = np.arange(1, 31)[10:][np.where(imCorr[10:] > 0.7)[0][0]] #determines the timepoint for movement
        # tMoveStop = np.arange(1, 31)[10:][np.where(imCorr[10:] > 0.7)[0][-1]]
    #         tMoveStop = np.arange(1,31)[tMove:][np.where(imCorr[tMove:]<=0.7)[0]]
    #         success = True
    #         if len(tMoveStop)>0:
    #             for k, g in groupby(enumerate(tMoveStop), lambda (i, x): i-x):
    #                 res = map(itemgetter(1), g)
    #                 if len(res)>=3:
    #                     success = False
    #                     break
    #             if success: print('{0}_{1} moves at t={2}'.format(emb.label, emb.date, tMove))
    #             else: print('{0}_{1} moves at t={2} and stops at t={3}'.format(emb.label, emb.date, tMove, tMoveStop[0]))
    #             return tMove, success
    #         else:

    # if tMoveStop - tMove > 6 or tMoveStop >= 30:
    #     print('{0}_{1} moves at t={2}'.format(emb.label, emb.date, tMove))
    #     return tMove + MS_MOVE_PRECISION, True
    # elif 6 >= tMoveStop - tMove > 3:
    #     print('{0}_{1} moves at t={2} and stops at t={3}'.format(emb.label, emb.date, tMove, tMoveStop))
    #     return tMove + MS_MOVE_PRECISION, False
    # else:
    #     print('{0}_{1} No movement'.format(emb.label, emb.date))
    #     return np.nan, False
    #
        # except:
    #     print('{0}_{1} No movement'.format(emb.label, emb.date))
    #     return np.nan, False


def getMovementGS1(emb):
    show = False
    th = 20
    y = emb.yellow
    if np.max(y) == 0: y = emb.green
    if np.max(y) == 0: y = emb.red
    if np.max(y) > 0:
        x = np.arange(y.size)
        #     dy = np.gradient(y)
        dy = np.abs(np.roll(y, -1) - y)
        dy[-1] = 100
        move = False
        #     dy = np.roll(y, -1)-y
        if show:
            fig = myFigure()
            fig2 = myFigure()
            fig.plot(x, dy)
            fig.title('product')
            fig2.plot(x, y)
            fig2.scatter(x, y)
        inds = np.where((y[:] > 50) * (y[:] < 200))[0]
        plateauN = 2
        plateau = 0
        position = 0
        flag = False
        if inds.size > 3:
            for i in inds[:-3]:
                if all(dy[i:i + plateauN] < th):
                    position = i + plateauN
                    flag = True
                    plateau = np.max(y[i:i + plateauN])
                elif all((y[i + plateauN:i + plateauN + 2] - plateau) > th):
                    if flag:
                        move = True
                        break
        peak = False
        if position > 0:
            endInd = position
        else:
            maxV = 0
            for val in y:
                if val > maxV: maxV = val
                if val < 0.8 * maxV and maxV > 0.8 * np.max(y):
                    peak = True
                    break
            if peak:
                endInd = np.argmin(np.abs(y - maxV))
            else:
                endInd = y.size - 1
        # if endInd<25:
        if show:
            print('getMovementGS1:', endInd + 1, move, plateau)
            fig.show()
        if move:
            return endInd + 1 + GLS_MOVE_PRECISION, move  # +1 because movement is in the next slide after the plateau
        else:
            return 31, move
    # else: return y.size, False
    else:
        return 31, False


def getDead1(emb):
    if emb.image is None: emb.loadImages()
    if emb.params['maxR'] > 3e8 or emb.params['maxG'] > 3e8: return False
    ims = np.float32(emb.image.images[:, 8, 2])
    imCorr = []
    for i in range(30):
        mtCoeff = cv2.matchTemplate(ims[i], ims[i + 1], 1)
        imCorr.append(mtCoeff[0])
    imCorr = np.array(imCorr)
    #     fig = myFigure()
    #     fig.plot(range(imCorr.size), imCorr)
    #     if np.mean(imCorr[10:])<0.02: print('DEAD embryo {0}'.format(emb.folder))
    #     print(np.mean(imCorr[10:]))
    return np.mean(imCorr[5:10]) < 0.02


if __name__ == '__main__':
    from emb_handler import loadEmbs

    folder = 'Emb6/'

    folders = [FOLDER_IN + folder]
    for folder in folders:
        embs = loadEmbs(folder)
        for emb in embs:
            emb.refresh()
            # emb.loadImages()
            tMove, movement = getMovement2(emb)
            print tMove, movement
            # print(emb.dead)
            #         plt.show()
