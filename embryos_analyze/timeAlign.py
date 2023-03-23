"""
Created on Jan 17, 2017

@author: renat

Collection of functions related to time alignment
"""

import numpy as np
from lmfit import Parameters, minimize, report_fit, Parameter
from scipy.interpolate import interp1d
from varLookup import FOLDER_IN, printLog, AG_AVG_MS, SG_AVG_MS, AR_AVG_MS, SR_AVG_MS, T0G_M_T0R_AVG_MS,\
    T0R_M_T0G_AVG_GLS, T0Y_M_T0R_AVG_GLS, AR_AVG_GLS, AY_AVG_GLS, SR_AVG_GLS, SY_AVG_GLS
from myMath import getSubIndex, getSig
from myFigure import myFigure


def getAlignParams(emb, curves, show=False, shiftVary=True, scaleVary=True, tScaleVary=False, verb=False,
                   firstPass=False):
    """
    Fits multiple embryo curves to the average curves allowing for shift and scale in time (one for all curves),
    y scaling (multiple) and a point of deviation (multiple)
    INPUT:
    :param emb: Embryo object
    :param curves: ist of curve names
    :param show: show alignment graph
    :param shiftVary: allow shift in x
    :param scaleVary: allow scale in y
    :param tScaleVary: allow scale in x
    :param verb: print progress
    :param firstPass: use boundaries (min/max values) for the first pass
    RETURN:
    shift and scale for x values
    list of indices at point of deviation, one per each curve
    list of y scalings and its standard deviations, one per each curve
    """
    scaleSpots = False  # do not allow scaling of spots in y
    if not shiftVary and not scaleVary and not tScaleVary:  # check that something could be varied otherwise raise error
        printLog('getAlignParams: can not vary shift and scaling during fitting to average curves')
        raise
    nCurves = len(curves)
    xRefList = []  # time for reference (control average) curves
    yRefList = []  # values for reference (control average) curves
    yRefErrList = []  # errorbars for reference (control average) curves
    fyRefList = []  # method of interpolation of reference values between time points
    fyRefErrList = []  # method of interpolation of reference errors between time points
    xIniList = []  # unaligned (initial) time points of the embryo curves
    yIniList = []  # unaligned (initial) values points of the embryo curves
    nMaxList = []  # maximum number of points in a curve that could be used (less then 31, because some curves start
    # later in time)
    yRefMeanList = []  # average value of the reference curve (used for normalization to 1)
    startPointList = []  # list of start points for each curve
    nLambda = emb.nLambda  # lambda value set in varLookup and stored in each embryo
    tauScale = emb.tauScale  # tau value (NOT USED ANYMORE, was used to enforce t scaling to 1)
    for curveName in curves:  # check that for all embryo curve a reference average curve exists
        if curveName not in emb.avgCurves:
            printLog('wrong curve name')
            raise
        else:  # read the average curves and populate reference lists
            xRef, yRef, yRefErr = emb.avgCurves[curveName]
            # yRefErr = np.ones_like(yRefErr)  # Use this for avg calculator MS curves only
            xRefList.append(xRef)
            yRefErrList.append(yRefErr)
            yRefList.append(yRef)
            fyRefErrList.append(interp1d(xRef, yRefErr, 'linear', fill_value=np.nan))
            fyRefList.append(interp1d(xRef, yRef, 'linear', fill_value=np.nan))
            yRefMeanList.append(np.mean(yRef))
        xIni, yIni = emb.getCurve(curveName)  # read embryo curves
        xIniList.append(xIni)
        yIniList.append(yIni)
        if 31 - sum(np.isnan(yIni)) > 0:  # check how many points are not nan, this is the maximum number of point
            # that we can use
            nMaxList.append(31 - sum(np.isnan(yIni)))
        else:
            nMaxList.append(1)
        #         nMaxList.append(31.)
        startPointList.append(emb.startPoint[curveName])
    xIniList = np.array(xIniList)
    yIniList = np.array(yIniList)
    nMaxList = np.array(nMaxList)
    yRefMeanList = np.array(yRefMeanList)
    startPointList = np.array(startPointList)
    nMin = 3  # minimum number of point for an embryo curve to overlap with the reference curve

    def getParams(nListTmp, m0=None):
        """
        fits parameters (time shift [xShift], time scale [xScale], and y scale [yScale]) for fixed cutoffs)
        :param nListTmp: list of end points to use for each of the curves
        :param m0: initial values to use for the parameters (we use always None)
        :return: lmfit minimizer (an object that stores all fitted parameters and residuals)
        """

        def J(pars):
            """
            function that calculates residuals between reference curves and embryo curves for given set of parameters
            :param pars: parameters (xShift, xScale, set yScales for each curve, and set of Ns (points of diviation)
            :return: residuals
            """
            xsh = pars['xShift'].value  # read values from pars object (lmfit Parameters class)
            xScale = pars['xScale'].value
            yScale = np.array([pars['yScale{0}'.format(i)].value for i in range(nCurves)])
            N = np.array([int(pars['N{0}'.format(i)].value) for i in range(nCurves)])
            Jcurve = 31. / nMaxList  # weighting for each curve contribution. Makes sure that curves with more points
            # do not contribute more to alignment.
            x = (xIniList - xsh) * xScale  # scale and shift initial time with time shift and time scale
            y = (yIniList.T * yScale).T  # scale curve values with yScale
            yre = np.ones_like(x) * np.inf  # create numpy arrays for reference values to avoid for loops for speed
            yr = np.ones_like(x) * np.nan
            useMatrix = np.ones_like(x)  # 0/1 array that defines if a value is used for residuals calculation
            wUseM = np.ones_like(x)  # array of weights for each residual
            for i in range(
                    nCurves):  # loop over all curves
                # (use loop because interpolation doesn't work for regions outside defined zone and need truncations)
                inds = getSubIndex(x[i], xRefList[i][0], xRefList[i][-1])  # find indices of embryo curve points
                # that overlap with the reference curves (do not use points that have no reference curve)
                yr[i][inds] = fyRefList[i](x[i][inds])  # interpolate values and errors from reference curves to
                # correspond to time points of embryo curve
                yre[i][inds] = fyRefErrList[i](x[i][inds])
                useMatrix[i][N[i]:] = np.ones(useMatrix[i].size - N[i]) * np.nan
                wUseM[i][N[i]:] = np.ones(wUseM[i].size - N[i]) * np.nan  # make all values after after cutoff nans
                if startPointList[i] > 0:
                    wUseM[i][:startPointList[i]] = np.ones(startPointList[i]) * np.nan  # make all values before start
                    # point nans
            #             w = 1./(yre/yr)**2*wUseM
            w = 1. / (yre) ** 2 * wUseM  # use error bars to calculate weights for the residuals
            w = (w.T / np.nansum(w, axis=1)).T  # normalize weights to sum to one
            #             w=np.ones_like(w)
            Jcurve = Jcurve * np.nansum(wUseM,
                                        axis=1)  # multiply by number of used points because we use normalized weights
            # that sum to one. If less points are used the individual weights for each point are higher,
            # because all weight sum to one. To correct for this effect we multiply by the number of points used
            val = ((w * ((yr - y).T / yRefMeanList).T ** 2).T * Jcurve).T  # calculate weighted residuals
            #             val = ((w*((yr-y)/yr)**2).T*Jcurve).T#R2
            val = val * useMatrix  # only specifica values may be used (before point of truncation)
            if np.sum(np.isnan(val)) >= val.size - nMin:
                return 1e8  # check that there are at least 3 values that used, otherwise return BIG number
            return val[~np.isnan(val)].ravel() + tauScale * (1. / emb.tScale - xScale) ** 2  # we add a term to
            # all residuals that forces xScale to be closer to 1. But we do not use it anymore, because tauScale = 0.

        p = Parameters()
        spotsCurves = []
        if not scaleSpots:  # check if spots are allowed to be scaled in y
            spotsCurves = ['spotsG', 'spotsR', 'spotsY']
        np.random.seed(1)
        if firstPass:  # use different adjustment criteria for the first pass and the second
            p['xShift'] = Parameter(value=0., min=-3, max=3, vary=shiftVary)
            p['xScale'] = Parameter(value=1., min=0.5 / emb.tScale, max=1.3 / emb.tScale, vary=tScaleVary)
        else:
            p['xShift'] = Parameter(value=0., min=-1, max=1, vary=shiftVary)
            p['xScale'] = Parameter(value=1., min=0.8, max=1.2, vary=tScaleVary)
        for i in range(len(curves)):
            p['N{0}'.format(i)] = Parameter(value=nListTmp[i], vary=False)  # the cutoff point is not varied
            p['yScale{0}'.format(i)] = Parameter(value=1., min=0.5, max=2.,
                                                 vary=scaleVary and curves[i] not in spotsCurves)  # if spots scaling
            # is not allowed (scaleSpots = False) then do not scale spots curves in y
        minzer = minimize(J, params=p, method='cobyla', options={'maxiter': 50})  # fit parameters
        # by minimizing residuals
        return minzer

    def getCosts(nListM, m0=None):
        """
        calculates total cost function for given cutoffs
        :param nListM: 2d array of multiple sets of list of cutoffs for the set of curves.
        1s dimension is going through various sets of catofs, second dimension is going through cutoff for each curve.
        :param m0: initial values for fitting parameters
        :return: cost and lmfit minimizer object that hold fitted parameters and residuals
        """
        costs = np.ones(nListM.shape[0]) * 1.e10  # set initial cost high
        minzer = None
        for j in range(nListM.shape[0]):  # loop through various sets of cutoffs to determine the one with best cost
            if all(nListM[j] >= startPointList + nMin) and all(nListM[j] <= 31):  # check that all cutoffs are in
                # reasonable range (at least nMin (3) points and less or equal to 31.
                minzer = getParams(nListM[j])  # fit parameters
                costs[j] = np.sum(minzer.residual) + sum(
                    [nLambda * 31. / nMaxList[i] * (1. - 1. * (nListM[j][i] - startPointList[i]) / nMaxList[i]) ** 2 for
                     i in range(nCurves)])  # calculate cost for all curves
        return costs, minzer

    stepRange = 5  # the range to change the cutoffs in each minimization iteration.
    nList = np.ones(nCurves) * (30 - stepRange)  # start with cutoffs for all curves as 30 - steprRange (=25)
    nListNew = np.ones(nCurves)
    nList[np.where(nList <= startPointList + nMin)] = (startPointList + nMin + stepRange)[
        np.where(nList <= startPointList + nMin)]  # set values in nList above startPoint + nMin and below 31
    nList[np.where(nList > 31)] = 31 * np.ones_like(nList)[np.where(nList > 31)]
    nList0 = np.ones(nCurves) * 0
    nList1 = np.ones(nCurves) * 0
    nIter = 0
    m0 = None
    while any(nList != nList0) and any(nList != nList1) and nIter < 5 and stepRange > 0:
        # stop if the values in nList don't change or became as in previous step (optimum achieved or oscillations)
        # or made 5 iterations
        nList1 = nList0  # nList in previous step
        nList0 = np.copy(nList)  # nList in this step
        steps = np.ones(nCurves)
        for j in range(nCurves):  # vary cutoff for curve j only keeping the rest of the cutoffs fixed.
            # Determine the best cutoff that minimizes cost function.
            if startPointList[j] + nMin < nList[j] <= 31:
                nListVary = np.tile(nList, stepRange * 2 + 1).reshape(stepRange * 2 + 1, nList.size)  # make a matrix
                # of nLists in each row and 2*stepRange +1 rows. This matrix will the set of cutoffs to consider.
                nListVary[:, j] += np.arange(-stepRange, stepRange + 1, 1)  # modify cutoff values of curve j in each
                # set to change up and down in range of stepRange.
                costs, minzer = getCosts(nListVary, m0)  # calculate costs for each cutoff set and their residuals
                steps[j] = np.abs(np.argmin(costs) - int(costs.size / 2))  # determine cutoff change for curve j that
                # produces minimual cost
                nListNew[j] = nListVary[np.argmin(costs)][j]  # assign the best cutoff to be used in next iteration
                m0 = minzer  # used to use as the starting point for next iteration, but not anymore.
        stepRange = min(5, int(np.max(steps)))
        nList = np.copy(nListNew)  # continue next iteration with new cutoff values for each curve
        nIter += 1
    minzer = getParams(nList)  # read out all optimized parameters
    xShift = minzer.params['xShift'].value
    xScale = minzer.params['xScale'].value

    yScaleList, yScaleListErr = [], []
    for i in range(nCurves):
        yScaleList.append(minzer.params['yScale{0}'.format(i)].value)
        yScaleListErr.append(minzer.params['yScale{0}'.format(i)].stderr)
    if show:  # show figure if desired
        print('successfully minimized = ', minzer.success)
        print(report_fit(minzer.params))
        print(nList)
        for i in range(len(curves)):
            N = nList[i]
            print('initial t0 = {0}, refined t0 = {1}'.format(emb.t0, emb.t0 + xShift))
            fig = myFigure()
            fig.errorbar(xRefList[i][::10], yRefList[i][::10], yRefErrList[i][::10], color='k')
            #             fig.plot(xRef,yRef,color='k')
            fig.plot(xIniList[i], yIniList[i], color='r', label='initial')
            fig.plot(xScale * (xIniList[i] - xShift), yIniList[i] * yScaleList[i], color='g', label='aligned')
            fig.plot(xScale * (xIniList[i][:N] - xShift), yIniList[i][:N] * yScaleList[i], 'k--')
            fig.legend(2)
            fig.title(curves[i])
        fig.show()
    if verb: printLog('timeAlign.getAlignParams: {cn} nList={n}'.format(cn=curves, n=nList))
    return xShift, nList, yScaleList, xScale


def alignIntens(emb):
    """
    preliminary alignment using fits to total intensity
    :param emb: MSEmbryo object
    :return: None
    """
    t0 = []  # set of time shifts
    w = []  # confidence weights
    ts = []  # set of time scales
    if emb.checkGoodParam('mG') and emb.checkGoodParam('sG') and emb.checkGoodParam('aG'):  # if green signal is good
        t0.append(emb.params['mG'] - 2. * emb.params['sG'])  # t0 at the intercept of tangent line at timepoint m to
        # green signal with baseline
        ts.append(emb.params['aG'] / AG_AVG_MS * SG_AVG_MS / emb.params['sG'])  # tscale is estimated from slope
        # relative to control. (Note that there is correction to s values arising from differences in plateau values)
        w.append(np.exp(-(emb.params['aG'] / AG_AVG_MS - 1) ** 2))  # confidence is gaussian with optimum at average.
    if emb.checkGoodParam('mR') and emb.checkGoodParam('sR') and emb.checkGoodParam('aR'):  # same as for green but
        # with shift by the average distance between red and green
        t0.append(emb.params['mR'] - 2. * emb.params['sR'] + T0G_M_T0R_AVG_MS)
        ts.append(emb.params['aR'] / AR_AVG_MS * SR_AVG_MS / emb.params['sR'])
        w.append(np.exp(-(emb.params['aR'] / AR_AVG_MS - 1) ** 2))
    w = np.array(w)
    t0 = np.array(t0)
    ts = np.array(ts)
    w = w / sum(w)
    if t0.size > 0:
        emb.t0 = t0[0]  # always use the first t0 (determined by green intensity if green is good, otherwise by red)
        #         emb.t0 = t0[-1]#sum(w*t0)
        if 6.5e8 > np.nanmax(emb.getCurve('tIntR')[1]) > 4.5e8:  # check if red intensity reaches this value
            emb.tScale = np.mean((ts[-1], 1))  # use half way scale to what is predicted by red.
            # (only determine direction of scale, magnitude might overshoot, so we reduce impact this way)


def alignHead(emb):
    """
    get t0 and tscale by aligning head curve only. We don't use it anymore.
    :param emb: MSEmbryo object
    :return:
    """
    t0 = []
    if emb.checkGoodParam('mSigHead') and emb.checkGoodParam('sSigHead') and emb.checkGoodParam('aSigHead'):
        t0.append(emb.params['mSigHead'] - 2. * emb.params['sSigHead'])
    if len(t0) > 0:
        emb.t0 = t0[0]
    else:
        alignIntens(emb)


def alignSpots(emb):
    """
    Preliminary alignment and scaling of embryo using sigmoidal fits to spots curves
    :param emb: GSEmbryo object
    :return:
    """
    t0, ts = [], []
    # find t0 as the intercept between tangent line at mR to red spots and 0.
    # tscale is calculated from s valuesof red and yellow.
    if not np.isnan(emb.params['mG']) and not np.isnan(emb.params['sG']) and emb.checkGoodParam('aG'):
        t0.append(emb.params['mG'] - 2. * emb.params['sG'] + T0R_M_T0G_AVG_GLS)
    if not np.isnan(emb.params['mR']) and not np.isnan(emb.params['sR']) and emb.checkGoodParam('aR'):
        t0.append(emb.params['mR'] - 2. * emb.params['sR'])
    if not np.isnan(emb.params['mY']) and not np.isnan(emb.params['sY']) and emb.checkGoodParam('aY'):
        t0.append(emb.params['mY'] - 2. * emb.params['sY'] - T0Y_M_T0R_AVG_GLS)
    if not np.isnan(emb.params['mG']) and not np.isnan(emb.params['sG']) and emb.checkGoodParam('aG') and \
            not np.isnan(emb.params['mR']) and not np.isnan(emb.params['sR']) and emb.checkGoodParam('aR') and \
            not np.isnan(emb.params['mY']) and not np.isnan(emb.params['sY']) and emb.checkGoodParam('aY'):
        ts.append(emb.params['aR'] / AR_AVG_GLS * SR_AVG_GLS / emb.params['sR'])
        ts.append(emb.params['aY'] / AY_AVG_GLS * SY_AVG_GLS / emb.params['sY'])
    if len(t0) > 0:  # check that there are any t0 values
        ts = np.array(ts)
        if ts.size > 0:
            emb.tScale = ts[np.argmin((1 - ts) ** 2)]  # use tscale closest to 1
        if len(t0) > 1:
            emb.t0 = t0[1]  # use t0 estimated from red (if red exists) or yellow
        elif len(t0) == 1:
            emb.t0 = t0[0]  # otherwise use the one from green


def alignTimeMS(emb):
    """
    time alignment of MS embryos
    :param emb: MSEmbryo object
    :return: time aligned MSEmbryo object
    """
    emb.tScale = 1.  # make sure to start fresh
    emb.t0 = 0.
    alignIntens(emb)  # preliminary alignment by intensity
    emb.tScale = 1.  # don't use time scale from intensity alignment (not reliable)
    if emb.curvesLoaded:  # if average curves are loaded
        for i in range(1):  # perform only one iteration of alignment
            curves = ['lengthR', 'MoI0R']  # use lengthR and MoI0R
            if np.nanmax(emb.getCurve('tIntG')[1]) > 3.5e8:  # if green intensity is good, use it
                curves += ['tIntG']
            if np.nanmax(emb.getCurve('tIntR')[1]) > 4.5e8:  # if red intensity is good, use it
                curves += ['tIntR']
            x0, nDev, yScale, xScale = getAlignParams(emb, curves, show=False, shiftVary=True, scaleVary=False,
                                                      tScaleVary=False, verb=True, firstPass=True)  # align using these
            # curves. Don't show alignment (show), allow time shift (shiftVary), do not scale in time (tScaleVary)
            # or y (scaleVary), print progress (verb), firstPass (first iteration).
            emb.t0 += x0
            emb.tScale *= xScale  # in this case xScale is 1 because we didn't allow it to vary in previous step

            curves = ['MoI0R', 'lengthR']  # in this situation use the same curves
            # if np.nanmax(emb.headInt) > 1.: curves += ['headInt']
            if np.nanmax(emb.getCurve('tIntR')[1]) > 4.5e8:
                curves += ['tIntR']
            if np.nanmax(emb.getCurve('tIntG')[1]) > 3.5e8:
                curves += ['tIntG']
            x0, nDev, yScale, xScale = getAlignParams(emb, curves, show=False, shiftVary=True, scaleVary=True,
                                                      tScaleVary=True, verb=True, firstPass=True)  # align using these
            # curves. Don't show alignment (show), allow time shift (shiftVary), scale in time (tScaleVary)
            # and y (scaleVary), print progress (verb), firstPass (first iteration).
            emb.t0 += x0  # update t0 and scale with new values
            emb.tScale *= xScale

    else:
        print('alignTimeMS: {0} curves not loaded'.format(emb.label))
    # used for all embryos processing prior to 03-07-18
    # emb.t0 = np.round(emb.t0, 1)
    # emb.tScale = np.round(emb.tScale, 2)
    # NOTE: we did not rerun all embryos to reflect the change, because the change is marginal for results
    emb.t0 = float(np.round(emb.t0, getSig(emb.t0) + 5))  # changed from round to float to stabilize borderline values
    # in t0, i.e. 0.49 vs 0.5
    emb.tScale = float(np.round(emb.tScale, getSig(emb.tScale) + 5))
    return emb


def alignTimeGLS(emb):
    """
    time alignment of GS embryos
    :param emb: GSEmbryo object
    :return: time aligned GSEmbryo object
    """
    emb.tScale = 1.  # make sure to start fresh
    emb.t0 = 0.
    alignSpots(emb)  # preliminary alignment based on spots intercept
    if emb.curvesLoaded:
        for i in range(2):  # do two iterations
            curves = ['spotsR']  # only use red spots
            x0, nDev, yScale, xScale = getAlignParams(emb, curves, show=False, shiftVary=True, scaleVary=False,
                                                      tScaleVary=False, verb=True, firstPass=i == 0)  # only time shift.
            # first pass is true for the first iteration and false for the second.
            emb.t0 += x0

            curves = ['spotsR', 'spotsY']  # use red and yellow spots
            x0, nDev, yScale, xScale = getAlignParams(emb, curves, show=False, shiftVary=False, scaleVary=False,
                                                      tScaleVary=True, verb=True, firstPass=i == 0)  # don't allow shift
            # in time, but allow time scale
            emb.tScale *= xScale
    # used for all embryos processing prior to 03-07-18
    # emb.t0 = np.round(emb.t0, 1)
    # emb.tScale = np.round(emb.tScale, 2)
    # NOTE: we did not rerun all embryos to reflect the change, because the change is marginal for results
    emb.t0 = float(np.round(emb.t0, getSig(emb.t0) + 5))  # changed from round to float to stabilize borderline values
    # in t0, i.e. 0.49 vs 0.5
    emb.tScale = float(np.round(emb.tScale, getSig(emb.tScale) + 5))
    return emb


if __name__ == '__main__':
    from Embryos import MSEmbryo
    folders = [FOLDER_IN + '20140319T140206/Emb{0}/'.format(i) for i in [13]]
    for folder in folders[:]:
        emb = MSEmbryo(folder)
        alignTimeMS(emb)
