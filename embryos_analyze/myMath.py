'''
Created on Jun 10, 2014

@author: renat

Collection of various mathematical functions
'''

import numpy as np
import scipy
import matplotlib.pyplot as plt
from numpy import pi
from scipy.optimize import curve_fit
from scipy.special import erf
import myFigure
from lmfit import Parameters, minimize, report_fit, Model, Parameter
import myFunc


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]


def fitLine(xdata, ydata, ysigma=None, xsigma=None):
    """
    Performs a linear fit to data.

    Parameters
    ----------
    xdata : An array of length N.
    ydata : An array of length N.
    sigma : None or an array of length N,
        If provided, it is the standard-deviation of ydata.
        This vector, if given, will be used as weights in the fit.

    Returns
    -------
    a, b   : Optimal parameter of linear fit (y = a*x + b)
    sa, sb : Uncertainties of the parameters
    """
    if xsigma is not None:
        w = 1.0 / (ysigma ** 2 + xsigma ** 2)
    elif ysigma is None:
        w = np.ones(len(ydata))  # Each point is equally weighted.
    else:
        w = 1.0 / (ysigma ** 2)

    w = w / sum(w)
    sw = sum(w)
    wx = w * xdata  # this product gets used to calculate swxy and swx2
    swx = sum(wx)
    swy = sum(w * ydata)
    swxy = sum(wx * ydata)
    swx2 = sum(wx * xdata)

    a = (sw * swxy - swx * swy) / (sw * swx2 - swx * swx)
    b = (swy * swx2 - swx * swxy) / (sw * swx2 - swx * swx)

    if ysigma is None:
        chi2 = sum(((a * xdata + b) - ydata) ** 2)
    else:
        chi2 = sum((((a * xdata + b) - ydata) / ysigma) ** 2)
    dof = len(ydata) - 2
    rchi2 = np.sqrt(chi2 / dof)

    sa = rchi2 * np.sqrt(sw / (sw * swx2 - swx * swx))
    sb = rchi2 * np.sqrt(swx2 / (sw * swx2 - swx * swx))
    return a, b, sa, sb


def fitSigm(x, y, yerr=None, base=True):
    def sigmoidal(x, a, b, m, s, r):
        return b - a * (1. + np.exp((x - m) / s)) ** (-r)

    show = False
    r0 = 1.  # exponent
    b0 = y[-1]  # np.max(y) #plateau at the end
    minInd = 0  # np.argmin(y[:y.size/2])
    a0 = b0 - np.min(y)  # b-a plateau at the beginning
    m0 = x[minInd:][np.argmin(np.abs(y[minInd:] - (b0 - a0 / 2.)))]  # x position of the middle intensity
    s0 = 3  # width
    p0 = (a0, b0, m0, s0, r0)
    params = Parameters()
    params.add('a', p0[0], min=0.5 * (max(y) - min(y)), max=1.2 * (max(y) - min(y)))
    if base:
        params.add('base', p0[1] - p0[0], min=0, max=0.3 * (max(y) - min(y)))
    else:
        params.add('base', 0, vary=False)
    params.add('b', p0[1], expr='base+a')
    params.add('m', p0[2], min=0, max=max(x))
    params.add('s', p0[3], min=0, max=max(x))
    params.add('r', p0[4], min=0, max=10, vary=False)
    mod = Model(sigmoidal)
    if yerr is not None:
        err = yerr[minInd:]
    else:
        err = np.ones_like(y)
    if show:
        fig = myFigure.myFigure()
        fig.errorbar(x[::10], y[::10], err[::10], join=False, color='k')
    try:
        res = mod.fit(y, x=x, params=params, weights=1. / err ** 2)
        popt = np.array([res.best_values['a'], res.best_values['b'], res.best_values['m'], res.best_values['s'],
                         res.best_values['r']])
        if show:
            fig.plot(x, sigmoidal(x, *popt), color='r')
            fig.plot(x, popt[0] / 4. / popt[3] * (x - popt[2]) + popt[1] - popt[0] / 2., 'k')
            fig.plot(x, sigmoidal(x, *p0), 'g')
            fig.ylim((0, popt[1] * 1.2))
    except Exception as e:
        print('could not fit: {0}'.format(str(e)))
        if show: fig.plot(x, sigmoidal(x, *p0), 'r')
        popt = np.ones(5) * np.nan
        perr = np.ones(5) * np.nan
    try:
        perr = np.array([res.params['a'].stderr, res.params['b'].stderr, res.params['m'].stderr, res.params['s'].stderr,
                         res.params['r'].stderr])
        perr[np.where(perr == 0)] = np.nan
    except:
        print('bad parameter errors')
        perr = np.nan * np.ones_like(popt)
    if show:
        print(report_fit(res))
        plt.show()
    return popt, perr


def fitSigmBump(x, y, yerr=None):
    def cost(x, a, m1, s1, b, m2, s2):
        return a / (1. + np.exp(-(x - m1) / s1)) + b * np.exp(-(x - m2) ** 2 / s2 ** 2)

    def getErr(y):
        dy = (np.roll(y, -1) - y)[:-1]
        ddy = (np.roll(dy, -1) - dy)[:-1]
        err = np.interp(np.arange(y.size), np.arange(ddy.size) + 1, np.abs(ddy))
        #     err= median_filter(err,3)
        errMin = 0.01
        err[np.isnan(err)] = 0
        err[np.where(err < errMin)] = errMin
        return err

    show = False
    minInd = 0
    maxV = 0
    peak = False
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]
    x = x[y > 0]
    y = y[y > 0]
    for endInd in range(y.size):
        val = y[endInd]
        if val > maxV and not peak: maxV = val
        if val > 1.1 * maxV and peak: break
        if val < 0.8 * maxV and val > 0.5:
            peak = True
    a0 = np.nanmedian(y[np.nonzero(y)][-5:])  # np.max(y) #plateau at the end
    b0 = (np.nanmax(y) - a0)
    dy = np.roll(y, -1) - y
    dy[np.isnan(dy)] = 0
    indMax = min(np.where(dy < -0.05)[0][0], endInd)
    if indMax < 5: indMax = 5
    m10 = x[minInd:][np.nanargmin(
        np.abs(y[minInd:indMax] - (a0 + np.nanmax(y)) / 2))] - 1  # x position of the middle intensity
    m20 = x[minInd:][np.nanargmin(
        np.abs(y[minInd:indMax] - (a0 + np.nanmax(y)) / 2))] + 3  # x position of the middle intensity
    s10 = 0.5  # width
    s20 = 2.  # width
    p0 = (a0, m10, s10, b0, m20, s20)
    params = Parameters()
    params.add('a', p0[0], min=0.9 * p0[0], max=1.5 * p0[0] + 1)
    params.add('m1', p0[1], min=p0[1] - 2, max=p0[1] + 2)
    params.add('s1', p0[2], min=0.5, max=max(x))
    params.add('b', p0[3], min=0.67 * p0[3], max=1.5 * p0[3] + 1)
    params.add('m2', p0[4], min=p0[4] - 4, max=p0[4] + 4)
    params.add('s2', p0[5], min=1, max=max(x))
    params.add('expr', value=1., expr='s2-2*s1', min=0.1)
    mod = Model(cost)
    if yerr is not None:
        err = yerr
    else:
        err = getErr(y)
    err = np.ones_like(err) * np.median(err)
    if show:
        fig = myFigure()
        fig.errorbar(x[minInd:endInd], y[minInd:endInd], err[minInd:endInd], color='k')
    try:
        res = mod.fit(y[minInd:endInd], x=x[minInd:endInd], params=params, weights=1. / err[minInd:endInd] ** 2)
        popt = np.array([res.best_values['a'], res.best_values['m1'], res.best_values['s1'], res.best_values['b'],
                         res.best_values['m2'], res.best_values['s2']])
        if show:
            fig.plot(x, cost(x, *popt), color='r')
            fig.plot(x, cost(x, *p0), 'g')
    except Exception as e:
        print('fitSigmBump: could not fit: {0}'.format(str(e)))
        if show: fig.plot(x, cost(x, *p0), 'r')
        popt = np.ones(6) * np.nan
        perr = np.ones(6) * np.nan
    try:
        perr = np.array(
            [res.params['a'].stderr, res.params['m1'].stderr, res.params['s1'].stderr, res.params['b'].stderr,
             res.params['m2'].stderr, res.params['s2'].stderr])
        perr[np.where(perr == 0)] = np.nan
    except:
        print('fitSigmBump: bad parameter errors')
        perr = np.nan * np.ones_like(popt)
    if show:
        try:
            print('fitSigmBump: successfully minimized = ', res.success)
            print('p0=', p0)
            print(report_fit(res))
        except:
            print('p0=', p0)
        plt.show()
    return popt, perr


def getSubIndex(a, left, right):
    if left > right: left, right = right, left
    if a[0] > a[-1]:
        return a.size - getSubIndex(a[::-1], left, right) - 1
    else:
        return np.where(a[np.where(a <= right)[0]] >= left)[0]


def getSig(x):
    x = abs(x)
    if x != np.inf and not np.isnan(x) and x != 0:
        return -int(np.floor(np.log10(x)))
    else:
        return 0


def calcMoment(points, center, axis):
    moment = 0.
    x0, y0, z0 = center
    for point in points:
        x, y, z = point
        if axis == 0:
            rSq = (y - y0) ** 2 + (z - z0) ** 2
        elif axis == 1:
            rSq = (x - x0) ** 2  # +(z-z0)**2 #use rotation independent projection
        elif axis == 2:
            rSq = (x - x0) ** 2  # +(y-y0)**2
        moment += rSq
    return moment


def getCenterOfMass(points):
    xc, yc, zc = 0., 0., 0.
    i = 0
    for point in points:
        if len(point) == 3:
            x, y, z = point
        else:
            x, y = point
            z = 0
        xc += x
        yc += y
        zc += z
        i += 1.
    xc /= i
    yc /= i
    zc /= i
    if len(point) == 3:
        return (xc, yc, zc)
    else:
        return (xc, yc)


class AvgCurve(object):
    def __init__(self):
        self.nMemb = 0

    def makeFirst(self, x, y):
        minX = np.min(x)
        maxX = np.max(x)
        self.dx = 0.1 * np.min((np.roll(x, -1) - x)[:-1])
        self.x = np.arange(np.floor(minX / self.dx), np.ceil(maxX / self.dx) + 1) * self.dx
        self.y = np.zeros_like(self.x)
        self.ysq = np.zeros_like(self.x)
        self.n = np.zeros_like(self.x)

    def add(self, x, y):
        if self.nMemb > 0:
            self.checkBounds(x, y)
        else:
            self.makeFirst(x, y)
        if x[0] <= self.x[-1]:
            indS = np.where(self.x >= x[0])[0][0]
        else:
            return
        indE = np.where(self.x >= x[-1])[0][0]
        yInt = np.interp(self.x[indS:indE], x, y)
        self.y[indS:indE] = (self.y[indS:indE] * self.n[indS:indE] + yInt) / (self.n[indS:indE] + 1.)
        self.ysq[indS:indE] = (self.ysq[indS:indE] * self.n[indS:indE] + yInt ** 2) / (self.n[indS:indE] + 1.)
        self.n[indS:indE] += 1
        self.nMemb += 1

    def checkBounds(self, x, y):
        minX, maxX = None, None
        if np.min(x) < self.x[0]: minX = np.min(x)
        if np.max(x) > self.x[-1]: maxX = np.max(x)
        if minX is not None or maxX is not None: self.extendX(minX, maxX)

    def extendX(self, minX, maxX):
        if minX is None: minX = self.x[0]
        if maxX is None: maxX = self.x[-1]
        xtmp = np.arange(np.floor(minX / self.dx), np.ceil(maxX / self.dx) + 1) * self.dx
        ytmp = np.zeros_like(xtmp)
        ysqtmp = np.zeros_like(xtmp)
        ntmp = np.zeros_like(xtmp)
        indS = int(np.floor(self.x[0] / self.dx) - np.floor(minX / self.dx))
        indE = indS + self.x.size
        xtmp[indS:indE] = self.x
        ytmp[indS:indE] = self.y
        ysqtmp[indS:indE] = self.ysq
        ntmp[indS:indE] = self.n
        self.x = xtmp
        self.y = ytmp
        self.ysq = ysqtmp
        self.n = ntmp

    def getAvg(self, nTh=1, stdErr=False, every=1):
        if nTh > self.nMemb: nTh = self.nMemb
        inds = np.where(self.n >= nTh)
        yerr = np.sqrt((self.ysq - self.y ** 2))
        if stdErr: yerr /= np.sqrt(self.n)
        return self.x[inds][::every], self.y[inds][::every], yerr[inds][::every]


if __name__ == '__main__':
    pass
