'''
Created on Aug 24, 2016

@author: renat, becky
'''

from myFigure import *
from varLookup import *
from emb_handler import loadEmbs, updateEmbsMulti, refreshEmbsMulti
from myFunc import readCSV, getStrVal
from myMath import nan_helper, fitLine
from lmfit import Parameters, minimize, report_fit, Parameter


def loadTestEmbryosControl():
    data = readCSV(FILENAME_ALIGNMENT_0000_MS, ',')[1:]  # skip first line
    data = [d[:9] + d[10:] for d in data]
    embs = []
    for d in data:
        folder = FOLDER_IN + 'cropped/EMBD0000/MS/{0}/{1}/'.format(d[2], d[4])
        embs.append(loadEmbs(folder)[0])
    #     for emb in embs:
    #         if not imageTester.checkEmbFocus(emb): print('check focus {d}/{l}'.format(l=emb.label, d=emb.date))
    #     embs = refreshEmbsMulti(embs)
    embs = updateEmbsMulti(embs)
    return embs, data


def loadTestEmbryosRNAi():
    data = readCSV(FILENAME_ALIGNMENT_RNAI_MS, ',')[1:]  # skip first line
    data = [d[:9] + d[10:] for d in data]
    d = data.pop(14)
    print('popped {0} {1}'.format(d[2], d[4]))
    embs = []
    dtmp = []
    for d in data:
        folder = FOLDER_IN + 'cropped/{0}/MS/{1}/'.format(d[3], d[4])
        if sum(np.isnan(d[6:12])) < 5 and not np.isnan(d[13]):
            #         if not np.isnan(d[13]):
            embs.append(loadEmbs(folder)[0])
            dtmp.append(d)

    # embs = refreshEmbsMulti(embs)
    embs = updateEmbsMulti(embs)
    return embs, dtmp


def reportMovementDetected(embs):
    tMove = []
    for emb in embs:
        if emb.movement: tMove.append((emb.params['tMove'] - emb.t0))
    print('movement statistics = {0}'.format(getStrVal(np.mean(tMove), np.std(tMove))))


def reportData(embs, data):
    timePoints = []
    t0s = []
    for emb, embInfo in zip(embs, data):
        timePoints.append((np.array(embInfo[6:12]) - emb.t0))
        t0s.append(emb.t0)
    timePoints = np.array(timePoints)
    print('avg time point spread = {0}'.format(np.std(np.nanmean(timePoints, axis=1))))


def findT0Obs(data):
    ''' t0 is the timepoint of average scoring as determined by manually scored benchmark timepoints '''
    return (np.nanmean(data[:, 6:12], axis=1))


def alignManualPoint(data, usePoints):
    nEmbs = data.shape[0]

    def cost(params):
        cost = np.array(
            [((data[i, j] - params['dt{0}'.format(i)]) * params['tScale{0}'.format(i)] - params['t{0}'.format(j)]) ** 2
             for j in usePoints for i in range(nEmbs)])
        nans, x = nan_helper(cost)
        return cost[~nans]

    p = Parameters()
    expr = '{n}'.format(n=nEmbs)
    for i in range(nEmbs):
        if i > 0:
            expr += '-tScale{0}'.format(i)
            p['tScale{0}'.format(i)] = Parameter(value=1., min=0.5, max=1.5, vary=True)
        if not np.isnan(data[i, 0]):
            p['dt{0}'.format(i)] = Parameter(value=data[i, 0], min=data[i, 0] - 5, max=data[i, 0] + 5)
        elif not np.isnan(np.nanmax(data[i])):
            p['dt{0}'.format(i)] = Parameter(value=0, min=-10, max=np.nanmax(data[i]))
    p['tScale0'] = Parameter(value=1., expr=expr, min=0.5, max=2.)
    tVals = [0., 5.2, 6.7, 8.2, 9.6, 12.7]
    for j in usePoints:
        p['t{0}'.format(j)] = Parameter(value=tVals[j], vary=(j != usePoints[0]), min=tVals[j] - 3, max=tVals[j] + 3)
    #         p['t{0}'.format(j)] = Parameter(value = tVals[j], vary = False , min = tVals[j]-3, max = tVals[j]+3)
    minzer = minimize(cost, params=p, method='cobyla')
    res = minzer.params
    print('successful fit ', minzer.success)
    print(report_fit(res))
    print('sum=', np.sum([res['tScale{0}'.format(i)].value for i in range(nEmbs)]))
    return np.array([res['dt{0}'.format(i)].value for i in range(nEmbs)]), \
           np.array([res['tScale{0}'.format(i)].value for i in range(nEmbs)]), \
           np.array([res['t{0}'.format(i)].value for i in usePoints])


def alignManualPointFixedLandmark(data, usePoints):
    nEmbs = data.shape[0]

    def cost(params):
        i = params['i'].value
        cost = np.array([((data[i, j] - params['dt{0}'.format(i)].value) * params['tScale{0}'.format(i)].value - params[
            't{0}'.format(j)].value) ** 2 for j in usePoints])
        nans, x = nan_helper(cost)
        if cost[~nans].size < 2:
            print('zero cost embryo {0}'.format(i))
            return [0., 0.]
        return cost[~nans]

    t0, tScale = [], []
    for i in range(nEmbs):
        p = Parameters()
        p['tScale{0}'.format(i)] = Parameter(value=1., min=0.5, max=2., vary=True)
        p['i'] = Parameter(value=i, vary=False)
        if not np.isnan(data[i, 0]):
            p['dt{0}'.format(i)] = Parameter(value=data[i, 0], min=data[i, 0] - 5, max=data[i, 0] + 5)
        elif not np.isnan(np.nanmax(data[i])):
            p['dt{0}'.format(i)] = Parameter(value=0, min=-10, max=np.nanmax(data[i]))
        tVals = [0., 5.2, 6.7, 8.2, 9.6, 12.7]
        for j in usePoints:
            p['t{0}'.format(j)] = Parameter(value=tVals[j], vary=False)
        minzer = minimize(cost, params=p)
        res = minzer.params
        t0.append(res['dt{0}'.format(i)].value)
        tScale.append(res['tScale{0}'.format(i)].value)
    return np.array(t0), np.array(tScale), np.array(tVals)


def plotDist(data, usePoints, t0obs, tScaleObs):
    fig = myFigure()
    for i in usePoints:
        d = (data[:, i] - t0obs) * tScaleObs
        nans, x = nan_helper(d)
        print('var(t{0})={1}'.format(i, np.std(d[~nans])))
        fig.hist(d[~nans], alpha=0.5)
    fig.xlim([-5, 20])


def plotTObsTa(embs, t0obs, tScalesObs, fig=None, figS=None, dtcList=None):
    if fig is None: fig = myFigure()
    if figS is None: figS = myFigure()
    figMoI0R = myFigure()
    figCoM0R = myFigure()
    if dtcList is not None:
        dtc = np.mean(dtcList)
    else:
        dtc = None
    t0a, tsa = [], []
    i = 0
    for emb in embs:
        figMoI0R = emb.show('MoI0R', fig=figMoI0R, setColor=False)
        figCoM0R = emb.show('CoM0R', fig=figCoM0R, setColor=False)
        t0a.append(emb.t0)
        tsa.append(emb.tScale)
        #         tsa.append(1)
        print('{3} {2} {0} t0={1}, tScale={ts}'.format(emb.label, emb.t0, emb.RNAi, emb.date, ts=emb.tScale), t0obs[i],
              tScalesObs[i], abs(tScalesObs[i] - emb.tScale))
        i += 1
    figMoI0R.close()
    figCoM0R.close()
    t0a = np.array(t0a)
    tsa = np.array(tsa)
    if dtc is None:
        fig.scatter(t0a, t0obs, color='g')
        figS.scatter(tsa, tScalesObs, color='g')
        a, b, sa, sb = fitLine(tsa, tScalesObs)
        print('lineFit={0}x+{1}'.format(getStrVal(a, sa), getStrVal(b, sb)))
        #         figSp = myFigure()
        #         figSp.hist(tScalesObs, bins=10, normed=True)
        #         x = np.linspace(np.min(tScalesObs), np.max(tScalesObs), 100)
        #         y = norm.pdf(x, np.mean(tScalesObs), np.std(tScalesObs))
        #         figSp.plot(x, y)
        #         figSp.show()
        qt = np.std(t0a - t0obs)
        dtScale = tScalesObs - tsa
        qts = np.sqrt(np.sum(dtScale ** 2) / dtScale.size)
        print('quality=', qt, qts)
        figS.plot(tsa, a * tsa + b, 'r--')
        figS.title('l={}, quality = {}'.format(nLambda_MS, qts))
        fig.title('l={}, quality = {}'.format(nLambda_MS, qt))
    else:
        fig.scatter(t0a, t0obs, color='b')
    dtList = t0obs - t0a
    dtScale = tScalesObs - tsa
    dt = np.mean(dtList)
    #     print('nLambda_GLS={1}; alignment quality = {0}'.format(q, nLambda_GLS))
    fig.plot(t0a, t0a + dt, color='k')
    if dtc is not None:
        print('dt difference is {0}'.format(abs(dt - dtc)))
        fig.plot(t0a, t0a + dtc, 'r--')
        #         figS.scatter(tsa, tScalesObs, color='b')
        fig.title('MS |dt-dtcont|={0}'.format(abs(dt - dtc)))
    fig.xlabel('t0 automated')
    fig.ylabel('t0 manual')
    figS.xlabel('tScale automated')
    figS.ylabel('tScale manual')
    figS.plot((0.5, 1.5), (0.5, 1.5), color='k')
    figS.ylim((0.5, 1.5))
    figS.xlim((0.5, 1.5))
    #     fig.ylim((-2,10))
    fig.xlim((-2, 14))
    if True or dtc is None:
        fig.save(
            FOLDER_IN + 'Automated_analysis/timeAlignment/MS/060917/t0_lam={l}_tau={t}_q={q}.svg'.format(l=nLambda_MS,
                                                                                                         q=qt,
                                                                                                         t=tauScale_MS))
        figS.save(FOLDER_IN + 'Automated_analysis/timeAlignment/MS/060917/tScale_lam={l}_tau={t}_q={q}.svg'.format(
            l=nLambda_MS, q=qts, t=tauScale_MS))
    return fig, figS, dtList, dtScale


def plotTDev(tdm, embs):
    fig = myFigure()
    tda = []
    for emb in embs:
        sc = [emb.cutoff[c] for c in ['spotsG', 'spotsR', 'spotsY']]
        tda.append(min(sc))
        #         tda.append(min(emb.cutoff.values()))
        print(
        '{l} {r} {d} tda={ta}, tdm={tm}, delta={dlt}'.format(l=emb.label, r=emb.RNAi, d=emb.date, tm=tdm[len(tda) - 1],
                                                             ta=tda[-1], dlt=abs(tdm[len(tda) - 1] - tda[-1])))
    fig.scatter(tda, tdm)
    fig.plot((0, 31), (0, 31), 'k')
    fig.ylabel('manual time of dev')
    fig.xlabel('automatic time of dev')
    fig.title('deviation validation for lambda={l}'.format(l=nLambda_GLS))
    return fig


def plotPredictedCurves(embs, t0obs, tScaleObs):
    fig = {}
    for curveName in AVG_CURVES_MS:
        fig[curveName] = myFigure()
        fig[curveName].title(curveName)
    for i in range(len(embs)):
        embs[i].t0 = t0obs[i]
        embs[i].tScale = tScaleObs[i]
        embs[i].setScaleCutoff()
        print('t0={t}, tScale={ts}'.format(t=embs[i].t0, ts=embs[i].tScale))
        for curveName in AVG_CURVES_MS:
            #             embs[i].scale[curveName] = 1.
            fig[curveName] = embs[i].show(curveName, fig=fig[curveName])
    fig[curveName].show()


def initialize():
    # embs, data = loadTestEmbryosControl()
    embs, data = loadTestEmbryosRNAi()
    data = np.array([d[6:12] for d in data])
    usePoints = range(data.shape[1])
    t0obs, tScaleObs, tdev = alignManualPointFixedLandmark(data, usePoints)
    #     plotPredictedCurves(embs, t0obs, tScaleObs)
    fig, figS, dtcList, dtScale = plotTObsTa(embs, t0obs, tScaleObs)
    #     plotDist(data, usePoints, t0obs, tScaleObs)

    #     embs, data = loadTestEmbryosRNAi()
    #     tdm = np.array([d[12] for d in data])
    #     data = np.array([d[6:12] for d in data])
    #     plotTDev(tdm, embs)
    #     usePoints = range(data.shape[1])
    #     t0obs, tScaleObs, tdev = alignManualPointFixedLandmark(data, usePoints)
    #     fig, figS, dtList, rtScale = plotTObsTa(embs, t0obs,tScaleObs, fig, figS, dtcList)
    #     dtc = np.mean(dtcList)
    #     qt = np.sqrt(np.sum((np.concatenate((dtcList,dtList))-dtc)**2)/(dtcList.size+dtList.size))
    #     qts = np.sqrt(np.sum(dtScale**2)/dtScale.size)
    #     fig.save(FOLDER_IN+'Automated_analysis/timeAlignment/MS/040417/t0_lam={l}_tau={t}_q={q}.svg'.format(l=nLambda_GLS, q=qt, t=tauScale_GLS))
    #     figS.save(FOLDER_IN+'Automated_analysis/timeAlignment/MS/040417/tScale_lam={l}_tau={t}_q={q}.svg'.format(l=nLambda_GLS, q=qts, t=tauScale_GLS))
    plt.show()


if __name__ == '__main__':
    initialize()
