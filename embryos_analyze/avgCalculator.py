"""
Created on Aug 19, 2016

@author: renat
"""
from varLookup import FOLDER_IN, printLog, AVG_CURVES_MS, AVG_CURVES_GLS, FILENAME_AVG_CURVES_MS_ADJTH, \
    FILENAME_AVG_CURVES_GLS
import emb_handler
from myFigure import myFigure
import os
import numpy as np
import pickle
from myFunc import getStrVal
from myMath import AvgCurve, fitSigm, getSubIndex
from lmfit import Parameters, Model
from lmfit.printfuncs import report_fit

global nEmb
# nEmb = 10
nEmb = 30


def getMSEmbs():
    embsTmp, embs = [], []
    folder = FOLDER_IN
    use_filepath = get_control_folders_ms()
    folders = [folder + d for d in use_filepath][:100]
    # folder = FOLDER_IN + 'cropped/EMBD0000/MS/'
    # folders = [folder + d + '/' for d in os.listdir(folder)][:20]

    i = 0
    for folder in folders:
        print('{0} folders left'.format(len(folders) - i))
        if folder not in badFolders:
            embsTmp += emb_handler.loadEmbs(folder)
        i += 1
    # embsTmp = emb_handler.updateEmbsMulti(embsTmp)
    for emb in embsTmp:
        if emb.movement == np.nan or emb.movement == None:
            print emb.label, emb.date
        if emb.movement:
            embs.append(emb)
        printLog(
            'getMSEmbs {date} {0}: t0={t0}, tScale={ts}, {1}'.format(emb.label, emb.scale, t0=emb.t0, date=emb.date,
                                                                     ts=emb.tScale))
    return embs


def getGLSEmbs():
    embsTmp, embs = [], []
    folder = FOLDER_IN
    use_filepath = get_control_folders_gls()
    folders = [folder + d for d in use_filepath][:100]
    # folder = FOLDER_IN + 'cropped/EMBD0000/GLS/'
    # folders = [folder + d + '/' for d in os.listdir(folder)][:20]
    badFolders = [FOLDER_IN + 'cropped/EMBD0000/GLS/20150129T130903/',
                  FOLDER_IN + 'cropped/EMBD0000/GLS/20140813T133845/']
    for folder in folders:
        if folder not in badFolders:
            embsTmp += emb_handler.loadEmbs(folder)
    # embsTmp = emb_handler.updateEmbsMulti(embsTmp)
    for emb in embsTmp:
        if emb.movement and emb.checkGoodParam('mG'):
            embs.append(emb)
        printLog('getGLSEmbs {date} {0}: t0={t0} {1}'.format(emb.label, emb.scale, t0=emb.t0, date=emb.date))
    return embs


def calcAvgCurvesMS(embs, show=True):
    curves = {}
    fig = {}
    for curveName in AVG_CURVES_MS:
        curves[curveName] = AvgCurve()
        if show:
            fig[curveName] = myFigure()
            fig[curveName].title(curveName)
    printLog('TOTAL EMBRYOS:{0}'.format(len(embs)))
    scales = {}
    nPoints = {}
    tScales = []
    for cn in AVG_CURVES_MS:
        scales[cn] = []
        nPoints[cn] = []
    for emb in embs:
        for curveName in AVG_CURVES_MS:
            x, y = emb.getCurve(curveName)
            curves[curveName].add(x[~np.isnan(y)], y[~np.isnan(y)])
            scales[curveName].append(emb.scale[curveName])
            tScales.append(emb.tScale)
            nPoints[curveName].append(emb.cutoff[curveName] - emb.startPoint[curveName])
            if show:
                fig[curveName] = emb.show(curveName, fig[curveName], setColor=False)
    if show:
        for curveName in AVG_CURVES_MS:
            x, y, yerr = curves[curveName].getAvg(nEmb)
            fig[curveName].errorbar(x[::10], y[::10], yerr[::10], color='k')
    # normalize the curves to have average scale 1 and dt 0
    popt = fitSigm(*curves['tIntG'].getAvg(nEmb))[0]
    dt = popt[2] - 2. * popt[3]
    printLog('calcAvgCurvesMS: dt = {0}'.format(dt))
    for curveName in AVG_CURVES_MS:
        w = np.array(nPoints[curveName]).astype(np.float)
        w /= sum(w)
        scale = np.sum(scales[curveName] * w)
        printLog('calcAvgCurvesMS: {cn} scale = {s}'.format(s=scale, cn=curveName))
        curves[curveName].x -= dt
        curves[curveName].x /= np.mean(tScales)
        curves[curveName].y /= scale
        curves[curveName].ysq /= scale ** 2
    if show:
        fig[curveName].show()
    return curves


def calcAvgCurvesGLS(embs, show=False):
    curves = {}
    fig = {}
    for curveName in AVG_CURVES_GLS:
        curves[curveName] = AvgCurve()
        if show:
            fig[curveName] = myFigure()
            fig[curveName].title(curveName)
    printLog('TOTAL EMBRYOS:{0}'.format(len(embs)))
    scales = {}
    nPoints = {}
    tScales = []
    for cn in AVG_CURVES_GLS:
        scales[cn] = []
        nPoints[cn] = []
    for emb in embs:
        for curveName in AVG_CURVES_GLS:
            x, y = emb.getCurve(curveName)
            curves[curveName].add(x[~np.isnan(y)], y[~np.isnan(y)])
            scales[curveName].append(emb.scale[curveName])
            nPoints[curveName].append(emb.cutoff[curveName] - emb.startPoint[curveName])
            tScales.append(emb.tScale)
            if show: fig[curveName] = emb.show(curveName, fig[curveName], setColor=False)
    if show:
        for curveName in AVG_CURVES_GLS:
            x, y, yerr = curves[curveName].getAvg(nEmb)
            fig[curveName].errorbar(x[::10], y[::10], yerr[::10], color='k')
    # normalize the curves to have average scale 1 and dt 0
    popt = fitSigmSpots(*curves['spotsR'].getAvg(nEmb))[0]
    dt = popt[2] - 2. * popt[3]
    printLog('calcAvgCurvesGLS: dt = {0}'.format(dt))
    for curveName in AVG_CURVES_GLS:
        w = np.array(nPoints[curveName]).astype(np.float)
        w /= sum(w)
        scale = np.sum(scales[curveName] * w)
        printLog('calcAvgCurvesGLS: {cn} scale = {s}, tScales={ts}'.format(s=scale, cn=curveName, ts=np.mean(tScales)))
        curves[curveName].x -= dt
        curves[curveName].x /= np.mean(tScales)
        curves[curveName].y /= scale
        curves[curveName].ysq /= scale ** 2
        if show: fig[curveName].legend(2)

    if show: fig[curveName].show()
    return curves


def saveAvgCurvesMS(curves):
    fileName = FILENAME_AVG_CURVES_MS_ADJTH
    curvesSave = {}
    for curveName in AVG_CURVES_MS:
        curvesSave[curveName] = curves[curveName].getAvg(nEmb)
    output = open(fileName, 'wb')
    pickle.dump(curvesSave, output, pickle.HIGHEST_PROTOCOL)


def saveAvgCurvesGLS(curves):
    fileName = FILENAME_AVG_CURVES_GLS
    curvesSave = {}
    for curveName in AVG_CURVES_GLS:
        curvesSave[curveName] = curves[curveName].getAvg(nEmb)
    output = open(fileName, 'wb')
    pickle.dump(curvesSave, output, pickle.HIGHEST_PROTOCOL)


def calcAvgParamsMS(embsMS):
    aR, aG, bR, bG = [], [], [], []
    sR, sG = [], []
    t0Gmt0R = []
    mg2move = []
    t02move = []
    figMS = myFigure()
    for emb in embsMS:
        if emb.checkGoodParam('aR') and emb.checkGoodParam('aG') \
                and emb.checkGoodParam('mR') and emb.checkGoodParam('mG') \
                and emb.checkGoodParam('sR') and emb.checkGoodParam('sG'):
            figMS = emb.showSigmFit(fig=figMS)
            aR.append(emb.params['aR'])
            aG.append(emb.params['aG'])
            bR.append(emb.params['bR'])
            bG.append(emb.params['bG'])
            sR.append(emb.params['sR'])
            sG.append(emb.params['sG'])
            t02move.append(emb.tMove - emb.t0)
            t0Gmt0R.append((emb.params['mG'] - 2. * emb.params['sG']) - (emb.params['mR'] - 2. * emb.params['sR']))
        if not np.isnan(emb.params['tMove']): mg2move.append(emb.params['tMove'] - emb.params['mG'])
    printLog('AR_AVG_MS={0}\nAG_AVG_MS={1}'.format(getStrVal(np.mean(aR), np.std(aR)), \
                                                   getStrVal(np.mean(aG), np.std(aG))))
    printLog('SR_AVG_MS={0}\nSG_AVG_MS={1}'.format(getStrVal(np.mean(sR), np.std(sR)), \
                                                   getStrVal(np.mean(sG), np.std(sG))))
    printLog('AVG_R_BASE_MS={0}\nAVG_G_BASE_MS={1}'.format(getStrVal(np.mean(bR) - np.mean(aR), np.std(aR)), \
                                                           getStrVal(np.mean(bG) - np.mean(aG), np.std(aG))))
    printLog('T0G_M_T0R_AVG_MS = {0}'.format(getStrVal(np.mean(t0Gmt0R), np.std(t0Gmt0R))))
    printLog('T0_2_MOVE_AVG_MS = {0}'.format(getStrVal(np.mean(t02move), np.std(t02move))))


def calcLengthRad(embs):
    length, radius = [], []
    for emb in embs:
        emb.refresh()
        l, r = emb.getSizeFromIm()
        length.append(l)
        radius.append(r)
    printLog('radius={0}'.format(getStrVal(np.mean(radius), np.std(radius))))
    printLog('length={0}'.format(getStrVal(np.mean(length), np.std(length))))


def calcAvgParamsGLS(embs):
    aG, aR, aY = [], [], []
    sG, sR, sY = [], [], []
    t0rmr0g, t0ymt0g, t0ymt0r = [], [], []
    for emb in embs:
        if emb.checkGoodParam('aG') and emb.checkGoodParam('aR') and emb.checkGoodParam('aY') \
                and emb.checkGoodParam('mG') and emb.checkGoodParam('mR') and emb.checkGoodParam('mY') \
                and emb.checkGoodParam('sR') and emb.checkGoodParam('sG') and emb.checkGoodParam('sY'):
            aR.append(emb.params['aR'])
            aG.append(emb.params['aG'])
            aY.append(emb.params['aY'])
            sR.append(emb.params['sR'])
            sG.append(emb.params['sG'])
            sY.append(emb.params['sY'])
            t0rmr0g.append(emb.params['mR'] - 2. * emb.params['sR'] - emb.params['mG'] + 2. * emb.params['sG'])
            t0ymt0g.append(emb.params['mY'] - 2. * emb.params['sY'] - emb.params['mG'] + 2. * emb.params['sG'])
            t0ymt0r.append(emb.params['mY'] - 2. * emb.params['sY'] - emb.params['mR'] + 2. * emb.params['sR'])
    printLog('AG_AVG_GLS = {0}'.format(getStrVal(np.mean(aG), np.std(aG))))
    printLog('AR_AVG_GLS = {0}'.format(getStrVal(np.mean(aR), np.std(aR))))
    printLog('AY_AVG_GLS = {0}'.format(getStrVal(np.mean(aY), np.std(aY))))
    printLog('SG_AVG_GLS = {0}'.format(getStrVal(np.mean(sG), np.std(sG))))
    printLog('SR_AVG_GLS = {0}'.format(getStrVal(np.mean(sR), np.std(sR))))
    printLog('SY_AVG_GLS = {0}'.format(getStrVal(np.mean(sY), np.std(sY))))
    printLog('T0R_M_T0G_AVG_GLS = {0}'.format(getStrVal(np.mean(t0rmr0g), np.std(t0rmr0g))))
    printLog('T0Y_M_T0G_AVG_GLS = {0}'.format(getStrVal(np.mean(t0ymt0g), np.std(t0ymt0g))))
    printLog('T0Y_M_T0R_AVG_GLS = {0}'.format(getStrVal(np.mean(t0ymt0r), np.std(t0ymt0r))))


def iterAvgCurvesMS(embs):
    maxIter = 10
    i = 0
    diff = 1
    curves0 = None
    while diff > 0.01 and i < maxIter:
        curves = calcAvgCurvesMS(embs, show=True)
        # saveAvgCurvesMS(curves) #FIXME
        for emb in embs:
            emb.loadAvgCurves()
        printLog('average curves calculated')
        embs = emb_handler.timeAlignMultiMS(embs)
        printLog('embryos time alignment complete')
        embs = emb_handler.setScalesEmbsMulti(embs)
        printLog('embryo scaling set')
        if curves0 is not None:
            diff = calcDiff(curves0, curves)
            printLog('difference calculated')
        curves0 = curves
        printLog('iterAvgCurvesMS diff={0}'.format(diff))
        i += 1


def iterAvgCurvesGLS(embs):
    maxIter = 10
    i = 0
    diff = 1
    curves0 = None
    while diff > 0.01 and i < maxIter:
        curves = calcAvgCurvesGLS(embs, show=True)
        # saveAvgCurvesGLS(curves) #FIXME
        for emb in embs:
            emb.loadAvgCurves()
        embs = emb_handler.timeAlignMultiGLS(embs)
        embs = emb_handler.setScalesEmbsMulti(embs)
        if curves0 is not None: diff = calcDiff(curves0, curves)
        curves0 = curves
        printLog('iterAvgCurvesGLS diff={0}'.format(diff))
        i += 1


def calcDiff(curves0, curves1):
    diff = np.array([])
    for key, val0 in curves0.iteritems():
        val1 = curves1[key]
        x0, y0, yerr = val0.getAvg(nEmb)
        x1, y1, yerr = val1.getAvg(nEmb)
        xmin = max(x0[0], x1[0])
        xmax = min(x0[-1], x1[-1])
        ind0 = getSubIndex(x0, xmin, xmax)
        y0 = y0[ind0]
        y1 = np.interp(x0[ind0], x1, y1)
        diff = np.concatenate((diff, (y0 - y1) ** 2 / np.mean(y1 ** 2)))
    return np.nanmean(diff)


def checkScalesMS(embs):
    from scipy.stats.stats import pearsonr

    scalesAll = {}
    radiusAll = []
    lengthAll = []
    fig = None
    for cn in AVG_CURVES_MS:
        scalesAll[cn] = []
    for emb in embs:
        for curveName in AVG_CURVES_MS:
            scalesAll[curveName].append(emb.scale[curveName])
        x, y = emb.getSizeFromIm()
        radiusAll.append(x)
        lengthAll.append(y)
    for cn in AVG_CURVES_MS:
        printLog('{cn} scale is {sc}'.format(cn=cn, sc=getStrVal(np.mean(scalesAll[cn]), np.std(scalesAll[cn]))))
    for i in range(len(AVG_CURVES_MS) - 1):
        curveName1 = AVG_CURVES_MS[i]
        for j in range(len(AVG_CURVES_MS[i + 1:])):
            curveName2 = AVG_CURVES_MS[i + 1 + j]
            pcc = pearsonr(scalesAll[curveName1], scalesAll[curveName2])[0]
            printLog('{0} and {1} pcc = {2}'.format(curveName1, curveName2, pcc))
            if abs(pcc) > 0.7:
                fig = myFigure()
                fig.scatter(scalesAll[curveName1], scalesAll[curveName2])
                fig.xlabel(curveName1)
                fig.ylabel(curveName2)
                fig.title('pcc = {0}'.format(pcc))
    for curveName in AVG_CURVES_MS:
        pcc = pearsonr(scalesAll[curveName], radiusAll)[0]
        printLog('{0} and {1} pcc = {2}'.format(curveName, 'radius', pcc))
        if abs(pcc) > 0.7:
            fig = myFigure()
            fig.scatter(scalesAll[curveName], radiusAll)
            fig.xlabel(curveName)
            fig.ylabel('radius')
            fig.title('pcc = {0}'.format(pcc))
        pcc = pearsonr(scalesAll[curveName], lengthAll)[0]
        printLog('{0} and {1} pcc = {2}'.format(curveName, 'length', pcc))
        if abs(pcc) > 0.7:
            fig = myFigure()
            fig.scatter(scalesAll[curveName], lengthAll)
            fig.xlabel(curveName)
            fig.ylabel('length')
            fig.title('pcc = {0}'.format(pcc))
    if fig is not None: fig.show()


def checkScalesGLS(embs):
    from scipy.stats.stats import pearsonr

    scalesAll = {}
    radiusAll = []
    lengthAll = []
    fig = None
    for cn in AVG_CURVES_GLS:
        scalesAll[cn] = []
    for emb in embs:
        for curveName in AVG_CURVES_GLS:
            scalesAll[curveName].append(emb.scale[curveName])
        x, y = emb.getSizeFromIm()
        radiusAll.append(x)
        lengthAll.append(y)
    for cn in AVG_CURVES_GLS:
        printLog('{cn} scale is {sc}'.format(cn=cn, sc=getStrVal(np.mean(scalesAll[cn]), np.std(scalesAll[cn]))))
    for i in range(len(AVG_CURVES_GLS) - 1):
        curveName1 = AVG_CURVES_GLS[i]
        for j in range(len(AVG_CURVES_GLS[i + 1:])):
            curveName2 = AVG_CURVES_GLS[i + 1 + j]
            pcc = pearsonr(scalesAll[curveName1], scalesAll[curveName2])[0]
            printLog('{0} and {1} pcc = {2}'.format(curveName1, curveName2, pcc))
            if abs(pcc) > 0.7:
                fig = myFigure()
                fig.scatter(scalesAll[curveName1], scalesAll[curveName2])
                fig.xlabel(curveName1)
                fig.ylabel(curveName2)
                fig.title('pcc = {0}'.format(pcc))
    for curveName in AVG_CURVES_GLS:
        pcc = pearsonr(scalesAll[curveName], radiusAll)[0]
        printLog('{0} and {1} pcc = {2}'.format(curveName, 'radius', pcc))
        if abs(pcc) > 0.7:
            fig = myFigure()
            fig.scatter(scalesAll[curveName], radiusAll)
            fig.xlabel(curveName)
            fig.ylabel('radius')
            fig.title('pcc = {0}'.format(pcc))
        pcc = pearsonr(scalesAll[curveName], lengthAll)[0]
        printLog('{0} and {1} pcc = {2}'.format(curveName, 'length', pcc))
        if abs(pcc) > 0.7:
            fig = myFigure()
            fig.scatter(scalesAll[curveName], lengthAll)
            fig.xlabel(curveName)
            fig.ylabel('length')
            fig.title('pcc = {0}'.format(pcc))
    if fig is not None: fig.show()


def fitSigmSpots(x, y, yerr):
    def sigmoidal(x, a, m, s, r):
        return a * (1. - (1. + np.exp((x - m) / s)) ** (-r))

    yerr = np.ones_like(y)
    show = True
    endInd = np.argmin((x - 9.) ** 2)  # Approximate frame of movement
    r0 = 1.  # exponent
    startInd = 0  # max(0,np.where(y>0)[0][0]-1)
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
    res = mod.fit(y[:endInd], x=x[:endInd], params=params, weights=1. / yerr[:endInd] ** 2)
    popt = np.array([res.best_values['a'], res.best_values['m'], res.best_values['s'], res.best_values['r']])
    perr = []
    if show:
        fig = myFigure()
        #             fig.errorbar(x[:endInd:10], y[:endInd:10], yerr[:endInd:10])
        fig.plot(x[:endInd:10], y[:endInd:10], color='r')
        fig.plot(x, sigmoidal(x, *popt), color='k')
        fig.plot(x, sigmoidal(x, *p0), color='g')
        print(report_fit(res))
        fig.show()
    return popt, perr


def get_control_folders_gls():
    dates = ['20140314T153735', '20140319T140206', '20140320T134454', '20140321T162110', '20140326T125744',
             '20140327T135219', '20140403T124950', '20140404T133853', '20140409T134416', '20140410T123736',
             '20140416T140401', '20140417T144336', '20140418T143933', '20140423T125529', '20140424T140159',
             '20140430T140422', '20140501T135409', '20140507T121259', '20140508T141139', '20140509T143953',
             '20140514T130911', '20140515T131046', '20140516T130559', '20140521T132613', '20140528T141502',
             '20140529T134135', '20140627T145639', '20140701T154139', '20140709T132509', '20140710T140504',
             '20140718T151611', '20140723T132003', '20140724T141838', '20140806T135105', '20140807T135038',
             '20140808T153637', '20140814T135456', '20140815T143554', '20140820T141612', '20140828T142637',
             '20140829T134745', '20140903T133504', '20140904T131050', '20140905T140506', '20140910T122848',
             '20140912T141649', '20140917T115255', '20140926T144329', '20141001T133638', '20141002T132022',
             '20141016T140752', '20141023T130422', '20141024T143550', '20141105T123638', '20141107T143327',
             '20141112T125319', '20141113T142111', '20141114T125207', '20141119T125538', '20141120T140815',
             '20141121T142152', '20141125T144501', '20141204T131945', '20141217T120303', '20141218T131231',
             '20141219T134142', '20150107T132211', '20150108T133738', '20150115T135137', '20150116T134901',
             '20150120T135640', '20150121T125804', '20150122T131308', '20150128T120821', '20150130T135638',
             '20150204T120204', '20150205T124453', '20150206T140216', '20150211T121018', '20150212T122807',
             '20150213T134117', '20150218T125636', '20150225T125230', '20150226T135321', '20150304T125310',
             '20150305T131306', '20150306T153126', '20150311T135645', '20150312T121507', '20150318T130640',
             '20150319T133546', '20150326T130210', '20150402T143618', '20150403T142853', '20150408T123836',
             '20150409T124850', '20150415T121757', '20150416T145930', '20150422T125031', '20150423T125304',
             '20150506T124235', '20150507T131111', '20150513T124944', '20150514T125824', '20150520T115933',
             '20150521T135040', '20150527T123035', '20150528T131406', '20150603T133456', '20150604T131219',
             '20150610T130133', '20150611T132858', '20150617T124551', '20150618T131713', '20150701T125127',
             '20150702T125423', '20150708T121623', '20150709T123422', '20150715T135249', '20150716T124801',
             '20150723T130801', '20150729T121201', '20150730T122827', '20150805T123025', '20150806T130215',
             '20150812T131345', '20150813T123805', '20150819T131249', '20150820T130314', '20150826T123412',
             '20150827T130101', '20150902T134023', '20150903T122831', '20150910T133240', '20150930T122117',
             '20151007T122635', '20151008T120416', '20151014T125902', '20151015T125611', '20151029T125426',
             '20151104T130721', '20151105T124752', '20151112T131012', '20151119T132259', '20151202T123250',
             '20151203T123705', '20160114T150421', '20160120T123135', '20160121T130634', '20160127T134318',
             '20160204T163336', '20160210T120517', '20160211T120752', '20160217T134018', '20160218T143452',
             '20160225T112241', '20160302T123313', '20160309T134346', '20160310T134056', '20160316T123005',
             '20160317T132308']
    use_filepaths = []
    for d in dates[:]:
        filepath = 'cropped/EMBD0000/GLS/{0}/'.format(d)
        use_filepaths.append(filepath)
    return use_filepaths


def get_control_folders_ms():
    dates = ['20140425T140657', '20140430T140422', '20140501T135409', '20140507T121259', '20140508T141139',
             '20140509T143953', '20140516T130559', '20140521T132613', '20140528T141502', '20140529T134135',
             '20140530T135628', '20140701T154139', '20140710T140504', '20140718T151611',
             '20140723T132003', '20140724T141838', '20140806T135105', '20140807T135038', '20140808T153637',
             '20140813T133845', '20140814T135456', '20140815T143554', '20140820T141612', '20140821T143207',
             '20140828T142637', '20140829T134745', '20140904T131050', '20140905T140506', '20140910T122848',
             '20140912T141649', '20140917T115255', '20140925T122500', '20140926T144329', '20141001T133638',
             '20141002T132022', '20141016T140752', '20141023T130422', '20141024T143550', '20141105T123638',
             '20141107T143327', '20141112T125319', '20141113T142111', '20141114T125207', '20141120T140815',
             '20141121T142152', '20141125T144501', '20141204T131945', '20141217T120303', '20141218T131231',
             '20150107T132211', '20150108T133738', '20150115T135137', '20150116T134901', '20150120T135640',
             '20150121T125804', '20150122T131308', '20150128T120821', '20150130T135638',
             '20150204T120204', '20150205T124453', '20150206T140216', '20150211T121018', '20150212T122807',
             '20150213T134117', '20150218T125636', '20150226T135321', '20150304T125310', '20150305T131306',
             '20150306T153126', '20150311T135645', '20150312T121507', '20150318T130640', '20150319T133546',
             '20150325T123625', '20150326T130210', '20150402T143618', '20150403T142853', '20150408T123836',
             '20150409T124850', '20150415T121757', '20150416T145930', '20150422T125031', '20150423T125304',
             '20150513T124944', '20150514T125824', '20150520T115933', '20150521T135040',
             '20150527T123035', '20150528T131406', '20150601T124426', '20150604T131219', '20150610T130133',
             '20150611T132858', '20150617T124551', '20150618T131713', '20150701T125127', '20150702T125423',
             '20150708T121623', '20150709T123422', '20150715T135249', '20150716T124801', '20150723T130801',
             '20150729T121201', '20150730T122827', '20150805T123025', '20150806T130215', '20150812T131345',
             '20150813T123805', '20150819T131249', '20150820T130314', '20150826T123412', '20150827T130101',
             '20150902T134023', '20150903T122831', '20150910T133240', '20150930T122117', '20151007T122635',
             '20151008T120416', '20151028T124731', '20151029T125426', '20151104T130721', '20151112T131012',
             '20151119T132259', '20151202T123250', '20151203T123705', '20160114T150421', '20160120T123135',
             '20160121T130634', '20160127T134318', '20160204T163336', '20160210T120517', '20160211T120752',
             '20160217T134018', '20160218T143452', '20160225T112241', '20160302T123313', '20160309T134346',
             '20160310T134056', '20160316T123005', '20160317T132308', '20160323T134818', '20160330T132134',
             '20160413T132939', '20160414T124411', '20160427T140158', '20160428T130512', '20160505T125834',
             '20160511T134745', '20160512T125616', '20160519T140603', '20160525T141941', '20160526T132118',
             '20160601T133958', '20160602T131501', '20160609T133750', '20160616T145853', '20160623T131147',
             '20160629T124228', '20160630T120211', '20160707T124255']
    use_filepaths = []
    for d in dates[:]:
        filepath = 'cropped/EMBD0000/MS/{0}/'.format(d)
        use_filepaths.append(filepath)
    return use_filepaths


if __name__ == '__main__':
    # embsGLS = getGLSEmbs()
    # embsGLS = emb_handler.updateEmbsMulti(embsGLS)
    # iterAvgCurvesGLS(embsGLS)
    # calcAvgParamsGLS(embsGLS)
    # checkScalesGLS(embsGLS)

    embsMS = getMSEmbs()
    embsMS = emb_handler.updateEmbsMulti(embsMS)
    iterAvgCurvesMS(embsMS)
    calcAvgParamsMS(embsMS)
    # checkScalesMS(embsMS)
