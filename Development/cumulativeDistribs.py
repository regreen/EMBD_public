'''
Created on Sep 8, 2016

@author: Admin
'''

from RNAiClass import RNAiClass, ControlClass
from myFigure import myFigure
from scipy.optimize.minpack import curve_fit
from scipy import stats
import csv
from sklearn.decomposition import PCA as sklearnPCA
from avgCalculator import getGLSEmbs, getMSEmbs
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from groupAssesment import getRNAiGroupsPheno
import numpy as np
import matplotlib.pyplot as plt
from myFunc import readCSV
from params_use_assign import get_params_use
from db_utils_embryos import get_all_params, initialize_db
from emb_handler import printLog
PARAM_NAMES_USE_GLS, PARAM_NAMES_USE_MS = get_params_use()


def plotCum(data, title, fig=None, color='k'):
    if fig is None: fig = myFigure()
    if data.size > 1:
        data = np.sort(data)
        y = np.arange(data.size)
        fig.plot(data, 1. * y / data.size, color=color)
        fig.title(title)
    return fig


def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def show_param_dist(p_name, strain, rnai=None, fig=None):
    """
    plots parameter distribution and fits it with gaussian
    :param p_name: parameter name
    :param strain: strain
    :param rnai: None (all rnais) or 0 (controls)
    :param fig: myFigure object
    :return: myFigure object
    """
    print('plotting {pn} distribution for {s}'.format(pn=p_name, s=strain))
    con, curs = initialize_db()
    ps_rnai = get_all_params(curs, p_name, strain, rnai=rnai).tolist()
    if fig is None:
        fig = myFigure()
    bins = fig.hist(ps_rnai, bins=100, normed=True)
    y, x = bins[0], bins[1]
    x = [(x[i] + x[i + 1]) / 2. for i in range(x.size - 1)]
    n = len(x)  # the number of data
    # mean = sum(x * y) / n  # note this correction
    # sigma = sum(y * (x - mean) ** 2) / n  # note this correction
    mean = np.mean(ps_rnai)
    sigma = np.std(ps_rnai)


    def gaus(x, a, x0, sigma):
        return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

    try:
        popt, pcov = curve_fit(gaus, x, y, p0=[1, mean, sigma])
        x = np.arange(0.8 * min(ps_rnai), 1.2 * max(ps_rnai), 0.01 * (1.2 * max(ps_rnai) - 0.8 * min(ps_rnai)))
        fig.plot(x, gaus(x, *popt), color='k')
    except:
        print('could not fit {pn}'.format(pn=p_name))
    fig.title(p_name)
    return fig


def getAllPCC(RNAi, c):
    folder = 'Z:/Automated_analysis/PAD/'
    with open(folder + 'EMBD_paramPCC_GLS.csv', 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        pns = list(PARAM_NAMES_USE_GLS)
        #         pns.reverse()
        spamwriter.writerow([''] + pns)
        im = np.zeros((len(PARAM_NAMES_USE_GLS), len(PARAM_NAMES_USE_GLS)))
        for i in range(len(PARAM_NAMES_USE_GLS)):
            pccList = []
            for j in range(len(PARAM_NAMES_USE_GLS)):
                pN1 = PARAM_NAMES_USE_GLS[i]
                pN2 = PARAM_NAMES_USE_GLS[j]
                pAll = []
                for r in RNAi:
                    p1 = r.paramsGLS[pN1]
                    p2 = r.paramsGLS[pN2]
                    if not np.isnan(p1) and not np.isnan(p2): pAll.append((p1, p2))
                pAll = np.array(pAll)
                if pAll.size > 4:
                    pccList.append(stats.pearsonr(pAll[:, 0], pAll[:, 1])[0])
                else:
                    pccList.append(0.)
            im[i] = np.array(pccList)
            # print ('{0} correlations = {1}'.format(pN1, pccList))
            #             pccList.reverse()
            spamwriter.writerow([pN1] + pccList)
    with open(folder + 'EMBD_paramPCC_MS.csv', 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        pns = list(PARAM_NAMES_USE_MS)
        spamwriter.writerow([''] + pns)
        im = np.zeros((len(PARAM_NAMES_USE_MS), len(PARAM_NAMES_USE_MS)))
        for i in range(len(PARAM_NAMES_USE_MS)):
            pccList = []
            for j in range(len(PARAM_NAMES_USE_MS)):
                pN1 = PARAM_NAMES_USE_MS[i]
                pN2 = PARAM_NAMES_USE_MS[j]
                pAll = []
                for r in RNAi:
                    p1 = r.paramsMS[pN1]
                    p2 = r.paramsMS[pN2]
                    if not np.isnan(p1) and not np.isnan(p2): pAll.append((p1, p2))
                pAll = np.array(pAll)
                if pAll.size > 4:
                    pccList.append(stats.pearsonr(pAll[:, 0], pAll[:, 1])[0])
                else:
                    pccList.append(0.)
            im[i] = np.array(pccList)
            # print ('{0} correlations = {1}'.format(pN1, pccList))
            spamwriter.writerow([pN1] + pccList)


# fig = myFigure()
#     fig.imshow(im, bw=False, colorbar=True)
#     fig.show()


def load_pcc():
    folder = 'Z:/Automated_analysis/PAD/'
    file_GLS = folder + 'EMBD_paramPCC_GLS.csv'
    data = readCSV(file_GLS, ',')[1:]
    tmp = []
    for i in range(len(data)):
        tmp.append((data[i][0] + '_GLS', np.max(np.abs(data[i][1:i + 1] + data[i][i + 2:]))))
    file_MS = folder + 'EMBD_paramPCC_MS.csv'
    data = readCSV(file_MS, ',')[1:]
    for i in range(len(data)):
        tmp.append((data[i][0] + '_MS', np.max(np.abs(data[i][1:i + 1] + data[i][i + 2:]))))
    tmp.sort(key=lambda x: x[1])
    for t in tmp:
        print(t)


def get_pca(RNAi):
    p_matrix = []  # this is a matrix of parameter values
    for r in RNAi:  # for each RNAi condition
        r.setOrigin()  # sets the origin
        p_values = []  # list to be filled with normalized parameter values for each RNAi condition
        for pN in PARAM_NAMES_USE_GLS:
            raw_pv = r.paramsGLS[pN]  # raw parameter value
            if np.isnan(raw_pv):
                raw_pv = 0  # sets nan parameters to 0
            pv = (raw_pv - r.controlOrigin[pN + '_GLS']) / r.rnaiStd[
                pN + '_GLS']  # calculates the normalized parameter value

            p_values.append(pv)  # appends normalized parameter values to p values list
        for pN in PARAM_NAMES_USE_MS:
            raw_pv = r.paramsMS[pN]
            if np.isnan(raw_pv):
                raw_pv = 0
            pv = (raw_pv - r.controlOrigin[pN + '_MS']) / r.rnaiStd[pN + '_MS']
            p_values.append(pv)
        p_matrix.append(p_values)  # appends p_values to P_matrix (needed for PCA analysis)
    p_list = [p + '_GLS' for p in PARAM_NAMES_USE_GLS] + [p + '_MS' for p in PARAM_NAMES_USE_MS]
    sklearn_pca = sklearnPCA(n_components=9)  # specifies the number of components to look at i.e. (n_components=2)
    X_pca = sklearn_pca.fit_transform(p_matrix)  # returns array-like, shape (n_samples, n_components)

    tot_alpha = np.sum(np.abs(sklearn_pca.components_), axis=0)
    #  tot_alpha = np.abs(sklearn_pca.mean)
    ta_ind = np.argsort(tot_alpha)
    # for ind in ta_ind:
    #     print p_list[ind], tot_alpha[ind]
    print list(np.array(p_list)[ta_ind])
    print list(tot_alpha[ta_ind])
    p_matrix = np.array(p_matrix)
    return p_matrix


#     comp_list = []
#     for comp in sklearn_pca.components_:
#         print 'start next component'
#         comp_ind = np.argsort(np.abs(comp))
#         for ind in comp_ind[::-1]:
#             print p_list[ind], comp[ind]

#     print sklearn_pca.components_
#     print sklearn_pca.explained_variance_ratio_

def get_pca_controlGLS(RNAi):
    '''the purpose of this function is to identify the parameters that give the most scatter in controls (least reliable parameters)'''
    GLSembryos = getGLSEmbs()
    raw_pv_matrix = []  # a list of all parameter values for all embryos (not normalized)
    #     p_matrix = []  # this is a matrix of parameter values
    #     p_values = []  # list to be filled with normalized parameter values for each RNAi condition
    r = RNAi[0]
    for emb in GLSembryos:
        raw_pv_values = []  # a list of all parameter values for 1 embryo
        for pN in PARAM_NAMES_USE_GLS:
            raw_pv = (emb.params[pN] - r.controlOrigin[pN + '_GLS']) / r.rnaiStd[pN + '_GLS']  # raw parameter value
            if np.isnan(raw_pv):
                raw_pv = 0  # sets nan parameters to 0
            raw_pv_values.append(raw_pv)
        raw_pv_matrix.append(raw_pv_values)
    raw_pv_matrix = np.array(raw_pv_matrix)
    mean_rpv = np.mean(raw_pv_matrix, axis=0)  # finds average for later normalization
    std_rpv = np.std(raw_pv_matrix, axis=0)  # finds std for later normalization

    p_matrix = raw_pv_matrix  # a list of all normalized parameters for all embryos
    # for row in raw_pv_matrix:
    #     corr_pvs = (row - mean_rpv)/std_rpv#calculates normalized values
    #     copy=np.nan_to_num(corr_pvs)#converts nans to 0
    #     p_matrix.append(copy)

    p_list = [p + '_GLS' for p in PARAM_NAMES_USE_GLS]

    #     np.nan_to_num(p_matrix, False)
    sklearn_pca = sklearnPCA(n_components=2)  # specifies the number of components to look at
    X_pca = sklearn_pca.fit_transform(p_matrix)  # returns array-like, shape (n_samples, n_components)

    tot_alpha = np.sum(np.abs(sklearn_pca.components_), axis=0)  # adds up the rank (alpha) value
    ta_ind = np.argsort(tot_alpha)
    # for ind in ta_ind:
    #     print p_list[ind], tot_alpha[ind]
    p_matrix = np.array(p_matrix)
    return p_matrix


#     comp_list = []
#     for comp in sklearn_pca.components_:
#         print 'start next component'
#         comp_ind = np.argsort(np.abs(comp))
#         for ind in comp_ind[::-1]:
#             print p_list[ind], comp[ind]

#     print sklearn_pca.components_
#     print sklearn_pca.explained_variance_ratio_

def get_pca_controlMS(RNAi):  # IN PROGRESS!!!###
    MSembryos = getMSEmbs()
    raw_pv_matrix = []  # a list of all parameter values for all embryos (not normalized)
    #     p_matrix = []  # this is a matrix of parameter values
    #     p_values = []  # list to be filled with normalized parameter values for each RNAi condition
    r = RNAi[0]
    for emb in MSembryos:
        raw_pv_values = []  # a list of all parameter values for 1 embryo
        for pN in PARAM_NAMES_USE_MS:
            raw_pv = (emb.params[pN] - r.controlOrigin[pN + '_MS']) / r.rnaiStd[pN + '_MS']  # raw parameter value
            if np.isnan(raw_pv):
                raw_pv = 0  # sets nan parameters to 0
            raw_pv_values.append(raw_pv)
        raw_pv_matrix.append(raw_pv_values)
    raw_pv_matrix = np.array(raw_pv_matrix)
    mean_rpv = np.mean(raw_pv_matrix, axis=0)  # finds average for later normalization
    std_rpv = np.std(raw_pv_matrix, axis=0)  # finds std for later normalization

    p_matrix = raw_pv_matrix  # a list of all normalized parameters for all embryos
    # for row in raw_pv_matrix:
    #     corr_pvs = (row - mean_rpv)/std_rpv
    #     copy=np.nan_to_num(corr_pvs)
    #     p_matrix.append(copy)

    p_list = [p + '_MS' for p in PARAM_NAMES_USE_MS]

    #     np.nan_to_num(p_matrix, False)
    sklearn_pca = sklearnPCA(n_components=2)  # specifies the number of components to look at
    X_pca = sklearn_pca.fit_transform(p_matrix)  # returns array-like, shape (n_samples, n_components)

    tot_alpha = np.sum(np.abs(sklearn_pca.components_), axis=0)
    ta_ind = np.argsort(tot_alpha)
    # for ind in ta_ind:
    #     print p_list[ind], tot_alpha[ind]
    p_matrix = np.array(p_matrix)
    return p_matrix


#     comp_list = []
#     for comp in sklearn_pca.components_:
#         print 'start next component'
#         comp_ind = np.argsort(np.abs(comp))
#         for ind in comp_ind[::-1]:
#             print p_list[ind], comp[ind]

#     print sklearn_pca.components_
#     print sklearn_pca.explained_variance_ratio_


def get_lda(rnai_matrix, cont_gls_matrix, cont_ms_matrix):
    '''compares control GLS, MS with RNAi GLS, MS and performs LDA between control and RNAi conditions to identify the parameters that give the best separation between these two conditions'''

    rnai_gls = rnai_matrix[:, :len(PARAM_NAMES_USE_GLS)]
    rnai_ms = rnai_matrix[:, len(PARAM_NAMES_USE_GLS):]

    gls_matrix = np.concatenate((rnai_gls, cont_gls_matrix))
    ms_matrix = np.concatenate((rnai_ms, cont_ms_matrix))

    gls_label = np.zeros(gls_matrix.shape[0])
    gls_label[:rnai_gls.shape[0]] = 1

    ms_label = np.zeros(ms_matrix.shape[0])
    ms_label[:rnai_ms.shape[0]] = 1

    sklearn_lda_gls = LDA()
    sklearn_lda_ms = LDA()
    X_lda_sklearn_gls = sklearn_lda_gls.fit_transform(gls_matrix, gls_label)
    X_lda_sklearn_ms = sklearn_lda_ms.fit_transform(ms_matrix, ms_label)

    p_list_gls = [p + '_GLS' for p in PARAM_NAMES_USE_GLS]
    ev = np.ravel(sklearn_lda_gls.scalings_)
    ev /= np.linalg.norm(ev)
    comp_ind = np.argsort(np.abs(ev))
    print('-------------------------GLS-------------------------')
    for ind in comp_ind[::-1]:
        print p_list_gls[ind], ev[ind]
    print list(np.array(PARAM_NAMES_USE_GLS)[comp_ind][::-1])
    print list(np.abs(ev[comp_ind])[::-1])

    p_list_ms = [p + '_MS' for p in PARAM_NAMES_USE_MS]
    ev = np.ravel(sklearn_lda_ms.scalings_)
    ev /= np.linalg.norm(ev)
    comp_ind = np.argsort(np.abs(ev))
    print('-------------------------MS-------------------------')
    for ind in comp_ind[::-1]:
        print p_list_ms[ind], ev[ind]
    print list(np.array(PARAM_NAMES_USE_MS)[comp_ind][::-1])
    print list(np.abs(ev[comp_ind])[::-1])
    fig = plot_scikit_lda(X_lda_sklearn_gls, gls_label)
    fig.title('GLS')
    fig = plot_scikit_lda(X_lda_sklearn_ms, ms_label)
    fig.title('MS')
    fig.show()


def get_lda_by_groups():
    '''takes in our manually curated groups as RNAi class objects and performs LDA on the groups to identify the parameters that give the best separation between these groups'''
    rlist_test, rlist_train, simGroups_test, simGroups_train, labels, labels_test, defect_label, defect_label_test = getRNAiGroupsPheno()
    all_matrix = []
    g_label = []  # group labels
    count = 0
    for group in rlist_train:  # this is a list of numbers representing RNAi conditions i.e. [21,302,26,5]
        gr_objects = []  # RNAi objects corresponding to the list of numbers in group
        for i in group:  # iterates through each RNAi # within group and adds to gr_objects
            r = RNAiClass(i, verb=True)
            gr_objects.append(r)
        p_matrix = get_pca(
            gr_objects)  # takes a list of RNAi objects and returns an assembled matrix of parameters for the group
        all_matrix.append(p_matrix)  # appends to all_matrix
        g_label.append(np.full(p_matrix.shape[0], count))  # generates a label array for each group- increases count
        count += 1

    all_matrixC = np.concatenate(all_matrix)
    g_labelC = np.concatenate(g_label)

    sklearn_lda_groups = LDA(n_components=8)
    X_lda_sklearn_groups = sklearn_lda_groups.fit_transform(all_matrixC, g_labelC)
    wt = sklearn_lda_groups.transform(np.zeros_like(all_matrixC))[0]
    p_list = [p + '_GLS' for p in PARAM_NAMES_USE_GLS] + [p + '_MS' for p in PARAM_NAMES_USE_MS]
    evs = sklearn_lda_groups.scalings_.T
    w = sklearn_lda_groups.explained_variance_ratio_
    print('LDA explained variance = {}'.format(w))
    ev_sum = np.zeros_like(evs[0])
    for i in range(w.size):
        ev = evs[i]
        ev /= np.linalg.norm(ev)
        ev *= w[i]
        ev_sum += np.abs(ev)
        # print ev
    comp_ind = np.argsort(ev_sum)
    # for ind in comp_ind[::-1]:
    #     print p_list[ind], ev[ind]
    p_gls = []
    p_ms = []
    v_gls = []
    v_ms = []
    for i in range(comp_ind.size):
        p_name = np.array(p_list)[comp_ind][::-1][i]
        if p_name[-3:] == 'GLS':
            p_gls.append(p_name[:-4])
            v_gls.append(np.abs(ev_sum[comp_ind])[::-1][i])
        else:
            p_ms.append(p_name[:-3])
            v_ms.append(np.abs(ev_sum[comp_ind])[::-1][i])

    print('GLS params: {0}\nvals: {1}'.format(p_gls, v_gls))
    print('MS params: {0}\nvals: {1}'.format(p_ms, v_ms))
    # print list(np.array(p_list)[comp_ind][::-1])
    # print list(np.abs(ev_sum[comp_ind])[::-1])

    # '''plot histogram for 1D data'''
    # fig = plot_scikit_lda(X_lda_sklearn_groups, g_labelC)
    # fig.title('LDA by groups')
    # fig.show()
    '''plot data by 2 LDA components'''
    plt.figure()
    colors = ['sienna', 'firebrick', 'orangered', 'gold', 'y', 'green', 'turquoise', 'c', 'dodgerblue', 'b', 'm',
              'purple', 'gray', 'k', 'rosybrown', 'darkkhaki', 'lightslategray', 'orchid', 'darkmagenta']
    lw = 2
    X_r2 = X_lda_sklearn_groups
    y = g_labelC
    target_names = labels
    for color, i, target_name in zip(colors, range(len(labels)), target_names):
        plt.scatter(X_r2[y == i, 0], X_r2[y == i, 1], alpha=1.0, color=color,
                    label=target_name)
    # plt.scatter(X_r2[0],X_r2[1])
    plt.scatter(wt[0], wt[1], color='k',
                label='WT', s=100, marker='*')
    plt.legend(loc='lower right', shadow=False, scatterpoints=1)
    plt.title('LDA of groups')
    plt.show()


def plot_scikit_lda(X, y):
    fig = myFigure()
    fig.hist(X[y == 0], alpha=0.5, label='WT', normed=True)
    fig.hist(X[y == 1], alpha=0.5, label='RNA', normed=True)
    fig.legend()
    return fig



if __name__ == '__main__':
    # load_pcc()
    # for pn in PARAM_NAMES_USE_GLS:
    #     fig = show_param_dist(pn, 'GLS')
    #     fig = show_param_dist(pn, 'GLS', rnai=0, fig=fig)
    #     fig.show()
    # for pn in PARAM_NAMES_USE_MS:
    #     fig = show_param_dist(pn, 'MS')
    #     fig = show_param_dist(pn, 'MS', rnai=0, fig=fig)
    #     fig.show()
    # printLog('______________START cumulativeDistribs_____________________')
    # # c = ControlClass(verb=True)
    # RNAi = []
    # for i in range(10, 25):
    #     try:
    #         print('loading rnai={}'.format(i))
    #         r = RNAiClass(i, verb=True)
    #         RNAi.append(r)
    #     except Exception, e:
    #         printLog('!!!!!!!!!ERROR OCCURED RNAi {0} :{1}'.format(i,str(e)))
    # # get_pca(RNAi)
    # get_lda(get_pca(RNAi), get_pca_controlGLS(RNAi), get_pca_controlMS(RNAi))
    # print "All Finished"

    get_lda_by_groups()

    #     getAllPCC(RNAi, c)
    #     KST = {}
    #     for pN in PARAM_NAMES:
    #         if isinstance(c.paramsEmbsGLS[pN],list) and len(c.paramsEmbsGLS[pN])>0:
    #             p = np.array([])
    #             for r in RNAi:
    #                 rps = np.array(r.paramsEmbsGLS[pN])
    #                 nans, x = nan_helper(rps)
    #                 rps = rps[~nans]
    #                 if rps.size>0: p = np.concatenate((p,rps))
    #             cPars = np.array(c.paramsEmbsGLS[pN])
    #             if cPars.size>1 and p.size>1:
    #                 print('{n} quality (K-S test)={kst}'.format(n=pN, kst=stats.ks_2samp(cPars,p)[0]))
    #                 KST[pN+'_GLS']=stats.ks_2samp(cPars,p)[0]
    #             else:
    #                 KST[pN+'_GLS']=0
    #             fig = plotCum(cPars, pN)
    #             fig = plotCum(p, pN, fig, 'r')
    #             if cPars.size>1 and p.size>1:fig.title('{n} K-S={ks}'.format(n=pN, ks=stats.ks_2samp(cPars,p)[0]))
    #         else:
    #             print('parameter is nan {p}'.format(p=pN))
    #             KST[pN+'_GLS']=0
    #
    #         if isinstance(c.paramsEmbsMS[pN],list) and len(c.paramsEmbsMS[pN])>0:
    #             p = np.array([])
    #             for r in RNAi:
    #                 rps = np.array(r.paramsEmbsMS[pN])
    #                 nans, x = nan_helper(rps)
    #                 rps = rps[~nans]
    #                 if rps.size>0: p = np.concatenate((p,rps))
    #             cPars = np.array(c.paramsEmbsMS[pN])
    #             fig = plotCum(cPars, pN)
    #             if cPars.size>1 and p.size>1:
    #                 print('{n} quality (K-S test)={kst}'.format(n=pN, kst=stats.ks_2samp(cPars,p)[0]))
    #                 KST[pN+'_MS']=stats.ks_2samp(cPars,p)[0]
    #             else:
    #                 print('{n} can not calculate K-S test control#={0}, rnai#={1}'.format(cPars.size, p.size,n=pN))
    #                 KST[pN+'_MS']=0
    #             fig=plotCum(p, pN, fig, 'r')
    #             fig.title('{n} K-S={ks}'.format(n=pN, ks=stats.ks_2samp(cPars,p)[0]))
    #         else:
    #             print('parameter is nan {p}'.format(p=pN))
    #             KST[pN+'_MS']=0
    #     print('n controls MS = {0}; n controls GLS = {1}'.format(len(c.paramsEmbsMS['aR']),len(c.paramsEmbsGLS['aR'])))
    #     print('KST=',KST)
    # #     plt.show()
    #
