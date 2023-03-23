'''
Created on Feb 16, 2017

@author: Renat
Check parameter distributions, compare to histogram to gaussian and cumulative between controls and RNAi 
'''

from RNAiClass import RNAiClass, ControlClass
from varLookup import FOLDER_IN
from myFigure import myFigure
import numpy as np
from scipy.optimize import curve_fit
import csv
from params_use_assign import get_params_use
PARAM_NAMES_USE_GLS, PARAM_NAMES_USE_MS = get_params_use()

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def showMSParamDist(rnai, pName, fig=None):
    print('plotting {pn} distribution'.format(pn=pName))
    ps = rnai.paramsEmbsGLS[pName]
#     ps = []
#     for emb in embs:
#         if not np.isnan(emb.params[pName]):
#             ps.append(emb.params[pName])
    if not fig: fig = myFigure()
    bins = fig.hist(ps, bins=20, normed=True)
    y, x = bins[0], bins[1]
    x = [(x[i]+x[i+1])/2. for i in range(x.size-1)]
    n = len(x)                          #the number of data
    mean = sum(x*y)/n                   #note this correction
    sigma = sum(y*(x-mean)**2)/n        #note this correction
     
    def gaus(x,a,x0,sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))
     
    try:
        popt,pcov = curve_fit(gaus,x,y,p0=[1,mean,sigma])
        x = np.arange(0.8*min(ps),1.2*max(ps), 0.01*(1.2*max(ps)-0.8*min(ps)))
        fig.plot(x, gaus(x,*popt), color='k')
    except: print('could not fit {pn}'.format(pn=pName))
    fig.title(pName)
    return fig

if __name__ == '__main__':
#     control = ControlClass()
#     control.refresh_params_data()
    

    #     embs = control.MSUseEmbs
    #     with open(FOLDER_IN+'Automated_analysis/MS_params.csv', 'wb') as csvfile:
    #         spamwriter = csv.writer(csvfile, delimiter=',',
    #                             quotechar='|', quoting=csv.QUOTE_MINIMAL)
    #         spamwriter.writerow(['date','emb#']+PARAM_NAMES_USE_MS)
    #         for emb in embs:
    #             ps = []
    #             for pN in PARAM_NAMES_USE_MS:
    #                 ps.append(emb.params[pN])
    #             spamwriter.writerow([emb.date, emb.label]+ps)
        
    '''writes all GLS parameters to csv file- this is used to pull out automated nuclear counts
    to be added into manual binary profiles '''
        
    with open(FOLDER_IN+'Automated_analysis/testGLS_params.csv', 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow(['date','RNAi','emb#']+PARAM_NAMES_USE_GLS)
#             useList = [1,2,4,10]
        for i in range(1,440):
            r1 = RNAiClass(i) 
            print r1.RNAi
            r1.refresh_params_data()
#                 r1.load_params_data()

            embs = r1.GLSUseEmbs
            print embs

            for emb in embs:
                ps = []
                for pN in PARAM_NAMES_USE_GLS:
                    ps.append(emb.params[pN])
                spamwriter.writerow([emb.date, emb.RNAi, emb.label]+ps)
#     for pN in PARAM_NAMES_USE_MS:
#         fig = showMSParamDist(embs, pN)
#     fig.show()
    control = ControlClass()
    control.refresh_params_data()
#     embs = control.GLSUseEmbs
#     with open(FOLDER_IN+'Automated_analysis/MS_params.csv', 'wb') as csvfile:
#         spamwriter = csv.writer(csvfile, delimiter=',',
#                             quotechar='|', quoting=csv.QUOTE_MINIMAL)
#         spamwriter.writerow(['date','emb#']+PARAM_NAMES_USE_MS)
#         for emb in embs:
#             ps = []
#             for pN in PARAM_NAMES_USE_MS:
#                 ps.append(emb.params[pN])
#             spamwriter.writerow([emb.date, emb.label]+ps)
#     embs = control.GLSUseEmbs
#     with open(FOLDER_IN+'Automated_analysis/GLS_params.csv', 'wb') as csvfile:
#         spamwriter = csv.writer(csvfile, delimiter=' ',
#                             quotechar='|', quoting=csv.QUOTE_MINIMAL)
#         spamwriter.writerow(['date','emb#']+PARAM_NAMES_USE_MS)
#         for emb in embs:
#             ps = []
#             for pN in PARAM_NAMES_USE_MS:
#                 ps.append(emb.params[pN])
#             spamwriter.writerow([emb.date, emb.label]+ps)
    for pN in ['MoI1Rscale_GLS']:
        fig = showMSParamDist(control, pN)
    fig.show()
