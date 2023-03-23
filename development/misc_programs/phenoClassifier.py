'''
Created on Feb 20, 2017

@author: Renat
machine learning method trained with embryos scored by manual descriptor
to identify params that contribute to common phenotypes 
'''
import embdFunc, os, myFunc, csv
from varLookup import *
from avgCalculator import getMSEmbs, getGLSEmbs
from sklearn.neural_network import MLPClassifier
from sklearn import preprocessing, metrics, tree, svm
import numpy as np
from RNAiClass import ControlClass
from myFunc import readCSV


def readEmbsFromCSV(fileName, nEmbs=None, skip=0):
    data = readCSV(fileName, ',')[1 + skip:]  # skip first line
    embs = []
    if nEmbs is None or nEmbs > len(data): nEmbs = len(data)
    foldersMS, foldersGLS = [], []
    for d in data[:nEmbs]:
        if d[2] == 'EMBD0000':
            folder = FOLDER_IN + 'cropped/{rna}/{strain}/{date}/{emb}/'.format(strain=d[0], date=d[1], rna=d[2],
                                                                               emb=d[3])
        else:
            folder = FOLDER_IN + 'cropped/{rna}/{strain}/{emb}/'.format(strain=d[0], date=d[1], rna=d[2], emb=d[3])
        if not os.path.exists(folder):
            print('embryo from csv file does not exist {f}'.format(f=folder))
            raise
        if d[0] == 'MS':
            foldersMS.append(folder)
        else:
            foldersGLS.append(folder)
    pool = mp.Pool(processes=nCores)
    if len(foldersMS) > 0:
        embsMS = pool.map(loadEmbs, foldersMS)
    else:
        embsMS = []
    if len(foldersGLS) > 0:
        embsGLS = pool.map(loadEmbs, foldersGLS)
    else:
        embsGLS = []
    pool.close()
    pool.join()
    return embsMS + embsGLS

def initiate():
    np.random.seed(71076)
#     np.random.seed(12483)
#     np.random.seed(40794)
    disc = '06' 
    strain = 'MS'
    fileName = FOLDER_IN + 'MachineLearningTrainingFiles/{disc}_{strain}_trainingset.csv'.format(disc=disc, strain=strain)
    embs = embdFunc.readEmbsFromCSV(fileName)
    fileName5 = FOLDER_IN + 'MachineLearningTrainingFiles/{disc}_{strain}_trainingset.csv'.format(disc='05', strain=strain)
    embs5 = embdFunc.readEmbsFromCSV(fileName5)
    embs = embdFunc.updateEmbsMulti(embs)
    embs5 = embdFunc.updateEmbsMulti(embs)
    if strain=='MS':
        p_use = PARAM_NAMES_USE_MS
        embsC = getMSEmbs()
    else:
        p_use = PARAM_NAMES_USE_GLS
        embsC = getGLSEmbs()
    inp = []
    inpTD = []
    res = []
    labels = []
    for emb in embs+embsC+embs5:
        x = []
        xtd = []
        for p in p_use:
            val = emb.params[p]
            if p[-5:]=='td_MS':
                if np.isnan(val): val=31
                xtd.append(val)
            else:
                if np.isnan(val):
                    print('{date} {emb} parameter {p}=NaN'.format(date=emb.date, emb=emb.label, p=p))
                    val=0
                x.append(val)
        inp.append(x)
        inpTD.append(xtd)
        res.append(int(emb in embs))
        labels.append(emb.label)
    inp = np.array(inp)
    inpTD = np.array(inpTD)
    res = np.array(res)
    labels = np.array(labels)
    inp = preprocessing.StandardScaler().fit_transform(inp)
    inpTD = preprocessing.MinMaxScaler().fit_transform(inpTD)
    inpAll = np.hstack((inp,inpTD))
    clf = MLPClassifier(solver='lbfgs', alpha=1., hidden_layer_sizes=(15), random_state=1)
#     clf = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(16,8,4,2), random_state=1)
    inds = np.arange(res.size)
    np.random.shuffle(inds)
    inpAll = inpAll[inds]
    res = res[inds]
    labels = labels[inds]
    testSize = int(0.5*inds.size)
    clf.fit(inpAll[testSize:], res[testSize:])
    pred = []
    for i in range(testSize):
        pred.append(clf.predict(inpAll[i])[0])
    pred = np.array(pred)
    f1Cl = metrics.f1_score(res[:testSize],pred, average='binary')
    randBi = np.random.choice(2,testSize)
    f1Rnd = metrics.f1_score(res[:testSize],randBi, average='binary')
    f10 = metrics.f1_score(res[:testSize],np.zeros(testSize), average='binary')
    f11 = metrics.f1_score(res[:testSize],np.ones(testSize), average='binary')
    print('Classifier F1 metric performance: {cl}\nrandom: {r}\nall 0: {z}\nall 1: {o}'.format(cl=f1Cl, r=f1Rnd, z=f10, o=f11))
    print('precision = {p}\nrecall = {r}'.format(p=metrics.precision_score(res[:testSize],pred),r=metrics.recall_score(res[:testSize],pred)))
#         print('{l} is {r} predicted {pred}'.format(pred=clf.predict(inpAll[i])[0], r=res[i], l=labels[i]))

def multiClassMS():
#     np.random.seed(71076)
    c = ControlClass()
    pMS = c.paramsMS
    pMSstd = c.paramsMSstd
    
    fileParams = FOLDER_IN + 'MachineLearningTrainingFiles/MS_params_use.csv'
    fileName = FOLDER_IN + 'MachineLearningTrainingFiles/all_phenotypes_MS.csv'
    params = {}
    if not os.path.exists(fileParams):
        embs = embdFunc.readEmbsFromCSV(fileName)
        params['date']=[]
        params['rnai']=[]
        params['emb#']=[]
        for key in PARAM_NAMES_USE_MS:
            params[key]=[]
        for emb in embs:
            print(emb.label, emb.RNAi)
            params['date'].append(emb.date)
            params['rnai'].append(emb.RNAi)
            params['emb#'].append(emb.label)
            for key in PARAM_NAMES_USE_MS:
                val = emb.params[key]
                if key[-5:]=='td_MS' and np.isnan(val):val=31.
                elif key=='tMove' and np.isnan(val):val=31.
                elif key=='movement' and val: val=1
                elif key=='movement' and not val: val=0
                elif np.isnan(val): val=0.
                params[key].append(val)
        with open(fileParams, 'w') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['date','rnai','emb#']+PARAM_NAMES_USE_MS)
            for i in range(len(embs)):
                row = []
                for key in ['date','rnai','emb#']+PARAM_NAMES_USE_MS:
                    row.append(params[key][i])
                writer.writerow(row)
        nEmbs = len(embs)
        del embs
    
    data = myFunc.readCSV(fileParams, ',')
    for i in range(len(data[0])):
        params[data[0][i]] = []
    for d in data[1:]:
        for i in range(len(data[0])):
            key = data[0][i]
            val = d[i]
            if key in PARAM_NAMES_USE_MS:
#                 m = pMS[key]
#                 sig = pMSstd[key]
#                 if not np.isnan(m) and not np.isnan(sig) and sig>0: val = (val-m)/sig
#                 elif key!='movement': print('Can not normalize {k}'.format(k=key))
                params[key].append(val)
            if key in ['date','rnai','emb#']: params[key].append(val)
    for i in range(len(data[0])):
        params[data[0][i]] = np.array(params[data[0][i]])
    nEmbs = len(data)-1
    del data
    data = myFunc.readCSV(fileName, ',')
    Dn = data[0][4:]
    for j in Dn:
        params['D{d}'.format(d=j)]=[]
    data = data[1:]
    if len(data)!=nEmbs: raise 'different number of embryos {n1} vs {n2}'.format(n1 = len(data), n2=nEmbs)
    y = []
    X = []
    Xmm = []
    for i in range(len(data)):
        d = data[i]
        if params['date'][i]!=d[1] or params['rnai'][i]!=d[2] or params['emb#'][i]!=d[3]:
            print('wrong embryo {i}: {d1} {r1} {e1} vs {d}'.format(i=i, d1=params['date'][i], r1=params['rnai'][i], e1=params['emb#'][i], d=d[1:4]))
            raise 
        x, xmm = [], []
        for pn in PARAM_NAMES_USE_MS:
            if pn[-5:]!='td_MS' and pn!='tMove' and pn!='movement': x.append(params[pn][i])
            else: xmm.append(params[pn][i])
        X.append(x)
        Xmm.append(xmm)
        y.append(d[4:])
    y = np.array(y).astype(np.int)
#     ytmp = np.max(y[:,3:5],axis=1)
#     y[:,3]=ytmp
#     y[:,4]=ytmp
#     dn=12
#     y=y[:,dn].reshape(-1,1)
    X = np.array(X)
    Xmm = np.array(Xmm)
    X = preprocessing.StandardScaler().fit_transform(X)
#     X = preprocessing.MinMaxScaler().fit_transform(X)
    Xmm = preprocessing.MinMaxScaler().fit_transform(Xmm)
    X = np.hstack((X,Xmm))
#     clfs = [MLPClassifier(activation='logistic', solver='lbfgs', alpha=0.1, hidden_layer_sizes=(15), random_state=i) for i in range(50)]
    clfs = [svm.SVC(kernel='rbf', c=1, gamma=1, random_state=i) for i in range(50)]
#     clfs = [tree.DecisionTreeClassifier(random_state=i) for i in range(50)]
#     clfs = [MLPClassifier(activation='logistic', solver='lbfgs', alpha=0.1, hidden_layer_sizes=(15,50), random_state=1)]
    inds = np.arange(y.shape[0])
    np.random.seed(71076)
    np.random.shuffle(inds)
    X = X[inds]
    y = y[inds]
    testSize = int(0.5*inds.size)
    Xtr = X[testSize:]
    ytr = y[testSize:]
    for i in range(len(clfs)):
        inds = np.arange(ytr.shape[0])
        np.random.shuffle(inds)
        clfs[i].fit(Xtr[inds], ytr[inds])
    pred = []
    for i in range(testSize):
        preds = np.array([clfs[j].predict(X[i])[0] for j in range(len(clfs))])
        p = np.around(np.mean(preds, axis=0))
        pred.append(p)
    pred = np.array(pred)
#     p = 1.*np.sum((pred== y[:testSize])*pred)/np.sum(pred)
#     r = 1.*np.sum((pred== y[:testSize])*pred)/np.sum(y[:testSize])
#     print('overall precision={p}, recall={r}, f1={f}'.format(p=p,r=r,f=2.*p*r/(p+r)))
    f1Cl = metrics.f1_score(y[:testSize],pred, average=None)
    randBi = np.random.choice(2,y[:testSize].shape).astype(np.float)
    f1Rnd = metrics.f1_score(y[:testSize],randBi, average=None)
    f11 = metrics.f1_score(y[:testSize],np.ones(y[:testSize].shape).astype(np.float), average=None)
    print('Classifier F1 metric performance: {cl}\nrandom: {r}\nall 1: {o}'.format(cl=f1Cl, r=f1Rnd, o=f11))
    print('precision = {p}\nrecall = {r}'.format(p=metrics.precision_score(y[:testSize],pred, average=None),r=metrics.recall_score(y[:testSize],pred, average=None)))
    print('random precision = {p}\nrecall = {r}'.format(p=metrics.precision_score(y[:testSize],randBi, average=None),r=metrics.recall_score(y[:testSize],randBi, average=None)))
    print('embryos show phenotypes in test: {tt}\nembryos predicted to show: {pr}'.format(tt=np.sum(y[:testSize],axis=0),pr=np.sum(pred,axis=0)))
    print('embryos show phenotype predicted right: {pr}'.format(pr=np.sum(pred*(y[:testSize]==pred),axis=0)))
    print('embryos show phenotype in training: {tr}'.format(tr=np.sum(y[testSize:],axis=0)))

def multiClassGLS():
    fileParams = FOLDER_IN + 'MachineLearningTrainingFiles/GLS_params_use.csv'
    fileName = FOLDER_IN + 'MachineLearningTrainingFiles/all_phenotypes_GLS.csv'
    params = {}
    if not os.path.exists(fileParams):
        embs = embdFunc.readEmbsFromCSV(fileName)
        params['date']=[]
        params['rnai']=[]
        params['emb#']=[]
        for key in PARAM_NAMES_USE_GLS:
            params[key]=[]
        for emb in embs:
            print(emb.label, emb.RNAi)
            params['date'].append(emb.date)
            params['rnai'].append(emb.RNAi)
            params['emb#'].append(emb.label)
            for key in PARAM_NAMES_USE_GLS:
                val = emb.params[key]
                if key[-6:]=='td_GLS' and np.isnan(val):val=31.
                elif key=='tMove' and np.isnan(val):val=31.
                elif key=='movement' and val: val=1
                elif key=='movement' and not val: val=0
                elif np.isnan(val): val=0.
                params[key].append(val)
        with open(fileParams, 'w') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['date','rnai','emb#']+PARAM_NAMES_USE_GLS)
            for i in range(len(embs)):
                row = []
                for key in ['date','rnai','emb#']+PARAM_NAMES_USE_GLS:
                    row.append(params[key][i])
                writer.writerow(row)
        nEmbs = len(embs)
        del embs
    else:
        data = myFunc.readCSV(fileParams, ',')
        for i in range(len(data[0])):
            params[data[0][i]] = []
        for d in data[1:]:
            for i in range(len(data[0])):
                params[data[0][i]].append(d[i])
        for i in range(len(data[0])):
            params[data[0][i]] = np.array(params[data[0][i]])
        nEmbs = len(data)-1
        del data
    data = myFunc.readCSV(fileName, ',')
    Dn = data[0][4:]
    for j in Dn:
        params['D{d}'.format(d=j)]=[]
    data = data[1:]
    if len(data)!=nEmbs: raise 'different number of embryos {n1} vs {n2}'.format(n1 = len(data), n2=nEmbs)
    y = []
    X = []
    Xmm = []
    for i in range(len(data)):
        d = data[i]
        if params['date'][i]!=d[1] or params['rnai'][i]!=d[2] or params['emb#'][i]!=d[3]:
            print('wrong embryo {i}: {d1} {r1} {e1} vs {d}'.format(i=i, d1=params['date'][i], r1=params['rnai'][i], e1=params['emb#'][i], d=d[1:4]))
            raise 
        x, xmm = [], []
        for pn in PARAM_NAMES_USE_GLS:
            if pn[-6:]!='td_GLS' and pn!='tMove' and pn!='movement': x.append(params[pn][i])
            else: xmm.append(params[pn][i])
        X.append(x)
        Xmm.append(xmm)
        y.append(d[4:])
    y = np.array(y).astype(np.int)
#     ytmp = np.max(y[:,3:5],axis=1)
#     y[:,3]=ytmp
#     y[:,4]=ytmp
#     dn=12
#     y=y[:,dn].reshape(-1,1)
    X = np.array(X)
    Xmm = np.array(Xmm)
    X = preprocessing.StandardScaler().fit_transform(X)
#     X = preprocessing.MinMaxScaler().fit_transform(X)
    Xmm = preprocessing.MinMaxScaler().fit_transform(Xmm)
    X = np.hstack((X,Xmm))
#     clfs = [MLPClassifier(activation='logistic', solver='lbfgs', alpha=0.1, hidden_layer_sizes=(15), random_state=i) for i in range(50)]
    clfs = [tree.DecisionTreeClassifier(random_state=i) for i in range(50)]
#     clfs = [MLPClassifier(activation='identity', solver='lbfgs', alpha=0.1, hidden_layer_sizes=(15,50), random_state=2)]
    inds = np.arange(y.shape[0])
    np.random.seed(71076)
    np.random.shuffle(inds)
    X = X[inds]
    y = y[inds]
    testSize = int(0.5*inds.size)
    Xtr = X[testSize:]
    ytr = y[testSize:]
    for i in range(len(clfs)):
        inds = np.arange(ytr.shape[0])
        np.random.shuffle(inds)
        clfs[i].fit(Xtr[inds], ytr[inds])
    pred = []
    for i in range(testSize):
        preds = np.array([clfs[j].predict(X[i])[0] for j in range(len(clfs))])
        p = np.around(np.mean(preds, axis=0))
        pred.append(p)
    pred = np.array(pred)
#     p = 1.*np.sum((pred== y[:testSize])*pred)/np.sum(pred)
#     r = 1.*np.sum((pred== y[:testSize])*pred)/np.sum(y[:testSize])
#     print('overall precision={p}, recall={r}, f1={f}'.format(p=p,r=r,f=2.*p*r/(p+r)))
    f1Cl = metrics.f1_score(y[:testSize],pred, average=None)
    randBi = np.random.choice(2,y[:testSize].shape).astype(np.float)
    f1Rnd = metrics.f1_score(y[:testSize],randBi, average=None)
    f11 = metrics.f1_score(y[:testSize],np.ones(y[:testSize].shape).astype(np.float), average=None)
    print('Classifier F1 metric performance: {cl}\nrandom: {r}\nall 1: {o}'.format(cl=f1Cl, r=f1Rnd, o=f11))
    print('precision = {p}\nrecall = {r}'.format(p=metrics.precision_score(y[:testSize],pred, average=None),r=metrics.recall_score(y[:testSize],pred, average=None)))
    print('random precision = {p}\nrecall = {r}'.format(p=metrics.precision_score(y[:testSize],randBi, average=None),r=metrics.recall_score(y[:testSize],randBi, average=None)))
    print('embryos show phenotypes in test: {tt}\nembryos predicted to show: {pr}'.format(tt=np.sum(y[:testSize],axis=0),pr=np.sum(pred,axis=0)))
    print('embryos show phenotype predicted right: {pr}'.format(pr=np.sum(pred*(y[:testSize]==pred),axis=0)))
    print('embryos show phenotype in training: {tr}'.format(tr=np.sum(y[testSize:],axis=0)))

# def updateAllEmbs():
#     fileName = FOLDER_IN + 'MachineLearningTrainingFiles/all_phenotipes_MS.csv'
#     for i in range(660, 1800, 50):
#         printLog('PhenoClassifier: MS loading lines {i}'.format(i=i))
#         embs = embdFunc.readEmbsFromCSV(fileName, nEmbs=50, skip=i)
#         del embs
#     embs = embdFunc.refreshEmbsMulti(embs)
#     fileName = FOLDER_IN + 'MachineLearningTrainingFiles/all_phenotipes_GLS.csv'
#     for i in range(1560, 1900, 50):
#         printLog('PhenoClassifier: GLS loading lines {i}'.format(i=i))
#         embs = embdFunc.readEmbsFromCSV(fileName, nEmbs=50, skip=i)
#         del embs
#     embs = embdFunc.refreshEmbsMulti(embs)

if __name__ == '__main__':
#     multiClassGLS()
    multiClassMS()
#     X = [[0., 0.], [1., 1.]]
#     y = [0, 1]
#     clf = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(5, 2), random_state=1)
#     clf.fit(X, y)