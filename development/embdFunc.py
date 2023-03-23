"""
Created on Sep 4, 2015

@author: Admin
"""
from myFunc import readCSV
from myMath import getSig
from varLookup import PARAM_NAMES
from fileLookup import FILENAME_GENE_NAMES, FILENAME_GENE_WB_INFO
import numpy as np
import git
import pandas as pd
global embsT
embsT = []


def getNDDistance(v1, v2, pNamesUse):
    ''' calculates distance between two vectors in N dimensional space. The vector component that are nan treated as dimensionality reduced '''

    if isinstance(v2, str) and v2 == 'control':
        v2 = np.ones_like(v1) * np.nan  # if v2 is control, set it to the origin
        for i in range(len(PARAM_NAMES)):  # iterates through the index of parameter names to analyze
            if PARAM_NAMES[i] in pNamesUse[0]:
                v2[i] = 0.
            if PARAM_NAMES[i] in pNamesUse[1]:
                v2[len(PARAM_NAMES) + i] = 0.
    if v1.size != v2.size:
        print('embdFunc.getNDDistance: wrong dimensionality of the supplied vectors')
        raise
    N = len(pNamesUse[0]) + len(pNamesUse[1])
    vOverlap = np.invert(np.isnan(
        v1 * v2))  # overlapping vector, checks whether a parameter (coordinate) is not nan for each considered group (RNAi)
    if sum(vOverlap) == 0:
        return np.inf  # if there are no overlapping coordinates, returns infinity
    i1 = np.sum(np.invert(np.isnan(v1))) - np.sum(
        vOverlap)  # sums all coordinates that are not equal to nan for v1 and subtracts overlap
    i2 = np.sum(np.invert(np.isnan(v2))) - np.sum(
        vOverlap)  # sums all coordinates that are not equal to nan for v2 and subtracts overlap
    dist = np.sqrt(1. * N / (N - i1 - i2)) * np.linalg.norm((v1 - v2)[
                                                                vOverlap])  # imposes weighting to penalize for non-overlapping parameters and caluculates Euclidean Distance between 2 points
    return dist


def getNDAvgVector(vectors):
    ''' calculates average between list of vectors in N dimensional space. The vector component that are zero treated as dimensionality reduced '''
    vectors = np.array(vectors)  # gets list of numbers (vectors) from RNAi class and converts it to an array
    avg = np.nanmean(vectors, axis=0)  # average of everything that is not nan
    avg[np.where(np.isnan(avg))] = np.nan  # fix...
    return avg

def merge_lookup_files():
    gene_key= pd.read_csv(FILENAME_GENE_NAMES)
    wb_info= pd.read_csv(FILENAME_GENE_WB_INFO)
    merge_file= pd.merge(gene_key,wb_info,how='left')
    print(merge_file.head(5))
    merge_file.to_csv('Z:/merge_file.csv')

def get_interacting_genes(embdN, printflag=True):
    '''
    Retrieves a list of genes and evidence for interaction from CSV file- this information was mined from SimpleMine on wormbase 6/28/18
    :param embdN:
    :return: list of genes, i.e. (
    '''
    gene = getGeneName(embdN)
    data = pd.read_csv(FILENAME_GENE_WB_INFO)
    # print(data.head(5))
    gene_data = data[data.Sequence_Name == gene]
    interacting_genes = gene_data['Interacting_Gene'].tolist()
    ig_split = interacting_genes[0].split(',')
    if printflag: print('EMBD{0}: {1}, interacts with {2}'.format(embdN, getGeneCommonName(embdN), ig_split))
    return ig_split

def check_for_interaction(test_string, embdN):
    '''
    This function checks for an interaction with a particular gene or can be used to list all genes that have a specific
    type of interaction i.e. ('Physical', 'Regulatory', 'Genetic')
    :param test_string: can be a gene name 'mbk-1' or type of interaction 'Regulatory' (note that all interaction types
     are capatilized
    :param embdN: EMBD unique id number
    :return:
    '''
    test = test_string  # this can be a gene name (i.e. 'mbk-1') or can be a type of interaction ('Physical', 'Regulatory', 'Genetic')
    interacting = get_interacting_genes(embdN, printflag=False)
    print('all interactions for {0}:\n {1}'.format(embdN, interacting))
    print('interactions matching string:')
    print [s for s in interacting if test in s]

def getGeneName(embdN):
    if isinstance(embdN, str): embdN = int(embdN[-4:])
    if embdN == 0: return 'control'
    data = readCSV(FILENAME_GENE_NAMES, ',')[1:]
    dataA = np.array(data)
    ind = np.where(dataA == 'EMB_P{n:04}'.format(n=embdN))[0]
    return data[ind[0]][1]


def getGeneCommonName(embdN):
    if isinstance(embdN, str): embdN = int(embdN[-4:])
    if embdN == 0: return 'control'
    data = readCSV(FILENAME_GENE_NAMES, ',')[1:]
    dataA = np.array(data)
    ind = np.where(dataA == 'EMB_P{n:04}'.format(n=embdN))[0]
    if ind.size > 0 and bool(data[ind[0]][2]):  # checks if anything is present in the cell
        return data[ind[0]][2]
    else:
        return 'unnamed'


def getEMBDNumber(geneName):
    if geneName == 'control': return 0
    data = readCSV(FILENAME_GENE_NAMES, ',')[1:]
    dataA = np.array(data)
    ind = np.where(dataA == geneName)[0]
    if ind.size > 1:
        #         print('gene duplicate: {g} {r1} & {r2}'.format(g=geneName, r1=data[ind[0]][0], r2=data[ind[1]][0]))
        ind = ind[0]
    elif ind != [None]:
        #         print int(data[ind][0][-4:])
        return int(dataA[ind][0][0][-4:])
    else:
        #         print 'gene not in EMBD list'
        return None


def getEMBDNumberList(geneNameList):
    geneNameList = np.array(geneNameList)
    data = readCSV(FILENAME_GENE_NAMES, ',')[1:1345]
    data = np.array(data)
    #     a = np.stack((data[:,1],np.arange(data[:,1].size)), axis=-1)
    index = np.lexsort((np.arange(data[:, 1].size), data[:, 1]))
    #     index = np.argsort(a, order=['name', 'EMBD_ID'])
    sorted_data = data[index]
    sorted_index = np.searchsorted(sorted_data[:, 1], geneNameList)
    gene_index = np.take(index, sorted_index, mode='clip')
    mask = data[gene_index, 1] != geneNameList
    resInd = np.ma.array(gene_index, mask=mask)
    result = []
    for i in range(mask.size):
        if not mask[i]:
            result.append(int(data[resInd[i]][0][-4:]))
        else:
            result.append(np.nan)
    return np.array(result)


# def readEmbsFromCSV(fileName, nEmbs=None, skip=0):
#     data = readCSV(fileName, ',')[1 + skip:]  # skip first line
#     embs = []
#     if nEmbs is None or nEmbs > len(data): nEmbs = len(data)
#     foldersMS, foldersGLS = [], []
#     for d in data[:nEmbs]:
#         if d[2] == 'EMBD0000':
#             folder = FOLDER_IN + 'cropped/{rna}/{strain}/{date}/{emb}/'.format(strain=d[0], date=d[1], rna=d[2],
#                                                                                emb=d[3])
#         else:
#             folder = FOLDER_IN + 'cropped/{rna}/{strain}/{emb}/'.format(strain=d[0], date=d[1], rna=d[2], emb=d[3])
#         if not os.path.exists(folder):
#             print('embryo from csv file does not exist {f}'.format(f=folder))
#             raise
#         if d[0] == 'MS':
#             foldersMS.append(folder)
#         else:
#             foldersGLS.append(folder)
#     pool = mp.Pool(processes=nCores)
#     if len(foldersMS) > 0:
#         embsMS = pool.map(loadMSEmbryo, foldersMS)
#     else:
#         embsMS = []
#     if len(foldersGLS) > 0:
#         embsGLS = pool.map(loadGLSEmbryo, foldersGLS)
#     else:
#         embsGLS = []
#     pool.close()
#     pool.join()
#     return embsMS + embsGLS


def get_param_names():
    """
    generates names for all parameters from varlookup.PARAM_NAMES
    :return: list of names
    """

    p_names = []
    for name in PARAM_NAMES:
        p_names.append(name + '_GLS')
        p_names.append(name + '_MS')
    return p_names


def split_params(params_all):
    """
    splits all params into lists of PARAM_NAMES GLS and PARAM_NAMES MS and the rest
    :param dict: dictionary of the parameters
    :return: three dictionaries: PARAM_NAMES GLS and PARAM_NAMES MS and the rest
    """
    params_gls = dict(zip(PARAM_NAMES, [np.nan for p in PARAM_NAMES]))
    params_ms = dict(zip(PARAM_NAMES, [np.nan for p in PARAM_NAMES]))
    other = {}
    for key in params_all:
        if key[-3:] == 'GLS' and key[:-4] in PARAM_NAMES:
            if params_all[key] is None:
                params_gls[key[:-4]] = np.nan
            else:
                params_gls[key[:-4]] = params_all[key]
        elif key[-3:] == '_MS' and key[:-3] in PARAM_NAMES:
            if params_all[key] is None:
                params_ms[key[:-3]] = np.nan
            else:
                params_ms[key[:-3]] = params_all[key]
        else:
            other[key] = params_all[key]
    return params_gls, params_ms, other


def get_rnai_revision_version():
    try:
        repo = git.Repo()
    except:
        path = "../"
        repo = git.Repo(path=path)
    version = list(repo.iter_commits())[0]
    # output = subprocess.Popen(["git", "describe", "--all", "--long"], stdout=subprocess.PIPE).communicate()[0]
    # version = output.split('-')[-1]
    return str(version)[:7]


def make_dict_float(mydict):
    for key in mydict:
        if not isinstance(mydict[key], int):
            mydict[key] = float(np.round(mydict[key], getSig(mydict[key])+5))
    return mydict


if __name__ == '__main__':
    interacting_genes = get_interacting_genes(16)
    # check_for_interaction('mbk-2', 63)




    pass