"""
Created on Sep 4, 2015

@author: Renat
"""

import os, sys, glob
import numpy as np
from myFunc import sort_nicely
import Embryos
from varLookup import printLog, nCores, FOLDER_IN
import multiprocessing as mp
import timeAlign
import traceback
import git
from functools import wraps
from db_utils_embryos import initialize_db, get_row


def emb_try_decorator(func):
    @wraps(func)
    def try_func(*arg, **kwargs):
        try:
            return func(*arg, **kwargs)
        except Exception as e:
            print(traceback.print_exc())
            if isinstance(arg, Embryos.Embryo):
                printLog('ERROR in {s}/{r}/{d}/{l} '.format(s=arg.strain, r=arg.RNAi, d=arg.date, l=arg.label))
                return arg
            elif isinstance(arg, str):
                printLog('ERROR in {} '.format(arg))
            return None
    return try_func


def parallelizer(func, my_list, cores=6):
    pool = mp.Pool(processes=cores)
    result = pool.map(func, my_list)
    pool.close()
    pool.join()
    return result


@emb_try_decorator
def loadGLSEmbryo(folder, verb=True, check_version=True):
    emb = Embryos.GSEmbryo(folder, verb, check_version)
    return emb


@emb_try_decorator
def loadMSEmbryo(folder, verb=True, check_version=True):
    emb = Embryos.MSEmbryo(folder, verb, check_version)
    return emb


def load_emb_by_id(myid, verb=True, check_version=True):
    conn, curs = initialize_db()
    row = get_row(myid, curs, 'embryos', ['strain', 'file_path'])
    if row['strain'] == 'GLS':
        return loadGLSEmbryo(FOLDER_IN + row['file_path'][3:], verb, check_version)
    else:
        return loadMSEmbryo(FOLDER_IN + row['file_path'][3:], verb, check_version)


def loadEmbs(folder, inds=None):
    """
    returns a list of Embryo objects loaded from the specified folder
    folder: folder with embryo folders inside (embryoFolderFlag = False) or a single embryo folder (embryoFolderFlag=True)
    inds: indices of embryos to load
    MS: whether to load MSEmbryo or GSEmbryos
    """
    if 'MS' in folder.split('/'):
        MS = True
    else:
        MS = False
    folders = []  # list of embryo folders
    if inds is not None: inds = np.array(inds)
    if folder.split('/')[-2][0] == 'E' or folder.split('/')[-2][
        0] == 'l':  # for the condition when running only a single embryo... adds only one item to list
        imNames = glob.glob('{0}/*.tif'.format(folder))
        if len(imNames) > 0: folders = [folder]
    else:  # for the condition where running multiple embryos... adds all embryo folders to list
        for d in os.listdir(folder):
            if os.path.isdir(os.path.join(folder, d)) and d[0] == 'E':
                imNames = glob.glob('{0}/*.tif'.format(os.path.join(folder, d)))
                if len(imNames) > 0:
                    folders.append(os.path.join(folder, d))
                else:
                    print('no tif files in {0}'.format(os.path.join(folder, d)))
        sort_nicely(folders)
    folders = np.array(folders)
    if inds is None:
        inds = np.arange(len(folders))  # if no index is specified, use all embryos
    else:
        inds = inds[np.where(
            inds < len(folders))]  # checks to make sure your indices are always below the number of available embryos
    folders = folders[inds]
    if folders.size < 3 and MS:
        embs = [loadMSEmbryo(f) for f in folders]
    elif folders.size < 3:
        embs = [loadGLSEmbryo(f) for f in folders]
    else:
        pool = mp.Pool(processes=nCores)
        if MS:
            embs = pool.map(loadMSEmbryo, folders)  # a list of Embryo objects (either all GS or all MS)
        else:
            embs = pool.map(loadGLSEmbryo, folders)
        pool.close()
        pool.join()
    embs = [emb for emb in embs if emb is not None]  # removes none type objects (i.e. may be present if wrong number of
    # files present and exception is raised)
    return embs  # a list of Embryo objects


@emb_try_decorator
def refreshEmb(emb):
    if emb.loaded:
        emb.refresh()
    return emb

@emb_try_decorator
def updateEmb(emb):
    if emb.loaded:
        emb.updateParams()
    return emb


@emb_try_decorator
def setEmbScales(emb):
    emb.setScaleCutoff()
    return emb

@emb_try_decorator
def timeAlignMSSingle(emb):
    timeAlign.alignTimeMS(emb)
    return emb


@emb_try_decorator
def timeAlignGLSSingle(emb):
    timeAlign.alignTimeGLS(emb)
    return emb


def refreshEmbsMulti(embs, cores=nCores):
    return parallelizer(refreshEmb, embs, cores)


def updateEmbsMulti(embs, cores=nCores):
    return parallelizer(updateEmb, embs, cores)


def setScalesEmbsMulti(embs, cores=nCores):
    return parallelizer(setEmbScales, embs, nCores)


def timeAlignMultiMS(embs, cores=nCores):
    return parallelizer(timeAlignMSSingle, embs, nCores)


def timeAlignMultiGLS(embs, cores=nCores):
    return parallelizer(timeAlignGLSSingle, embs, nCores)


def get_embryo_revision_version():
    try:
        path = "../embryos_analyze/"
        repo = git.Repo(path=path)
    except:
        path = "../../embryos_analyze/"
        repo = git.Repo(path=path)
    version = list(repo.iter_commits())[0]
    # output = subprocess.Popen(["git", "describe", "--all", "--long"], stdout=subprocess.PIPE).communicate()[0]
    # version = output.split('-')[-1]
    return str(version)[:7]


if __name__ == '__main__':
    emb = load_emb_by_id(7356)
    emb.showSigmFit()
    fig = emb.show('MoI0R')
    fig.show()
