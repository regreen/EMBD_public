'''
Created on Feb 28, 2017

@author: Becky, Renat
'''
from myFunc import readCSV, saveImagesMulti
import numpy as np
from embdFunc import getGeneName, getEMBDNumber, getGeneCommonName
from fileLookup import FILENAME_PAD, FOLDER_TOP_HITS_MOVIES
from sklearn import metrics
from myFigure import myFigure
from varLookup import DT0_GLS_2_MS
from RNAiClass import RNAiClass
import networkx as nx
from networkx.algorithms.approximation.clique import max_clique, clique_removal
from networkx.algorithms.approximation import node_connectivity
import networkx.algorithms.clique as clique
from networkx.algorithms.bipartite.centrality import betweenness_centrality
from myMath import getSubIndex
from dataValidation import getBinaryManualData
from numpy import argsort
from RNAiClass import RNAiClass
from Embryos import MSEmbryo, GSEmbryo
import cv2
import numpy as np
from PIL import ImageFont, ImageDraw, Image
from copy import copy
from db_utils_embryos import initialize_db, get_columns
from avgCalculator import getGLSEmbs, getMSEmbs



def PADFilter(rnaName, data, nTop):
    """
    Returns top ranked hits by pad
    :param rnaName: seed (EMBD string)
    :param data: data file
    :return: a list of lists for nTop hits. For a given seed gene (rnaName) each list consists of top hit embd name,
    PAD value and the number of overlapping dimensions('EMBD0076', 0.934496703417, 23.0)
    """
    if isinstance(rnaName, int):
        rnaName = 'EMBD{r:04}'.format(r=rnaName)
    nTop = nTop
    dataA = np.array(data)
    inds = np.where(dataA[:, 0] == rnaName)[0]
    dRNA = [data[i] for i in inds]
    pads = [p[6] for p in dRNA]  # this is PAD
    inds = np.argsort(pads)[::-1]
    res = [(dRNA[i][1], dRNA[i][6], dRNA[i][4]) for i in inds[1:nTop + 1]]
    return res


def CSIFilter(rnaName, data):
    dataA = np.array(data)
    inds = np.where(dataA[:, 0] == rnaName)[0]
    dRNA = [data[i] for i in inds]
    pads = [p[5] for p in dRNA]  # this is CSI
    inds = np.argsort(pads)[::-1]
    res = [(dRNA[i][1], dRNA[i][5], dRNA[i][4]) for i in inds[1:11]]
    return res


def distFilter(rnaName, data):
    dataA = np.array(data)
    inds = np.where(dataA[:, 0] == rnaName)[0]
    dRNA = [data[i] for i in inds]
    pads = [p[2] for p in dRNA]  # this is Euclidean distance
    inds = np.argsort(pads)[::-1]
    res = [(dRNA[i][1], dRNA[i][2], dRNA[i][4]) for i in inds[1:11]]
    return res


def getPADFromFile(rn1, rn2, data):
    ind = np.where(data == 'EMBD{n:04}'.format(n=rn1))[0]
    if ind.size > 0:
        ind2 = np.where(data[ind] == 'EMBD{n:04}'.format(n=rn2))[0]
    else:
        print('!!!!!!!!!!!!!!!!!!!ALARM', rn1)
    if ind2.size > 1: ind2 = ind2[0]
    val = float(data[ind][ind2][6])
    return val


def getCSIFromFile(rn1, rn2, data):
    ind = np.where(data == 'EMBD{n:04}'.format(n=rn1))[0]
    ind2 = np.where(data[ind] == 'EMBD{n:04}'.format(n=rn2))[0]
    if ind2.size > 1: ind2 = ind2[0]
    val = float(data[ind][ind2][5])
    return val


#     for d in data:
#         if d[0]=='EMBD{n:04}'.format(n=rn1) and d[1]=='EMBD{n:04}'.format(n=rn2): return d[6]
#         elif d[0]=='EMBD{n:04}'.format(n=rn2) and d[1]=='EMBD{n:04}'.format(n=rn1): return d[6]
#     return np.nan

def getInteractions(geneList):
    fileName = 'Z:/IntegratedInteractome.csv'
    data = readCSV(fileName, ',')[1:]
    res = []
    for d in data:
        if d[0] in geneList and d[1] in geneList and d[0] != d[1]:
            res.append((getEMBDNumber(d[0]), getEMBDNumber(d[1])))
    return np.array(res)


def makeNetwork(edgeData):
    '''loads genes as nodes and PAD values as edges using networkx'''
    G = nx.Graph()  # creates empty graph with no nodes or edges, can assign graph attributes i.e. G = nx.Graph(day="Friday")
    for i in range(1, 504):
        G.add_node(i)
    #     G.add_nodes_from([range(1,504)], color = 'blue') #add a list of nodes
    th = 0.74
    #     th2= 0.7
    edgeData = np.array(edgeData)
    inds = np.where(edgeData[:, 6].astype(float) > th)
    #     inds = getSubIndex(edgeData[:,6].astype(float), th2, th)
    edgeData = edgeData[inds]
    for row in edgeData:
        n1 = int(row[0][-4:])  # node1
        n2 = int(row[1][-4:])  # node2
        weight = float(row[6])  # PAD
        #         weight = float(row[5])#CSI
        if n1 != n2:
            G.add_edge(n1, n2, weight=weight)  # adds edges
        #     for i in range(1,3):
        #         print G[i]
    return G  # graph of all genes with PAD as edge value


def findCliques(G, popList):
    '''takes in graph(network) and identifies groups of genes that are interconnected'''
    #     print clique.graph_number_of_cliques(G)#returns the number of maximal cliques in G
    #     betweenness_centrality(G, nodes)
    #     print clique.cliques_containing_node(G,9)#returns a list of cliques containing the given node
    wtList, lowPenList, highPenList, arrestList, cfsList, featureFullList = getBinaryManualData()
    #     nbunch = [181]   #21, 181, 31, 34, 63, 55, 9, 426, 10]
    #     nbunch = range(1,504)

    nbunch = []
    for g in highPenList:
        if g not in wtList: nbunch.append(g)
    #     nbunch = highPenList

    goodCliques = []
    cliqueGenes = []
    # wtGenes = [3, 23, 41, 62, 66, 84, 85, 91, 92, 96, 127, 132, 148, 149, 174, 176, 179, 202, 206, 207, 213, 227, 228,
    #            236, 245, 249, 251, 252, 254, 257, 267, 271, 272, 275, 279, 280, 292, 294, 296, 299, 307, 310, 311, 314,
    #            319, 329, 331, 332, 347, 352, 367, 368, 373, 375, 378, 381, 383, 384, 391, 392, 394, 399, 407, 409, 414,
    #            416, 418, 421, 429, 437, 442, 454, 456, 458, 459, 465, 467, 468, 475, 476, 481, 482, 483, 484, 486, 487,
    #            493, 494, 495, 497, 498, 499]
    wtGenes = []
    for n in nbunch:
        if n not in popList and n not in wtGenes:
            #         print "neighbors of {0}".format(n), G.neighbors(n)
            dict_of_cliques = nx.cliques_containing_node(G, nodes=n, cliques=None)
            if len(dict_of_cliques[0]) > 1:
                if len(dict_of_cliques[0]) > 2 or (
                        dict_of_cliques[0][0] not in popList and dict_of_cliques[0][1] not in popList):
                    for clique in dict_of_cliques:
                        if clique not in goodCliques and np.intersect1d(clique, wtGenes).size < 1:  #
                            goodCliques.append(clique)
                            cliqueGenes += clique
                            print "clique of {0}".format(n), clique
    cliqueGenes = np.unique(cliqueGenes).tolist()
    #     print cliqueGenes
    print goodCliques


#     btwnCent= betweenness_centrality(G, nbunch)#dictionary of each node and the betweenness centrality of that node (measure of the importance of the node in the network)
# #     print btwnCent
# #     print sorted(btwnCent.values(), reverse=True)
#     for n in nbunch:
#         for key, value in btwnCent.items():
#             if key == n: 
#                 print key, value, G.neighbors(key)
# #         if value>0.02and value<0.14: 
# #             print key, value, G.neighbors(key)
#                 dict_of_cliques=nx.cliques_containing_node(G,nodes=key,cliques=None)
#                 print "clique of {0}".format(n), dict_of_cliques

#     print "max", max_clique(G)


def getPairsParams(seed, data):
    '''get PAD values and most correlated parameter names between seed and top ten neighbors (i.e. 156)'''

    rnaSeed = RNAiClass(seed)
    rnaSeed.setNDPosition()
    pairList = PADFilter('EMBD{r:04}'.format(r=seed),
                         data, 10)  # returns a list of top hits by PAD values- each as gene, PAD, dimensions
    rList = [RNAiClass(int(p[0][5:])) for p in pairList]
    print('Pairs for EMBD{r:04}'.format(r=seed))
    i = 0
    for p in pairList:
        print('{r}: PAD={pad}, nOverlap={o}, {par}'.format(r=p[0], pad=rnaSeed.getPAD(rList[i]), o=p[2],
                                                           par=rnaSeed.getSimParams(rList[i])))
        i += 1


def getCommonPairs(seedList, data):
    allPairList = []  # list of EMBD numbers that are paired with seed
    allPADvals = []  # list of PAD values that from genes paired with seed
    for seed in seedList:
        topGenes = PADFilter('EMBD{r:04}'.format(r=seed),
                             data, 10)  # returns list of lists- each list has gene, PAD and number of dimensions
        topGenes = np.array(topGenes)
        allPairList.append(topGenes[:, 0])  # appends genes
        allPADvals.append(topGenes[:, 1])  # appends PAD vals
    uniqueList, countList = np.unique(allPairList,
                                      return_counts=True)  # returns two lists- first list is all unique genes, second list is count for how many times genes are in list
    allPairList = np.array(allPairList)
    allPADvals = np.array(allPADvals)
    uniquePADvals = []  # list of PAD values for list of genes (reduced to unique list)
    for gene in uniqueList:
        uniquePADvals.append(np.mean(allPADvals[np.where(allPairList == gene)].astype(
            float)))  # finds the TRUE indicies where gene is in allPairList, retrieves the PAD values from those indices. If there are more than one pad value, appends the average
    index = np.lexsort((uniquePADvals,
                        countList))  # returns a list of indices corresponded to sorting by countList, then by uniquePADvals
    for i in index[::-1]:  # [::-1] slices through the list in reverse order
        print uniqueList[i], countList[i], uniquePADvals[i]
    return [int(uniqueList[i][-4:]) for i in index[::-1]]


#     np.argsort(uniqueList[:,1])

def get_best_match(RNAi, data):
    '''input- RNAi object- i.e EMBD0076
    Finds top matches using PAD (orCSI) by calling a CSV file, returns a list of ints [12,33,126] representing the genes with the closest phenotypes'''

    seed = RNAi[-4:]
    seedObj = RNAiClass(int(seed))
    top_hits = [
        seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
    top_hit_data = np.array(PADFilter(RNAi, data, 5))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
    top_hit_vals = np.array(PADFilter(RNAi, data, 5))[:, 1]  # retrieves the PAD values for top matches

    ## top_hit_data = np.array(CSIFilter(RNAi, data))[:,0] #retrieves the top matches i.e. 'EMBD0076' based on CSI
    ## top_hit_vals = np.array(CSIFilter(RNAi, data))[:,1] #retrieves the CSI values for top matches

    for top_hit in top_hit_data:
        top_hits.append(top_hit[-4:])  # slices the last four digits of string to remove the 'EMBD'
    print top_hits, seed, top_hit_vals
    return top_hits, seed, top_hit_vals


def show_best_match_movies(top_hits, seed, top_hit_vals, saveFlag=True, recip_saveFlag=False):
    '''input- get_best_match returns top_hits, seed, top_hit_vals. Top hits is a list of RNAi conditions [76, 23,12,112]
    and seed RNAi condition i.e. 187. Finds top matches using PAD (orCSI). Can also use make_reciprocal_match_movies
     function to populate arguments. This function assembles a time aligned, labeled movie of the
    representative embryo for each condition- assembles GLS and MS side-by-side
    :param top_hits: top hits is a list of strings ['0010', '0243', '0090', '0321', '501']
    :param seed: seed is an EMBD number in string format '0010'
    :param top_hit_vals: is a list of PAD (could be CSI) values in string format ['1.0', '0.80', '0.79', '0.7621', '0.7501']
    :param saveFlag: if True it saves the .tif file output for top_hits by EMBD number to the database folder
    :param recip_saveFlag: if True it saves the .tif file output for reciprocal top_hits, by EMBD number to the database folder
    '''

    from fileLookup import FILENAME_WEB_CONTENT
    import os

    seed_obj = RNAiClass(int(seed))
    rep_emb_all = []  # a list of two lists: rep embs for GLS and for MS
    rep_emb_gall = []  # contains the representative embryo object for the each of the top hits GLS
    rep_emb_mall = []  # contains the representative embryo object for the each of the top hits MS
    # repEmbNamesGLS = []
    # repEmbNamesMS= []

    for RNA in top_hits:   
        rnai = RNAiClass(int(RNA))
        rnai.getEmbryos()
        rnai.asignUse()
        # allE = [rnai.GLSUseEmbs, rnai.MSUseEmbs]
        rep_emb_g, rep_emb_m = rnai.getRepEmb(
            seed_obj)  # a list of representative embryo object [GLS,MS] **note that this gets the top hit emb that looks most like the seed gene, not the 'average' emb for this condition
        # repEmbGLS = repEmbG.folder.split('/')[2] + '_' + repEmbG.folder.split('/')[4]
        # repEmbMS = repEmbM.folder.split('/')[2] + '_' + repEmbM.folder.split('/')[4]
        rep_emb_gall.append(rep_emb_g)
        rep_emb_mall.append(rep_emb_m)
        # repEmbNamesGLS.append(repEmbGLS)
        # repEmbNamesMS.append(repEmbMS)

    rep_emb_all.append(rep_emb_gall)
    rep_emb_all.append(rep_emb_mall)
    for emblist in rep_emb_all:  # iterates through GLS embs, then through MS embs
        maxProjEmbs = []
        shape = []
        name = []
#         for emb in emblist:
        for i in range(len((emblist))):
            emb = emblist[i]
            if i >0:
                pad_to_seed = top_hit_vals[i-1][:5]
            else:
                pad_to_seed = str(1.000)
            embProjTime = []  # projections for each timepoint for a given embryo
            if emb.folder.split('/')[3] == 'GLS':
                name.append(rep_emb_g.folder.split('/')[2] + '_' + rep_emb_g.folder.split('/')[4])
            elif emb.folder.split('/')[3] == 'MS':
                name.append(rep_emb_m.folder.split('/')[2] + '_' + rep_emb_m.folder.split('/')[4])
            for t in range(-2,
                           27):  # goes through each timepoint, makes max proj, draws RNAi/embryo name,-- here time has been converted  (nearest slide used for dimensionless time)- displays time points starting 3 back from t0.
                if emb.strain == 'MS':
                    embProj, (xG, yG), (xR, yR) = emb.getMaxProj(
                        t - DT0_GLS_2_MS)  # gets maximum projection for the slide (t). To display similar stages, we shift t by delta t0 value for MS
                else:
                    embProj, (xG, yG), (xR, yR) = emb.getMaxProj(t)  # gets maximum projection for the slide (t)
                embProj = (255 * embProj).astype(np.uint8)  # converts to 8bit

                embProjIm = cv2.cvtColor(embProj, cv2.COLOR_BGR2RGB)  # Convert the image to RGB (OpenCV uses BGR)
                pil_im3 = Image.fromarray(embProjIm)  # Pass the image to PIL
                draw = ImageDraw.Draw(pil_im3)  #
                fontPil3 = ImageFont.truetype("arial.ttf", size=18, index=0)  # use a truetype font
                # fontPil4 = ImageFont.truetype("ariali.ttf", size=18, index=0)  # use a truetype font
                # embName = emb.folder.split('/')[2] + ' ' + emb.folder.split('/')[4]
                # embName = emb.folder.split('/')[2] + ' ' + getGeneCommonName((emb.folder.split('/')[2])[-4:])
                # if emb.strain == 'GLS':
                #     draw.text((10, 10), embName, font=fontPil3,
                #               fill=(255, 255, 255))  # Draw the text (embName) on the image

                embName = emb.folder.split('/')[2] + ' ' + getGeneCommonName((emb.folder.split('/')[2])[-4:])

                embName2 ='PAD = ' + pad_to_seed
                if emb.strain == 'GLS':
                    draw.text((5, 5), embName, font=fontPil3,
                              fill=(255, 255, 255))  # Draw the text (embName) on the image
                elif emb.strain == 'MS':
                    draw.text((5, 5), embName2, font=fontPil3,
                              fill=(255, 255, 255))  # Draw the text (embName) on the image
                embProj = cv2.cvtColor(np.array(pil_im3), cv2.COLOR_RGB2BGR)
                embProjTime.append(embProj)  # adds to list for each emb
            maxProjEmbs.append(
                embProjTime)  # appends each list for each embryo to maxProjEmbs for each RNAi for each strain
            embProjTime = np.array(embProjTime)  # converts to array (this is mostly for next line)
            shape.append(
                embProjTime.shape)  # appends shape of each emb (this is needed for layout of all embs in final output)
            emb.image = None  # deletes images to save on memory
        shape = np.array(shape)
        time = max(shape[:, 0])  # finds number of timepoints for the embryo with the longest timeseries
        width = max(shape[:, 2])  # finds width of the widest image
#         height = sum(shape[:, 1]) + (len(maxProjEmbs) - 1)  # finds the sum of all heights for each RNAi
#         max_height = max(shape[:, 1])
        max_height = 150
        print max_height
        height = max_height * len(maxProjEmbs)# finds the tallest emb for each RNAi and multiplies it by the number of embryos
#         height = max(shape[:, 1])* len(maxProjEmbs)# finds the tallest emb for each RNAi and multiplies it by the number of embryos
        offset = 0  # len(maxProjEmbs)+2 #not sure what is reasonable offset
        image = np.zeros(
            (time, height + offset, int(width + 0.5 * offset), 3))  # creates a black background to place max projects

        yTot = 0  # used to define spacing as each embryo is added
        for i in range(len(maxProjEmbs)):  # compiles all timepoints and stacks movies for each strain
            maxP = maxProjEmbs[i]
            t = shape[i][0]
            x = shape[i][2]  # defines length of projection
            y = max_height # defines height of background to place projection
            yI = shape[i][1]  # defines height of projection of image
            #                 eName= name[i]#i.e. emb1
            image[:, yTot:yTot + yI, int(0.25 * offset):x,
            :] = maxP  # redefines pixel values for specified coordinates [yRange,xRange] on black background to maxP image values
#             image[:, yTot:yTot + y, int(0.25 * offset):x,
#             :] = maxP  # redefines pixel values for specified coordinates [yRange,xRange] on black background to maxP image values
            yTot += y  # adds y to yTot so that next image is properly positioned below previous image

        image = np.array(image).astype(np.uint8)  # required because black background is not uint8
        for tp in image[:1]:
            fig = myFigure()
            fig.imshow(tp, colorbar=False)  # plots image onto figure
            fig.noAxis()
        # fig.show()
        if emb.folder.split('/')[3] == 'GLS':
            imageGLS = image
        else:
            imageMS = image

    widthComp = max(imageGLS.shape[2], imageMS.shape[2])  # finds the width of the widest set of images (GLS or MS)
    heightComp = max(imageGLS.shape[1],
                     imageMS.shape[1])  # finds the length of the widest set of images (GLS or MS)

    hOffset = 50
    imageComp = np.zeros((time, heightComp + hOffset, int(2 * widthComp + 0.5 * offset),
                          3))  # creates a black background to place max projects
    imageComp = np.array(imageComp).astype(np.uint8)  # makes black background 8 bit

    text_to_show = 'EMBD' + seed
    if getGeneName(seed)== getGeneCommonName(seed):  # checks to make sure we dont print CE ID twice when a gene is unnamed
        text_to_show2 = '{0}         '.format(getGeneName(seed))
    else:
        text_to_show2 = '{0}     {1}'.format(getGeneName(seed), getGeneCommonName(seed))
    cv2_im_rgb = cv2.cvtColor(imageComp[0], cv2.COLOR_BGR2RGB)  # Convert the image to RGB (OpenCV uses BGR)
    #             pil_im= Image.frombytes("L",imageComp.size, imageComp.tobytes())
    pil_im = Image.fromarray(cv2_im_rgb)  # Pass the image to PIL
    draw = ImageDraw.Draw(pil_im)  #
    fontPil = ImageFont.truetype("arial.ttf", size=22, index=0)  # use a truetype font
    fontPil2 = ImageFont.truetype("ariali.ttf", size=22,
                                  index=0)  # use a truetype font (arial italic is 'ariali' and arial bold is 'arialbd')
    draw.text((15, 10), text_to_show, font=fontPil, fill=(255, 255, 255))  # Draw the text
    draw.text((widthComp - 60, 10), text_to_show2, font=fontPil2, fill=(255, 255, 255))  # Draw the text
    #         pil_im.show()
    #         cv2_im_processed = cv2.cvtColor(np.array(pil_im), cv2.COLOR_RGB2BGR)#Get back the image to OpenCV

    cv2_im_processedAll = []
    for tmpt in range(-2, 27):
        im_pil_t = copy(pil_im)  # makes a copy of imageComp, so old timestamp is not there for the next timept
        draw = ImageDraw.Draw(im_pil_t)
        #             draw.text((2*widthComp-60,10), str((tmpt-1)*20)+'m', font=fontPil, fill=(0,0,0)) # Note- time displayed is GLS normalized time
        draw.text((2 * widthComp - 60, 10), str((tmpt) * 20) + 'm', font=fontPil,
                  fill=(255, 255, 255))  # Note- time displayed is GLS normalized time
        cv2_im_processed = cv2.cvtColor(np.array(im_pil_t), cv2.COLOR_RGB2BGR)  # Get back the image to OpenCV
        cv2_im_processedAll.append(cv2_im_processed)

    cv2_im_processedAll = np.array(cv2_im_processedAll)
    imageComp[:, :cv2_im_processedAll.shape[1], :cv2_im_processedAll.shape[2], :] = cv2_im_processedAll
    imageComp[:, hOffset:imageGLS.shape[1] + hOffset, :imageGLS.shape[2],
    :] = imageGLS  # re-assigns pixel values for GLS composite on black background to left edge
    imageComp[:, hOffset:imageMS.shape[1] + hOffset, widthComp:widthComp + imageMS.shape[2],
    :] = imageMS  # re-assigns pixel values for MS images to the right of GLS images on black background
    for tp in imageComp[:1]:
        fig = myFigure()
        fig.imshow(tp, colorbar=False)  # plots image onto figure
        fig.noAxis()
    # fig.show()
    # saveImagesMulti(imageComp, FOLDER_TOP_HITS_MOVIES + '{0}_tophits.tif'.format('EMBD' + seed))

    # new edits
    if saveFlag:
        folder = FILENAME_WEB_CONTENT + '{0}/'.format('EMBD' + seed)
        if not os.path.exists(folder):
            os.makedirs(folder)
        saveImagesMulti(imageComp, folder + '{0}_tophits.tif'.format('EMBD' + seed))
    if recip_saveFlag:
        folder = FILENAME_WEB_CONTENT + '{0}/'.format('EMBD' + seed)
        if not os.path.exists(folder):
            os.makedirs(folder)
        saveImagesMulti(imageComp, folder + '{0}_reciprocal_hits.tif'.format('EMBD' + seed))

def make_movies(data):
    """
    make movie of time aligned top hits
    :param data: rnai pairs pads read from a csv file
    :return:
    """

    RNAi_list = [i for i in range(445,504)]
    for RNAi in RNAi_list:
        RNAi = 'EMBD{r:04}'.format(r=RNAi)
        print("STARTING {0}-----------------------".format(RNAi))
        top_hits, seed, top_hit_vals = get_best_match(RNAi, data)
        show_best_match_movies(top_hits, seed, top_hit_vals, saveFlag=True, recip_saveFlag=False)

def get_data():
    """
    Retrieves all pairwise information for sets of genes. for each gene pair: (rnai1 label, rnai2 label, distance,
    distance to wt, overlapping dimensions, csi and pad from either mysql or csv

    :return: data: a list of lists. Each lists consists of rnai1 label, rnai2 label, distance,
    distance to wt, overlapping dimensions, csi and pad
    """

    '''get data from mysql'''
    # return get_data_from_mysql()

    # '''get data from csv'''
    return get_data_from_csv()


def get_data_from_mysql():
    conn, cursor = initialize_db()
    sql = "SELECT CONCAT('EMBD', LPAD(rnai1, 4, '0')), CONCAT('EMBD', LPAD(rnai2, 4, '0')), dist, dist2wt, dims, " \
          "CSI, PAD FROM pair_lookup"  # LPAD Left PADs the INT rnai1 value with 'EMBD' and zeros, up to 4 digits

    cursor.execute(sql)
    dataset = list(cursor.fetchall())
    return dataset


def get_data_from_csv():
    dataset = readCSV(FILENAME_PAD, ',')[3:]
    return dataset

def update_pad_rank(cursor, conn, data):
    """
    adds genes to the pad_rank database
    :param cursor: mysqldb cursor
    :data: data set containing gene pairs, PAD, distance, CSI values generated by get_data(), which
    either populates data from csv or mysql
    :return: None
    """
    from db_utils_embryos import update_row, insert_row
    table = 'pad_rank'
    columns = get_columns(cursor, table)  # retrieves columns from table as list
    for i in range(1, 504):
        print(i)
        top_hits = PADFilter(i, data, 51)
        count = 0
        colms = {}
        for col in columns[1:]:
            colms[col] = 0 #add a default value to populate key value pair
        for hit in top_hits:
            rnai = i
            hit_name = int(hit[0][4:])
            pad = hit[1]
            count+=1
            rank_col = 'rank'+'{}'.format(count)
            pad_col = 'pad'+'{}'.format(count)
            colms['rnai']=rnai
            colms[rank_col]=hit_name
            colms[pad_col]=pad

        sql = "SELECT id FROM {t} WHERE rnai={r}".format(t=table, r=i)
        cursor.execute(sql)
        vals = cursor.fetchall()
        if len(vals) > 0:
            update_row(int(vals[0][0]), colms, cursor, table)
        elif len(vals) == 0:
            insert_row(colms, cursor, table)
    conn.commit()
    conn.close()

def convert_tifstack2mp4(path, save_folder):
    '''
    Takes in a folder of multitif files and saves them as mp4 files. For each multitif stack, a temp folder is generated
    and the multitif is split into individual tifs (using cv2). These tifs are then converted to mp4 via a command line
    bash, employing ffmpeg.
    :param path: path to folder containing tif files i.e. 'Z:/Fiji_Combined/'
    :param save_folder: '50_gene_mp4_files/'
    :return: generates .mp4 files and saves to a separate folder
    '''
    import pims
    import cv2
    import subprocess
    import os
    from myFunc import clearFolder
    import time

    # path = 'Z:/Fiji_Combined/'
    # file = 'pha-4_EMBD0002_GLS.tif'
    # save_folder = '50_gene_mp4_files/'

    # files = [f for f in os.listdir(path)]

    for file in os.listdir(path):
        if file[-3:] == 'tif':
            filename_in = path+file
            name = file[:-4]

            path_out = path+save_folder+name+'.mp4'
            # path_out = path+save_folder+'test.mp4'
            print('starting {0}'.format(name))
            # print('path out is {0}'.format(path_out))

            data = pims.TiffStack(filename_in) #brings in multitif stack (i.e. 31,903,448,3)
            # dims = data.frame_shape[1:3]

            if not os.path.exists(path+'tmp'): os.makedirs(path+'tmp')  #makes temp folder for individual tifs
            else:
                print('file exist, clearing folder', path+'tmp')
                # clearFolder(path+'tmp') #for some reason clearing this folder doesnt work properly- commenting it out works around bug

            count = 0
            for frame in data[0]:  # writes individual tif files to tmp folder
                rgb_frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
                form_count = '{:02d}'.format(count)
                cv2.imwrite(path+'tmp/'+name+form_count+'.tif', rgb_frame)
                # print(path+'tmp/'+name+form_count+'.tif')
                count+=1

            path_to_tifs = path+'tmp/'+name+'%02d.tif'
            # path_to_tifs = path+'tmp/'+name+'.tif'

            print(path_to_tifs)
            # bash_command = "ffmpeg -f image2 -i {0} -r 4 -pix_fmt yuv420p -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" {1}".format(path_to_tifs,path_out)
            bash_command = "ffmpeg -f image2 -framerate 7 -i {0} -r 21 -pix_fmt yuv420p -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" {1}".format(path_to_tifs,path_out)

            # -f image 2 (tells it there are multiple images), -r 7 (frame rate), -pix_fmt...(formats pixel color that is QT supported), -vf... makes sure image height is not odd value
            print(bash_command)
            subprocess.Popen(bash_command, shell=True)

def lookup_yanai_cluster(embdN):
    '''
    Takes in embd number (int) and returns developmental gene expression cluster number from Yanai lab 2015 Nature paper
    :param embdN: EMBD number, i.e. 321
    :return: expression broad cluster number, fine cluster number
    '''
    from embdFunc import getGeneName

    data_broad = readCSV('Z:/Yanai_broad_clusters.csv', ',')
    data_fine = readCSV('Z:/Yanai_fine_clusters.csv', ',')
    gene = getGeneName(embdN)
    geneN = getGeneCommonName(embdN)
    # print(gene, geneN)

    data_broad = np.array(data_broad)
    data_fine = np.array(data_fine)
    # print(data_broad)[:, 1]

    inds_cl_b = np.where(geneN == data_broad[:,1])[0]
    if len(inds_cl_b) is 0:
        # print('GENE NOT PRESENT IN YANAI SET')
        return 'NA'
    else:
        inds_cl_b = inds_cl_b[0]
        # print(inds_cl_b)
        inds_cl_f = np.where(geneN == data_fine[:,1])[0][0]
        # print(inds_cl_f)
        if inds_cl_f or inds_cl_b:
            cl_broad = np.int(np.float(data_broad[inds_cl_b,0]))
            cl_fine = np.int(np.float(data_fine[inds_cl_f,0]))

            print('{0},{1},{2},{3}').format(embdN, geneN, cl_broad, cl_fine)
            return cl_broad, cl_fine

def check_yanai_cluster_by_EMBDgroup():
    '''
    checks EMBD genes in manual groups against the Yanai 2015 Nature paper gene expression groupings
    :return: prints broad cluster and fine cluster assignment for each gene present in both data sets
    '''
    rnaiList = [[21, 364, 408, 130], [264, 386, 31, 32, 255, 388, 422, 118, 359], [417, 64, 115, 281, 77],
                [63, 19, 501, 77], [181, 182, 357], [184, 185, 154, 363, 117, 45, 447, 108], [34, 16, 31, 264],
                [9, 398, 435, 45, 52, 38],
                [95, 4, 57, 5, 277], [67, 90, 235, 403, 503, 261, 404, 25],
                [26, 453, 235, 327, 489, 379, 420, 225, 289, 261, 350],
                [10, 15, 217, 18, 28, 439, 177, 291, 209], [110, 142, 98, 101],
                # note that 91 and 186 were removed from dev delay because these were originally put in the 'test' group and 91 is bad
                [383, 495, 498, 414, 375, 396, 321], [387, 385, 422, 31], [320, 58, 125, 288],
                [426, 197, 1, 76]]  # 320 group is sim in gls but morph is dif. 387 is wnt group, 426 is no/low markers

    labels = ['eCFS_R', 'sect rupt', 'eMEX', 'lMEX', 'lCFS_R', 'hiNeur', 'dorsRupt', 'dorsBend', 'rupture', 'gRupt',
              'gMix', '2xCrunch', 'dev delay', 'wt', 'CFS_wnt', 'CFS_lGlY', 'no_markers']

    yanai_clusterlist = []
    for j in range(len(rnaiList)):
        group = rnaiList[j]
        label = labels[j]
        yanai_clusterlist.append(label)
        for condition in group:
            print('condition = {0}').format(condition)
            cl_b, cl_f = lookup_yanai_cluster(condition)
            yanai_clusterlist.append((cl_b,cl_f))
    print(yanai_clusterlist)

def get_germlayer_yanai(embdN):
    from embdFunc import getGeneName

    data = readCSV('Z:/germlayer_lists_yanai.csv', ',')[1:]
    data = np.array(data)
    gene = getGeneName(embdN)
    geneN = getGeneCommonName(embdN)
    # print(gene, geneN)

    inds_endo = np.where(geneN == data[:, 0])[0]
    inds_ecto = np.where(geneN == data[:, 1])[0]
    inds_meso = np.where(geneN == data[:, 2])[0]
    # print(inds_endo,inds_ecto,inds_meso)
    if len(inds_endo) is 0 and len(inds_ecto) is 0 and len(inds_meso) is 0:
        # print('GENE NOT PRESENT IN YANAI SET')
        return 'NA'
    elif len(inds_endo)>0:
        print('{0},{1}, Endoderm').format(embdN, geneN)
        return 'ENDO'
    elif len(inds_ecto)>0:
        print('{0},{1}, Ectoderm').format(embdN, geneN)
        return 'ECTO'
    elif len(inds_meso)>0:
        print('{0},{1}, Mesoderm').format(embdN, geneN)
        return 'MESO'

def output_best_matches_to_csv(RNAi):
    '''input- RNAi object- i.e EMBD0076
    Finds top matches using PAD (orCSI) by calling a CSV file, returns a list of ints [12,33,126] representing the genes with the closest phenotypes'''

    seed = RNAi[-4:]
    seedObj = RNAiClass(int(seed))

    data_opt = readCSV('Z:/phenoDist_rank_04142022.csv', ',')[3:]
    data_GLS_params_only = readCSV('Z:/phenoDist_rank_GLS_ONLY05102022.csv', ',')[3:]
    data_MS_params_only = readCSV('Z:/phenoDist_rank_MS_ONLY05102022.csv', ',')[3:]
    data_all_params = readCSV('Z:/phenoDist_rank_ALL_PARAMS05102022.csv', ',')[3:]


    # top_hits_opt = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
    top_hit_data_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
    top_hit_vals_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 1]  # retrieves the PAD values for top matches

    # top_hits_GLS = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
    top_hit_data_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
    top_hit_vals_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 1]  # retrieves the PAD values for top matches

    # top_hits_MS = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
    top_hit_data_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
    top_hit_vals_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 1]  # retrieves the PAD values for top matches

    # top_hits_all_params = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
    top_hit_data_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
    top_hit_vals_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 1]  # retrieves the PAD values for top matches

    all_hits_matrix = [['rank', 'query unique_ID', 'query genename', 'hit unique_ID - optimal parameters', 'hit genename - optimal parameters', 'hit PAD value - optimal parameters', 'hit unique_ID - GLS parameters', 'hit genename - GLS parameters', 'hit PAD value - GLS parameters', 'hit unique_ID - MS parameters', 'hit genename - MS parameters', 'hit PAD value - MS parameters', 'hit unique_ID - all parameters', 'hit genename - all parameters', 'hit PAD value - all parameters']]
    for i in range(len(top_hit_data_opt)):
        rank = i + 1
        query = seed
        query_name = getGeneCommonName(seed)
        opt_hit =top_hit_data_opt[i][-4:]
        opt_hit_name = getGeneCommonName(opt_hit)
        opt_hit_val = top_hit_vals_opt[i]
        GLS_hit =top_hit_data_GLS[i][-4:]
        GLS_hit_name = getGeneCommonName(GLS_hit)
        GLS_hit_val = top_hit_vals_GLS[i]
        MS_hit =top_hit_data_MS[i][-4:]
        MS_hit_name = getGeneCommonName(MS_hit)
        MS_hit_val = top_hit_vals_MS[i]
        all_hit =top_hit_data_all_params[i][-4:]
        all_hit_name = getGeneCommonName(all_hit)
        all_hit_val = top_hit_vals_all_params[i]
        all_hits_matrix.append([rank, query, query_name, opt_hit, opt_hit_name,opt_hit_val, GLS_hit, GLS_hit_name,GLS_hit_val, MS_hit, MS_hit_name,MS_hit_val, all_hit, all_hit_name,all_hit_val])



        # all_hits_matrix.append(top_hit_data_opt[i][-4:])  # slices the last four digits of string to remove the 'EMBD'
        # hit_name = getGeneCommonName(top_hit_opt[-4:])
    all_hits_matrix = np.array(all_hits_matrix)
    print(all_hits_matrix)

    import pandas as pd
    savepath = 'Z://{0}_hits_matrix.csv'.format(RNAi)
    pd.DataFrame(all_hits_matrix).to_csv(savepath)

    return(all_hits_matrix)
    # print top_hits, seed, top_hit_vals
    # return top_hits, seed, top_hit_vals

def output_all_best_matches_to_csv_combined():

    '''
    for all conditions, finds top matches using PAD, writes to a CSV file: shows top matches with alternate parameter sets'''
    import time

    fd = open('Z:\\all_best_matches_{d}.csv'.format(d=time.strftime('%m%d%Y')),
              'w')
    for i in range(1,504):
        RNAi = 'EMBD{r:04}'.format(r=i)
        seed = RNAi[-4:]
        seedObj = RNAiClass(int(seed))

        data_opt = readCSV('Z:/phenoDist_rank_04142022.csv', ',')[3:]
        data_GLS_params_only = readCSV('Z:/phenoDist_rank_GLS_ONLY05102022.csv', ',')[3:]
        data_MS_params_only = readCSV('Z:/phenoDist_rank_MS_ONLY05102022.csv', ',')[3:]
        data_all_params = readCSV('Z:/phenoDist_rank_ALL_PARAMS05102022.csv', ',')[3:]


        # top_hits_opt = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
        top_hit_data_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_vals_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 1]  # retrieves the PAD values for top matches

        # top_hits_GLS = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
        top_hit_data_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_vals_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 1]  # retrieves the PAD values for top matches

        # top_hits_MS = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
        top_hit_data_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_vals_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 1]  # retrieves the PAD values for top matches

        # top_hits_all_params = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
        top_hit_data_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_vals_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 1]  # retrieves the PAD values for top matches
        fd.write('rank, query EMBD_ID, query genename, hit EMBD_ID - optimal parameters, hit genename - optimal parameters, hit PAD value - optimal parameters, hit EMBD_ID - GLS parameters, hit genename - GLS parameters, hit PAD value - GLS parameters, hit EMBD_ID - MS parameters, hit genename - MS parameters, hit PAD value - MS parameters+ \n')
        all_hits_matrix = [['rank', 'query EMBD_ID', 'query genename', 'hit EMBD_ID - optimal parameters', 'hit genename - optimal parameters', 'hit PAD value - optimal parameters', 'hit EMBD_ID - GLS parameters', 'hit genename - GLS parameters', 'hit PAD value - GLS parameters', 'hit EMBD_ID - MS parameters', 'hit genename - MS parameters', 'hit PAD value - MS parameters', 'hit EMBD_ID - all parameters', 'hit genename - all parameters', 'hit PAD value - all parameters']]
        for i in range(len(top_hit_data_opt)):
            rank = i + 1
            query = seed
            query_name = getGeneCommonName(seed)
            opt_hit =top_hit_data_opt[i][-4:]
            opt_hit_name = getGeneCommonName(opt_hit)
            opt_hit_val = top_hit_vals_opt[i]
            GLS_hit =top_hit_data_GLS[i][-4:]
            GLS_hit_name = getGeneCommonName(GLS_hit)
            GLS_hit_val = top_hit_vals_GLS[i]
            MS_hit =top_hit_data_MS[i][-4:]
            MS_hit_name = getGeneCommonName(MS_hit)
            MS_hit_val = top_hit_vals_MS[i]
            all_hit =top_hit_data_all_params[i][-4:]
            all_hit_name = getGeneCommonName(all_hit)
            all_hit_val = top_hit_vals_all_params[i]
            all_hits_matrix.append([rank, query, query_name, opt_hit, opt_hit_name,opt_hit_val, GLS_hit, GLS_hit_name,GLS_hit_val, MS_hit, MS_hit_name,MS_hit_val, all_hit, all_hit_name,all_hit_val])
            fd.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}\n'.format(rank, query, query_name, opt_hit, opt_hit_name, opt_hit_val, GLS_hit, GLS_hit_name, GLS_hit_val, MS_hit, MS_hit_name, MS_hit_val))
            if i==9:
                fd.write('\n')

            # all_hits_matrix.append(top_hit_data_opt[i][-4:])  # slices the last four digits of string to remove the 'EMBD'
            # hit_name = getGeneCommonName(top_hit_opt[-4:])
        # all_hits_matrix = np.array(all_hits_matrix)
        # print(all_hits_matrix)


        # import pandas as pd
        # savepath = 'Z:/WebContent/{0}/{0}_hits_matrix.csv'.format(RNAi)
        # pd.DataFrame(all_hits_matrix).to_csv(savepath)
    fd.close()
        # return(all_hits_matrix)
        # print top_hits, seed, top_hit_vals
        # return top_hits, seed, top_hit_vals

def output_all_best_matches_to_one_csv_combined():
    '''
    IN PROGRESS__for all conditions, finds top matches using PAD, writes to one CSV file for all conditions: shows top matches  and PAD (hit-1: 0.765) with alternate parameter sets'''
    import time

    # all_hits_matrix = [['query EMBD_ID', 'query genename', 'hit EMBD_ID - optimal parameters',
    #                     'hit genename - optimal parameters', 'hit PAD value - optimal parameters',
    #                     'hit EMBD_ID - GLS parameters', 'hit genename - GLS parameters',
    #                     'hit PAD value - GLS parameters', 'hit EMBD_ID - MS parameters', 'hit genename - MS parameters',
    #                     'hit PAD value - MS parameters', 'hit EMBD_ID - all parameters',
    #                     'hit genename - all parameters', 'hit PAD value - all parameters']]

    all_hits_matrix = [['query EMBD_ID', 'query genename', 'hits - optimal parameters', 'hits - GLS only parameters', 'hits - MS only parameters',
                        'hits - all parameters']]

    # fd = open('Z:\\Automated_analysis\\combined_best_matches_{d}.csv'.format(d=time.strftime('%m%d%Y')),
    #           'w')
    # fd.write(
    #     'query EMBD_ID, query genename, hit EMBD_ID - optimal parameters, hit genename - optimal parameters, hit PAD value - optimal parameters, hit EMBD_ID - GLS parameters, hit genename - GLS parameters, hit PAD value - GLS parameters, hit EMBD_ID - MS parameters, hit genename - MS parameters, hit PAD value - MS parameters, hit EMBD_ID - All parameters, hit genename - All parameters, hit PAD value - All parameters+ \n')

    for i in range(1,504):
        RNAi = 'EMBD{r:04}'.format(r=i)
        print(RNAi)
        seed = RNAi[-4:]
        seedObj = RNAiClass(int(seed))
        query = seed
        query_name = getGeneCommonName(seed)

        data_opt = readCSV('Z:/phenoDist_rank_04142022.csv', ',')[3:]
        data_GLS_params_only = readCSV('Z:/phenoDist_rank_GLS_ONLY05102022.csv', ',')[3:]
        data_MS_params_only = readCSV('Z:/phenoDist_rank_MS_ONLY05102022.csv', ',')[3:]
        data_all_params = readCSV('Z:/phenoDist_rank_ALL_PARAMS05102022.csv', ',')[3:]


        # top_hits_opt = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
        top_hit_data_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        # top_hit_data_opt_gene_names = []

        # for j in range(len(top_hit_data_opt)):
        #     hit_name = getGeneCommonName(top_hit_data_opt[j][-4:])
        #     top_hit_data_opt_gene_names.append(hit_name)

        top_hit_data_opt_gene_names= [getGeneCommonName(top_hit_data_opt[j][-4:]) for j in range(len(top_hit_data_opt))]
        top_hit_vals_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 1]  # retrieves the PAD values for top matches
        top_hit_opt_set = []
        for k in range(10):
            val = float(top_hit_vals_opt[k])
            if val>0.50:
                val = str(val)[0:5]
                top_hit_opt_set.append(top_hit_data_opt_gene_names[k]+': '+ val)

        top_hit_data_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_data_GLS_gene_names = [getGeneCommonName(top_hit_data_GLS[j][-4:]) for j in
                                       range(len(top_hit_data_GLS))]
        top_hit_vals_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 1]  # retrieves the PAD values for top matches
        top_hit_GLS_set = []
        for l in range(10):
            val = float(top_hit_vals_GLS[l])
            if val > 0.50:
                val = str(val)[0:5]
                top_hit_GLS_set.append(top_hit_data_GLS_gene_names[l] + ': ' + val)

        top_hit_data_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_data_MS_gene_names = [getGeneCommonName(top_hit_data_MS[j][-4:]) for j in
                                       range(len(top_hit_data_MS))]
        top_hit_vals_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 1]  # retrieves the PAD values for top matches
        top_hit_MS_set = []
        for m in range(10):
            val = float(top_hit_vals_MS[m])
            if val > 0.50:
                val = str(val)[0:5]
                top_hit_MS_set.append(top_hit_data_MS_gene_names[m] + ': ' + val)

        top_hit_data_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_data_all_gene_names = [getGeneCommonName(top_hit_data_all_params[j][-4:]) for j in
                                       range(len(top_hit_data_all_params))]
        top_hit_vals_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 1]  # retrieves the PAD values for top matches
        top_hit_all_set = []
        for n in range(10):
            val = float(top_hit_vals_all_params[n])
            if val > 0.50:
                val = str(val)[0:5]
                top_hit_all_set.append(top_hit_data_all_gene_names[n] + ': ' + val)

        # fd.write(
        #     '{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}\n'.format(query, query_name, top_hit_data_opt.tolist(), top_hit_data_opt_gene_names,
        #                                                                  top_hit_vals_opt.tolist(), top_hit_data_GLS.tolist(), top_hit_data_GLS_gene_names,
        #                                                                  top_hit_vals_GLS.tolist(), top_hit_data_MS.tolist(), top_hit_data_MS_gene_names, top_hit_vals_MS.tolist(), top_hit_data_all_params.tolist(), top_hit_data_all_gene_names, top_hit_vals_all_params.tolist()))

        # all_hits_matrix.append(
        #     [query, query_name, top_hit_data_opt.tolist(),top_hit_data_opt_gene_names,top_hit_vals_opt.tolist(),top_hit_data_GLS.tolist(),top_hit_data_GLS_gene_names,top_hit_vals_GLS.tolist(),top_hit_data_MS.tolist(),top_hit_data_MS_gene_names, top_hit_vals_MS.tolist(),top_hit_data_all_params.tolist(),top_hit_data_all_gene_names,top_hit_vals_all_params.tolist()])

        all_hits_matrix.append([query, query_name, top_hit_opt_set, top_hit_GLS_set, top_hit_MS_set,top_hit_all_set])

    import pandas as pd
    savepath = 'Z:/hits_matrix_v2_{d}.csv'.format(d=time.strftime('%m%d%Y'))
    pd.DataFrame(all_hits_matrix).to_csv(savepath)



        # all_hits_matrix = [['rank', 'query EMBD_ID', 'query genename', 'hit
        # all_hits_matrix = [['rank', 'query EMBD_ID', 'query genename', 'hit EMBD_ID - optimal parameters', 'hit genename - optimal parameters', 'hit PAD value - optimal parameters', 'hit EMBD_ID - GLS parameters', 'hit genename - GLS parameters', 'hit PAD value - GLS parameters', 'hit EMBD_ID - MS parameters', 'hit genename - MS parameters', 'hit PAD value - MS parameters', 'hit EMBD_ID - all parameters', 'hit genename - all parameters', 'hit PAD value - all parameters']]
        # for i in range(len(top_hit_data_opt)):
        #     rank = i + 1
        #     query = seed
        #     query_name = getGeneCommonName(seed)
        #     opt_hit =top_hit_data_opt[i][-4:]
        #     opt_hit_name = getGeneCommonName(opt_hit)
        #     opt_hit_val = top_hit_vals_opt[i]
        #     GLS_hit =top_hit_data_GLS[i][-4:]
        #     GLS_hit_name = getGeneCommonName(GLS_hit)
        #     GLS_hit_val = top_hit_vals_GLS[i]
        #     MS_hit =top_hit_data_MS[i][-4:]
        #     MS_hit_name = getGeneCommonName(MS_hit)
        #     MS_hit_val = top_hit_vals_MS[i]
        #     all_hit =top_hit_data_all_params[i][-4:]
        #     all_hit_name = getGeneCommonName(all_hit)
        #     all_hit_val = top_hit_vals_all_params[i]
        #     all_hits_matrix.append([rank, query, query_name, opt_hit, opt_hit_name,opt_hit_val, GLS_hit, GLS_hit_name,GLS_hit_val, MS_hit, MS_hit_name,MS_hit_val, all_hit, all_hit_name,all_hit_val])
        #     fd.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}\n'.format(rank, query, query_name, opt_hit, opt_hit_name, opt_hit_val, GLS_hit, GLS_hit_name, GLS_hit_val, MS_hit, MS_hit_name, MS_hit_val))
        #     if i==9:
        #         fd.write('\n')

            # all_hits_matrix.append(top_hit_data_opt[i][-4:])  # slices the last four digits of string to remove the 'EMBD'
            # hit_name = getGeneCommonName(top_hit_opt[-4:])
        # all_hits_matrix = np.array(all_hits_matrix)
        # print(all_hits_matrix)


        # import pandas as pd
        # savepath = 'Z:/WebContent/{0}/{0}_hits_matrix.csv'.format(RNAi)
        # pd.DataFrame(all_hits_matrix).to_csv(savepath)
    # fd.close()
        # return(all_hits_matrix)
        # print top_hits, seed, top_hit_vals
        # return top_hits, seed, top_hit_vals


def custom_updates_notes():
    seed = '1905'
    top_hits = ['1905','364', '326', '408', '115', '21', '130', '1903']
    top_hit_vals = ['0.67993', '0.63764', '0.61267', '0.58633', '0.57857', '0.56956', '0.55368']
    show_best_match_movies(top_hits, seed, top_hit_vals, saveFlag=True, recip_saveFlag=False)

    seed = '1903'
    top_hits = ['1903', '478', '359', '141', '89', '118', '13', '422', '115', '408']
    top_hit_vals = ['0.66098', '0.64050', '0.63266', '0.62765', '0.59764', '0.59212', '0.59015', '0.57881',
                    '0.57621']
    show_best_match_movies(top_hits, seed, top_hit_vals, saveFlag=True, recip_saveFlag=False)

    seed = '1904'
    top_hits = ['1904', '1', '426', '76', '281', '64', '434', '365']
    top_hit_vals = ['0.92980', '0.85668', '0.81408', '0.76214', '0.62551', '0.58568', '0.55706']
    show_best_match_movies(top_hits, seed, top_hit_vals, saveFlag=True, recip_saveFlag=False)

    seed = '0703'
    top_hits = ['703', '885', '364', '408', '1905', '130', '326', '77', '118', '1903']
    top_hit_vals = ['0.82921', '0.77982', '0.77140', '0.70873', '0.69040', '0.63525', '0.58980', '0.57257', '0.57257']
    show_best_match_movies(top_hits, seed, top_hit_vals, saveFlag=True, recip_saveFlag=False)

    seed = '0885'
    top_hits = ['885', '703', '364', '408', '1905', '326', '130', '1360', '21', '359']
    top_hit_vals = ['0.82921', '0.77310', '0.73323', '0.70939', '0.63357', '0.62595', '0.58663', '0.57472', '0.57277']
    show_best_match_movies(top_hits, seed, top_hit_vals, saveFlag=True, recip_saveFlag=False)

    seed = '1406'
    top_hits = ['1406', '674', '6', '440', '133', '441', '45', '9', '357', '55']
    top_hit_vals = ['0.72692', '0.69791', '0.68559', '0.67087', '0.66879', '0.66488', '0.63687', '0.63102',
                    '0.62884']
    show_best_match_movies(top_hits, seed, top_hit_vals, saveFlag=True, recip_saveFlag=False)

    seed = '0674'
    top_hits = ['674', '6', '1406', '397', '952', '1360', '356', '9', '133', '357']
    top_hit_vals = ['0.74037', '0.72692', '0.71087', '0.69899', '0.67270', '0.66870', '0.65705', '0.65645',
                    '0.65542']
    show_best_match_movies(top_hits, seed, top_hit_vals, saveFlag=True, recip_saveFlag=False)

    seed = '0952'
    top_hits = ['952', '866', '674', '57', '6', '441', '435', '137', '133', '45']
    top_hit_vals = ['0.73255', '0.69899', '0.65246', '0.65169', '0.64878', '0.63979', '0.63593', '0.63434', '0.63263']
    show_best_match_movies(top_hits, seed, top_hit_vals, saveFlag=True, recip_saveFlag=False)

    seed = '0866'
    top_hits = ['866', '952', '16', '32', '6', '386', '31', '674', '264', '435']
    top_hit_vals = ['0.73255', '0.66751', '0.66485', '0.65980', '0.65376', '0.64882', '0.63847', '0.63841', '0.63562']
    show_best_match_movies(top_hits, seed, top_hit_vals, saveFlag=True, recip_saveFlag=False)

    seed = '1360'
    top_hits = ['1360', '674', '356', '32', '866', '422']
    top_hit_vals = ['0.67270', '0.65182', '0.61506', '0.58917', '0.58917']
    show_best_match_movies(top_hits, seed, top_hit_vals, saveFlag=True, recip_saveFlag=False)

def initialize():
    '''run'''

    '''make tophits movie'''
    # data = get_data()
    # print(data[:10])
    # make_movies(data)  # define EMBD range inside function


    # print('here')
    # print PADFilter('EMBD{r:04}'.format(r=1), data, 3)
    # print PADFilter(1, data, 3)
    '''update PAD rank'''
    # conn, curs = initialize_db()
    # data = get_data()
    # update_pad_rank(curs, conn, data)

    '''make mp4 movie'''
    # path = 'Z:/Fiji_Combined/fix/'
    # save_folder = '50_gene_mp4_files/'
    # convert_tifstack2mp4(path, save_folder)

    '''get top hits'''
    # print PADFilter('EMBD{r:04}'.format(r=417), data, 10)
    # print(getGeneCommonName(388))


    "make csv of best matches"
    # for r in range(1,504):
    #     RNAi = 'EMBD{r:04}'.format(r=r)
    #     output_best_matches_to_csv(RNAi)

    # output_all_best_matches_to_csv_combined()
    output_all_best_matches_to_one_csv_combined()


    # pilot_list = range(1,504)
    # pilot_list = [2,3,4,5,6,7,9,12,13,14,15,16,18,22,25,26,27,29,38,39,41,42,43,52,53,54,55,57,62,63,68,74,76,77,80,82,385,386,387,388]
    # yanai_list = []
    # for i in pilot_list:
    #     cl_b, cl_f = lookup_yanai_cluster(i)
    #     # print(i,cl_b,cl_f)
    #     if cl_b != 'N':
    #         layer = get_germlayer_yanai(i)
    #         yanai_list.append((i,getGeneCommonName(i),layer,cl_b,cl_f))
    #         # print('...processing')
    # yanai_list =np.array(yanai_list)
    # print(yanai_list[np.where(yanai_list[:,2]== 'ENDO')[0]])
    # print(yanai_list[np.where(yanai_list[:,2]== 'MESO')[0]])
    # print(yanai_list[np.where(yanai_list[:,2]== 'ECTO')[0]])
    # print(yanai_list[np.where(yanai_list[:,2]== 'NA')[0]])
    # print(len(yanai_list))
        # get_germlayer_yanai(i)

    # check_yanai_cluster_by_EMBDgroup()
    # get_germlayer_yanai(10)
    # lookup_yanai_cluster(10)

    # crunchlist = [10, 15, 217, 18, 28, 439, 177, 291, 209]
    # for i in crunchlist:
    #     print(getGeneCommonName(i))

    # data = get_data_from_csv()
    # seedList = [398,320,440,9,435]
    # getCommonPairs(seedList, data)


if __name__ == '__main__':
    # find_cliques_by_seed()
    initialize()

    '''gets PAD value for gene pair'''
    # r1= RNAiClass(5)
    # print r1.paramsMS['tMove']
    #     r2=RNAiClass(95)
    #     print r1.getDiffParams(r2)
    #
    #     fileName = 'Z:\\Automated_analysis\\phenoDistCalc\\phenoDist_Spots.csv'
    #     data = readCSV(fileName, ',')[3:]


    # embs = getGLSEmbs()
    # ts_list = []
    # for emb in embs:
    #     ts_list.append(emb.params['tScale'])
    # ts_list = np.array(ts_list)
    # # emb = embs[np.argmin(np.abs(1-ts_list))]
    # diff = np.abs(1-ts_list)
    # emb_ind = np.argsort(diff)
    # for i in emb_ind[:11]:
    #     print i, embs[i].label, embs[i].date, embs[i].params['tScale']

    # print(emb.params['tScale'], emb.date, emb.label)

    '''vet groups and identify new group members'''
    # fileName = 'Z:\\Automated_analysis\\phenoDistCalc\\phenoDist_07242017.csv'
    # fileName = 'Z:\\Automated_analysis\\phenoDistCalc\\phenoDist_rank_02282020.csv'
    # data = readCSV(fileName, ',')[3:]
    # seedList = [383, 495]
    # getCommonPairs(seedList, data)

#     seedList = range(1,504)
#     for l in seedList:
#         r1 = RNAiClass(l)
# #         if r1.paramsGLS['aY']> 128:
#         print(l, 'MS aR', r1.paramsMS['aR'],  'GLS aG', r1.paramsGLS['aG'], 'GLS aR', r1.paramsGLS['aR'], 'GLS aY', r1.paramsGLS['aY'])

#     seed = 397
#     seedList = [34, 16, 31, 264, 9, 398, 435, 45, 52, 38, 95, 4, 57, 5, 277]
#     data = get_data()
#     print PADFilter('EMBD{r:04}'.format(r=356), data, 10)
#     for seed in seedList:
#         getPairsParams(seed, data)


#     rnaiList = [[19,31,63,77,118,181,264,357], [9,397,357,181,53],[5,90,274,268,435,255,264], [7,52,54,390,447,398],[10,217,291,412,15,59,288,398],[177,195,239,330,61,241,433,462,493,390]]
#     for rn in rnaiList:
#         for r in rn:
#             print(PADFilter('EMBD{r:04}'.format(r=r), data, 10))

#     geneList = [getGeneName(i) for i in range(2,504)]
#     interact = getInteractions(geneList)
#     allInterGenes = interact.ravel()
#     subData = []
#     for d in data:
# #         if int(d[0][-4:]) in allInterGenes and int(d[1][-4:]) in allInterGenes: subData.append([int(d[0][-4:]), int(d[1][-4:]), d[2]]) #distance
#         if int(d[0][-4:]) in allInterGenes and int(d[1][-4:]) in allInterGenes: subData.append([int(d[0][-4:]), int(d[1][-4:]), d[6]]) #PAD
# #         if int(d[0][-4:]) in allInterGenes and int(d[1][-4:]) in allInterGenes: subData.append([int(d[0][-4:]), int(d[1][-4:]), d[5]]) #CSI
#     subData = np.array(subData)
#     
# #     actual = np.random.choice(2,subData.shape[0])#Random
# 
#     actual = np.zeros(subData.shape[0])
#     for pair in interact:
#         ind0 = np.where(subData[:,0]==pair[0])[0]
#         ind1 = np.where(subData[:,1]==pair[1])[0]
#         ind = np.intersect1d(ind0, ind1)
#         actual[ind]=1
#         ind0 = np.where(subData[:,0]==pair[1])[0]
#         ind1 = np.where(subData[:,1]==pair[0])[0]
#         ind = np.intersect1d(ind0, ind1)
#         actual[ind]=1
#     f1,rec,pr = [], [], []
#     maxD = max(subData[:,2][np.where(np.isfinite(subData[:,2]))])
# #     x = np.arange(0.,maxD,maxD/100.) #distance
#     x = np.arange(0.,1.,0.01)
#     for pth in x:
# #         pred = (subData[:,2]<=pth) #distance
#         pred = (subData[:,2]>=pth)
#         f1.append(metrics.f1_score(actual,pred))
#         rec.append(metrics.recall_score(actual,pred))
#         pr.append(metrics.precision_score(actual,pred))
#     f1 = np.array(f1)
#     figF = myFigure()
#     figF.title('F1 metric')
#     figF.plot(x, f1)
#     figF.xlabel('pad threshold')
#     figP = myFigure()
#     figP.title('precision')
#     figP.plot(x, pr)
#     figP.xlabel('pad threshold')
#     figR = myFigure()
#     figR.title('recall')
#     figR.plot(x, rec)
#     figR.xlabel('pad threshold')
#     padTh = x[np.argmax(f1)]
# #     pred = (subData[:,2]>=padTh) #distance
#     pred = (subData[:,2]>=padTh)
#     print('for PAD threshold = {p} F1 metric performance: {f}'.format(p=padTh, f=max(f1)))
#     print('precision = {p}\nrecall = {r}'.format(p=metrics.precision_score(actual,pred),r=metrics.recall_score(actual,pred)))
#     figF.show()
