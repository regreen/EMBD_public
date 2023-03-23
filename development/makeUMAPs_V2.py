import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import umap
import pandas as pd

all_params= False  #if false, runs optimal params dataset
dist = False  #runs UMAP on Euclidean dist. array
pad = True  #runs UMAP on PAD array
run_query = True


for i in range(1,2):
    if run_query:
        query = i
        gene_list = pd.read_csv('Z:/gene_key.csv', delimiter=',',encoding='latin-1')
        gene_list = pd.DataFrame(gene_list)
        gene_list = np.array(gene_list)
        gene_list_b = list(gene_list)
        row = gene_list[np.where(gene_list[:,0]=='EMB_P{n:04}'.format(n=query))]
        row_list = list(row)
        common_name = row_list[0][2]
        print(common_name)

    sns.set(style='white', context='poster', rc={'figure.figsize':(14,10)})
    np.random.seed(71076)
    # np.random.seed(92507)

    rnaiList = [i for i in range(1,504)] # a list of EMBD numbers to generate RNAi objects from
    labels = [str(x) for x in rnaiList]



    # print(labels)
    if dist:
        if all_params:
            data = np.load('Z:/Dist_Matrix_042122_allParams.npy')
        else:
            data = np.load('Z:/Dist_Matrix_042122_optParams.npy')
    elif pad:
        if all_params:
            data = np.load('Z://PAD_Matrix_042122_allParams.npy')
        else:
            data = np.load('Z://PAD_Matrix_042122_optParams.npy')


    groupColors = groupColors = ['sienna', 'orangered', 'firebrick', 'darkorange', 'gold', 'y', 'green', 'turquoise', 'teal',
                       'dodgerblue', 'b']
    plot_labels = ['Early Arrest, \n High Ectoderm', 'Early Arrest, \n High Mesoderm',
                       'Mid Arrest, \n High Mesoderm',
                       'Early Arrest Sectored  \n Germ-Layers', 'Elongation Defect,  \n High Neurons',
                       'Enclosure/Early Elong.  \n Dorsal Defect', 'Enclosure/Early Elong. \n Ventral Rupture',
                       'Variable \n Defects', 'Late Elongation \n Arrest', 'Dev. Delay', 'Wild-type-like']
    # groupLists = [[1,2], [3,4], [18], [5], #test list-- shorter than full dataset
    #                 [6], [7,8,9], [10],
    #                 [11],
    #                 [12, 13, 14], [15],
    #                 [16]]
    # groupLists = [[21], [31, 32], [64,77],[63, 19],
    #                 [45], [9, 45, 52, 38], [95, 4, 57, 5],
    #                 [26],
    #                 [10, 15, 18, 28], [98, 91],
    #                 [72]]#use this set to test out without having to wait for full dataset

    groupLists = [[21, 364, 408, 130], [417, 64, 115, 281, 77], [77, 63, 19, 501],
                      [264, 386, 31, 32, 255, 388, 422, 118, 359],
                      [181, 182, 357, 184, 185, 154, 363, 117, 45, 447, 108], [9, 398, 435, 45, 52, 38],
                      [95, 4, 57, 5, 277],
                      [26, 453, 235, 327, 489, 379, 420, 225, 289, 261, 350],
                      [10, 15, 217, 18, 28, 439, 177, 291, 209], [110, 142, 98, 101, 186],
                      [383, 495, 498, 414, 375, 396, 321]]

    if run_query:
        query = [[query]]
        groupLists += query  # adds the query onto the last position in the list
        groupLabels = ['eCFS_R',  'eMEX', 'lMEX', 'sect rupt', 'hiNeur', 'dorsBend', 'rupture', 'gMix', '2xCrunch', 'dev delay', 'wt', 'query']
        groupColors = ['sienna',  'orangered','firebrick','darkorange' ,'gold', 'y', 'green', 'turquoise', 'teal', 'dodgerblue', 'b', 'deeppink']
        plot_labels = ['Early Arrest, \n High Ectoderm', 'Early Arrest, \n High Mesoderm', 'Mid Arrest, \n High Mesoderm',
                   'Early Arrest Sectored  \n Germ-Layers', 'Elongation Defect,  \n High Neurons',
                  'Enclosure/Early Elong.  \n Dorsal Defect', 'Enclosure/Early Elong. \n Ventral Rupture',
                  'Variable \n Defects', 'Late Elongation \n Arrest', 'Dev. Delay', 'Wild-type-like', 'Query']

    fit = umap.UMAP()
    points = fit.fit_transform(data)

    fig = plt.figure(figsize=(12,12))
    for point in points:
        plt.scatter(point[0], point[1], color="grey", marker=".", alpha=0.5)

    for g in range(len(groupLists)):
        man_group_list = groupLists[g] #individual manual group
        lab = plot_labels[g]
        inds1= [] # this will be populated with the index position for each member of the man group within the array
        for p in man_group_list:
            inds1.append(labels.index(str(p)))
    #     print(inds1)

        for ind in inds1:
            plt.scatter(points[ind][0], points[ind][1], color=groupColors[g], marker=".",
                            alpha=0.5) #plots the individual members of a group in a specific color
            plt.scatter(np.mean(points[inds1], axis=0)[0], np.mean(points[inds1], axis=0)[1], label=lab, marker="*",
                            alpha=0.5, color=groupColors[g])
            if lab=='Query':
                query_label = '{n}(RNAi)'.format(n=common_name)
    #             query_label = 'EMBD{n:04}(RNAi)'.format(n=query[0][0])
                plt.text(np.mean(points[inds1]+0.5, axis=0)[0], np.mean(points[inds1], axis=0)[1], query_label,
                         color=groupColors[g], size=14, style = 'italic', weight='bold')
            elif lab=='Late Elongation \n Arrest':
                plt.text(np.mean(points[inds1], axis=0)[0], np.mean(points[inds1]-0.6, axis=0)[1], lab,
                         color=groupColors[g], size=11)
            elif lab=='Mid Arrest, \n High Mesoderm':
                plt.text(np.mean(points[inds1], axis=0)[0]+0.1, np.mean(points[inds1]+0.3, axis=0)[1], lab,
                         color=groupColors[g], size=11)
            elif (g % 2)==0:
                plt.text(np.mean(points[inds1]-0.2, axis=0)[0], np.mean(points[inds1]-0.5, axis=0)[1], lab,
                         color=groupColors[g], size=11)
            else:
                plt.text(np.mean(points[inds1]+0.2, axis=0)[0], np.mean(points[inds1]+0.1, axis=0)[1], lab,
                         color=groupColors[g], size=11)

    if all_params:
        if dist:
            plt.title('UMAP of Euclidean distance with ALL parameters')
            plt.savefig('Z:/UMAP_EMBD_042222_Dist_ALL.svg')
        elif pad:
            plt.title('UMAP of PAD with ALL parameters')
            plt.savefig('Z:/UMAP_EMBD_042222_PAD_ALL.svg')
    else:
        if dist:
            plt.title('UMAP of Euclidean distance with optimal parameters')
            plt.savefig('Z:/UMAP_EMBD_042222_Dist_optv2.svg')
        elif pad:
            plt.title('UMAP of PAD with optimal parameters')
            plt.savefig('Z:/UMAP_EMBD_042222_PAD_optv2_EMBD{n:04}.svg'.format(n=query[0][0]))
    plt.close()