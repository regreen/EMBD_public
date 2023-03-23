"""
Created on Feb 28, 2017
This file contains programs to make figures for EMBD paper 2
@author: Becky
"""

import matplotlib.pyplot as plt
from Embryos import MSEmbryo, GSEmbryo
from varLookup import FOLDER_IN
import numpy as np
from embdFunc import getGeneCommonName
from RNAiClass import RNAiClass, ControlClass
import pandas as pd
import matplotlib


def make_MS_fig_with_subplots(embd, embryo, saveFlag=False, showFlag=False):
    """
    Takes in EMBD number and embryo number and makes a figure consisting of subplots showing all measured features for a
    given embryo
    :param embd: EMBD number (i.e. 0423)
    :param embryo: embryo number (i.e. 3)
    :param saveFlag keyword enables saving to Z://MS_embs
    :param showFlag keyword enables display of each embryo individually. Set to False when plotting multiple embs on 1
     plot, since "show" call resides outside of loop.
    :return:
    """
    no_yscale = True  # allows you to plot without y-scale (if true, plots without yscale-- if false, plots with y scale)
    if no_yscale:
        print("----Plotting WITHOUT Y-Scale----")
    folders = [FOLDER_IN + 'cropped/EMBD{n:04}/MS/Emb{x}/'.format(n=embd, x=embryo)]
    for folder in folders[:]:
        emb = MSEmbryo(folder, check_version=False)  # loads embryo object

    subplots = ['tIntR', 'tIntG', 'headInt', 'lengthR', 'MoI0G', 'MoI1G', 'MoI0R', 'MoI1R', 'CoM0G', 'CoM1G', 'CoM0R',
                'CoM1R']  # defines curves to plot
    labels = [r'Total Intensity Neurons $(I^R_{tot}$)', r'Total Intensity Epidermis $(I^G_{tot}$)', r'Area of Head Epidermis $(\bar A^G_{head})$', r'Neuronal Length $(L^R)$',
              r'Epidermal Shape- Short $(\bar K^G_{short})$',
              r'Epidermal Shape- Long $(\bar K^G_{long})$', r'Neuronal Shape- Short $(\bar K^R_{short})$',
              r'Neuronal Shape- Long $(\bar K^R_{long})$',
              r'Epidermal Position- Short $(\bar D^G_{short})$', r'Epidermal Position- Long $(\bar D^G_{long})$',
              r'Neuronal Position- Short $(\bar D^R_{short})$',
              r'Neuronal Position- Long $(\bar D^R_{long})$']
    # labels = ['Total Intensity Neurons (red)', 'Total Intensity Epidermis (green)', 'Head Intensity', 'Neuronal Length',
    #           'Epidermal Shape- Short Axis (MoI0G)',
    #           'Epidermal Shape- Long Axis (MoI1G)', 'Neuronal Shape- Short Axis (MoI0R)',
    #           'Neuronal Shape- Long Axis (MoI1R)',
    #           'Epidermal Position- Short Axis (Com0G)', 'Epidermal Position- Long Axis (Com1G)',
    #           'Neuronal Position- Short Axis (Com0R)',
    #           'Neuronal Position- Long Axis (Com1R)']
    xtick_labels = ['-200', '-100', '0', '100', '200', '300', '400', '500', '600', '700']
    columns = 3
    if len(subplots) % columns > 0:  # rounds up to include sufficient rows to accomodate all plots
        rows = len(subplots) / columns + 1
    else:
        rows = len(subplots) / columns

    if plt.get_fignums():  # checks to see if figure exists already, so multiple embs can be plotted at once
        fig = plt.figure(num=1)
    else:
        fig = plt.figure(figsize=(12, 12))  # creates figure

    name = getGeneCommonName(embd)
    fig.suptitle("{g}(RNAi)\n EMBD{n:04} Emb{x}  Morphogenesis Strain".format(g=name, n=embd, x=embryo), fontsize=12, style='italic',
                 y=.97)  # adds title to the top of all subplots

    for i in range(len(subplots)):  # toggles through subplot index and plots as it goes
        time_thres = 550
        plt.subplot(rows, columns, i + 1)  # defines the subplot to draw on
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
        plt.title(labels[i], fontsize=12, y=1.05)  # makes subplot title
        plt.xlabel('Time (minutes)', fontsize=9)  # labels x-axis
        ax = plt.subplot(rows, columns, i + 1)
        ax.tick_params(axis='x', bottom=True, top=False, left=False, right=False, labelsize=8) # get rid of ticks on top and right
        ax.tick_params(axis='y', bottom=False, top=False, left=True, right=False, labelsize=8)

        ax.yaxis.get_offset_text().set_fontsize(8)
        x, y = emb.getCurve(subplots[i])
        x = 20 * x
        if no_yscale:
            y=y/emb.scale[subplots[i]]

        x = x[np.where(x < time_thres)]
        y = y[np.where(x < time_thres)]
        plt.plot(x, y, 'r', label=labels[i])  # plots data in red

        plt.plot(x[:emb.cutoff[subplots[i]]], y[:emb.cutoff[subplots[i]]], 'k--',
                 dashes=(3, 2))  # plots fit-region with dashed line
        # xt = np.arange(x[0] / emb.tScale + emb.t0, x[-1] / emb.tScale + emb.t0, 0.1)
        # # fig.plot((xt - self.t0) * self.tScale,
        # #          self.scale['headInt'] * cost(xt, self.params['aSigHead'], self.params['mSigHead'],
        # #                                       self.params['sSigHead'], \
        # #                                       self.params['aGaussHead'], self.params['mGaussHead'],
        # #                                       self.params['sGaussHead']))
        if subplots[i] in emb.avgCurves:
            xR, yR, yerrR = emb.avgCurves[subplots[i]]
            xR = xR * 20
            xR = xR[np.where(xR < time_thres)]  # makes the plots without late timepoints
            yR = yR[np.where(xR < time_thres)]
            yerrR = yerrR[np.where(xR < time_thres)]
            upper = np.add(yR, yerrR)
            lower = np.subtract(yR, yerrR)
        plt.plot(xR, yR, 'k')
        plt.fill_between(xR, upper, lower, color='0.8')  # plots with shaded error bars
        # plt.errorbar(xR[::10], yR[::10], yerrR[::10], color='k')  # plots with normal looking error bars
        plt.subplots_adjust(wspace=0.35, hspace=0.5)
    if saveFlag:
        plt.savefig('Z:/EMBD{n:04}_Emb{x}_MS.svg'.format(n=embd, x=embryo))
    if showFlag:
        plt.show()


def make_GLS_fig_with_subplots(embd, embryo, saveFlag=False, showFlag=False):
    """
    Takes in EMBD number and embryo number and makes a figure consisting of subplots showing all measured features for a
    given embryo
    :param embd: EMBD number (i.e. 0423)
    :param embryo: embryo number (i.e. 3)
    :return:
    """
    no_yscale = True  # allows you to plot without yscale
    if no_yscale:
        print("----Plotting WITHOUT Y-Scale----")
    folders = [FOLDER_IN + 'cropped/EMBD{n:04}/GLS/Emb{x}/'.format(n=embd, x=embryo)]
    for folder in folders[:]:
        emb = GSEmbryo(folder, check_version=False)  # loads embryo object
    curves = emb.getCurveNames()
    print(curves)

    subplots = ['spotsG', 'spotsR', 'spotsY', 'MoI0G', 'MoI1G', 'MoI0R', 'MoI1R', 'MoI0Y', 'MoI1Y', 'CoM0R', 'CoM1R',
                'CoM0Y', 'CoM1Y']  # defines curves to plot
    labels = [r'Nuclear Count- Endoderm/Pharynx $(N^G_{nuc})$', r'Nuclear Count- Ectoderm $(N^R_{nuc})$',
              r'Nuclear Count- Mesoderm $ (N^Y_{nuc})$', r'Endoderm/Pharynx Shape- Short $(\bar K^G_{short})$',
              r'Endoderm/Pharynx Shape- Long $(\bar K^G_{long})$', r'Ectoderm Shape- Short $(\bar K^R_{short})$',
              r'Ectoderm Shape- Long $(\bar K^R_{long})$',
              r'Mesoderm Shape- Short $(\bar K^Y_{short})$', r'Mesoderm Shape- Long $(\bar K^Y_{long})$',
              r'Ectoderm Position - Short $(\bar D^R_{short})$',
              r'Ectoderm Position- Long $(\bar D^R_{long})$', r'Mesoderm Position - Short $(\bar D^Y_{short})$',
              r'Mesoderm Position- Long $(\bar D^Y_{long})$']

    # labels = ['Number of Endoderm/Pharynx Nuclei (green)', 'Number of Ectoderm Nuclei (red)',
    #           'Number of Mesoderm Nuclei (yellow)', 'Endo/Pharynx Shape- Short Axis (MoI0G)',
    #           'Endo/Pharynx Shape- Long Axis (MoI1G)', 'Ectoderm Shape- Short Axis (MoI0R)',
    #           'Ectoderm Shape- Long Axis (MoI1R)',
    #           'Mesoderm Shape- Short Axis (MoI0Y)', 'Mesoderm Shape- Long Axis (MoI1Y)',
    #           'Ectoderm Position - Short Axis (Com0R)',
    #           'Ectoderm Position- Long Axis (Com1R)', 'Mesoderm Position - Short Axis (Com0Y)',
    #           'Mesoderm Position- Long Axis (Com1Y)']

    xtick_labels = ['-200', '-100', '0', '100', '200', '300', '400', '500', '600', '700']
    columns = 3

    if len(subplots) % columns > 0:  # rounds up to include sufficient rows to accommodate all plots
        rows = len(subplots) / columns + 1
    else:
        rows = len(subplots) / columns

    if plt.get_fignums():
        fig = plt.figure(num=1)
    else:
        fig = plt.figure(figsize=(12, 12))  # creates figure
    name = getGeneCommonName(embd)
    fig.suptitle("{g}(RNAi)\n EMBD{n:04} Emb{x}  Germ Layer Strain".format(g=name, n=embd, x=embryo), fontsize=12, style='italic',
                 y=.99)  # adds title to the top of all subplots

    Gcurves = ['MoI0G', 'MoI1G']  # this block is in place to allow for dropping skewed values before markers turn on.
    Rcurves = ['MoI0R', 'MoI1R', 'CoM0R', 'CoM1R']
    Ycurves = ['MoI0Y', 'MoI1Y', 'CoM0Y', 'CoM1Y']
    fixGcurves = False
    fixGinds = 0
    fixRcurves = False
    fixRinds = 0
    fixYcurves = False
    fixYinds = 0


    for i in range(len(subplots)):  # toggles through subplot index and plots as it goes
        time_thres = 550 # switch back to 550
        plt.subplot(rows, columns, i + 1)  # defines the subplot to draw on
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)
        plt.title(labels[i], fontsize=12, y=1.00)  # makes subplot title

        plt.xlabel('Time (minutes)', fontsize=9)  # labels x-axis
        ax = plt.subplot(rows, columns, i + 1)
        # ax.tick_params(axis='x', labelsize=8)
        # ax.tick_params(axis='y', labelsize=8)
        ax.tick_params(axis='x', bottom=True, top=False, left=False, right=False, labelsize=8) # get rid of ticks on top and right
        ax.tick_params(axis='y', bottom=False, top=False, left=True, right=False, labelsize=8)
        ax.yaxis.get_offset_text().set_fontsize(8)

        x, y = emb.getCurve(subplots[i])
        x = 20 * x
        if no_yscale:
            y = y/emb.scale[subplots[i]]
        x = x[np.where(x < time_thres)]
        y = y[np.where(x < time_thres)]


        if subplots[i] in ['spotsG', 'spotsR','spotsY']:  # this is a fix put in place to deal with many plots where the early timepoints had massively skewed spots counts (i.e 1000)
            if np.any(y[0:10]>400):
                inds = np.where(y[0:10]>400)
                y = pd.DataFrame(y)
                y.iloc[[inds]] = np.nan
                if subplots[i] == 'spotsG':
                    fixGcurves = True
                    fixGinds = inds
                elif subplots[i] == 'spotsR':
                    fixRcurves = True
                    fixRinds = inds
                elif subplots[i] == 'spotsY':
                    fixYcurves = True
                    fixYinds = inds
        if subplots[i] in Gcurves and fixGcurves:
            y = pd.DataFrame(y)
            y.iloc[fixGinds] = np.nan
        if subplots[i] in Rcurves and fixRcurves:
            y = pd.DataFrame(y)
            y.iloc[fixRinds] = np.nan
        if subplots[i] in Ycurves and fixYcurves:
            y = pd.DataFrame(y)
            y.iloc[fixYinds] = np.nan



        plt.plot(x, y, 'r', label=labels[i])  # plots data in red
        plt.plot(x[:emb.cutoff[subplots[i]]], y[:emb.cutoff[subplots[i]]], 'k--',
                 dashes=(3, 2))  # plots fit-region with dashed line
        # xt = np.arange(x[0] / emb.tScale + emb.t0, x[-1] / emb.tScale + emb.t0, 0.1)
        # # fig.plot((xt - self.t0) * self.tScale,
        # #          self.scale['headInt'] * cost(xt, self.params['aSigHead'], self.params['mSigHead'],
        # #                                       self.params['sSigHead'], \
        # #                                       self.params['aGaussHead'], self.params['mGaussHead'],
        # #                                       self.params['sGaussHead']))
        if subplots[i] in emb.avgCurves:
            xR, yR, yerrR = emb.avgCurves[subplots[i]]
            xR = xR * 20
            xR = xR[np.where(xR < time_thres)]  # makes the plots without late timepoints
            yR = yR[np.where(xR < time_thres)]
            yerrR = yerrR[np.where(xR < time_thres)]
            upper = np.add(yR, yerrR)
            lower = np.subtract(yR, yerrR)



        plt.plot(xR, yR, 'k')
        plt.fill_between(xR, upper, lower, color='0.8')  # plots with shaded error bars
        # plt.errorbar(xR[::10], yR[::10], yerrR[::10], color='k')  # plots with normal looking error bars
        plt.subplots_adjust(wspace=0.35, hspace=0.5)
    if saveFlag:
        plt.savefig(
            'Z:/EMBD{n:04}_Emb{x}_GLS.svg'.format(n=embd, x=embryo))
    if showFlag:
        plt.show()


def plot_allEmbs_allSubplots(EMBD, strain='GLS', saveAllFlag=False, showFlag=False):
    from fileLookup import FILENAME_WEB_CONTENT
    import os
    '''
    plots all used individual embryos on figure with subplots for every measured feature
    :param strain: 'GLS' or 'MS' callable by keyword
    :param EMBD: RNAi#. i.e. 315
    :param saveAllFlag: keyword param to save figure of all used embryos relative to controls
    :return: makes a figure with subplots
    '''
    # from RNAiClass import RNAiClass


    if EMBD == 0:
        r = ControlClass()
        r.defineFolders()
        r.getEmbryos()
        r.asignUse()
    else:
        r = RNAiClass(EMBD)
        r.getEmbryos()
        r.asignUse()

    usedEmbsMS = r.MSUseEmbs  # retrieves a list of embryo objects that are used (considered by automated analysis)
    usedEmbsGLS = r.GLSUseEmbs
    # allMovingEmbsGLS = r.GLSMoveEmbs
    # allNoMoveEmbsGLS = r.GLSNoMoveEmbs
    # allMovingEmbsMS = r.MSMoveEmbs
    # allNoMoveEmbsMS = r.MSNoMoveEmbs
    name = getGeneCommonName(EMBD)
    print(EMBD)
    if strain == 'GLS':
        for emb in usedEmbsGLS:
            print('plotting {0}'.format(emb.folder.split('/')[4]))
            emb_number = int(emb.folder.split('/')[4][3:])  # retrieves the folder number (i.e. embryo number)
            make_GLS_fig_with_subplots(EMBD, emb_number, saveFlag=True, showFlag=False)
        plt.suptitle("{g}(RNAi)\n EMBD{n:04} Germ Layer Strain (all embryos)".format(g=name, n=EMBD), fontsize=14,
                     style='italic', y=.97)  # adds title to the top of all subplots
        if showFlag:
            plt.show()
        if saveAllFlag:
            folder = FILENAME_WEB_CONTENT + 'EMBD{n:04}'.format(n=EMBD)
            if not os.path.exists(folder):
                os.makedirs(folder)
            plt.savefig(
                # 'Z:/GLS_embs/EMBD{n:04}_all_GLS.svg'.format(n=EMBD))
                FILENAME_WEB_CONTENT + 'EMBD{n:04}/EMBD{n:04}_all_GLS.svg'.format(n=EMBD))
                # FILENAME_WEB_CONTENT + 'EMBD{n:04}/EMBD{n:04}_all_GLS_yScaled.svg'.format(n=EMBD))


            plt.close()
    if strain == 'MS':
        for emb in usedEmbsMS:
            print('plotting {0}'.format(emb.folder.split('/')[4]))
            emb_number = int(emb.folder.split('/')[4][3:])  # retrieves the folder number (i.e. embryo number)
            make_MS_fig_with_subplots(EMBD, emb_number, saveFlag=True, showFlag=False)
        plt.suptitle("{g}(RNAi)\n EMBD{n:04} Morphogenesis Strain (all embryos)".format(g=name, n=EMBD), fontsize=14,
                     style='italic', y=.97)  # adds title to the top of all subplots
        if saveAllFlag:
            folder = folder = FILENAME_WEB_CONTENT + 'EMBD{n:04}'.format(n=EMBD)
            if not os.path.exists(folder):
                os.makedirs(folder)
            plt.savefig(FILENAME_WEB_CONTENT + 'EMBD{n:04}/EMBD{n:04}_all_MS.svg'.format(n=EMBD))
            # plt.savefig(FILENAME_WEB_CONTENT + 'EMBD{n:04}/EMBD{n:04}_all_MS_yScaled.svg'.format(n=EMBD))
        if showFlag:
            plt.show()
        plt.close()


def get_adjusted_average_curve(EMBD, curveName, strain):
    '''Takes in RNAi number and curveName and returns average curve (x, y) and error for plotting
    :param curveName is curve name, such as 'spotsG'
    :param EMBD is EMBD number, i.e. 0423
    :param strain: either 'GLS' or 'MS'
    '''
    no_yscale = False #allows you to generate values without y scaling
    r = RNAiClass(EMBD)
    r.getEmbryos()
    r.asignUse()
    if strain == 'GLS':
        usedEmbs = r.GLSUseEmbs  # retrieves a list of embryo objects that are used (considered by automated analysis)
    elif strain == 'MS':
        usedEmbs = r.MSUseEmbs  # retrieves a list of embryo objects that are used (considered by automated analysis)

    x_all = []
    y_all = []
    # figure = plt.figure()
    for emb in usedEmbs:
        x, y = emb.getCurve(curveName)
        x_all.append(x)
        if no_yscale:
            y = y / emb.scale[curveName]
            y_all.append(y)
        else:
            y_all.append(y)
        # print y[0:9]
        # print x[0:9]
        # plt.plot(x,y)

    y_avg = np.nanmean(y_all, axis=0)
    y_std = np.nanstd(y_all, axis=0)
    x_avg = np.nanmean(x_all, axis=0)
    x_all = emb.time

    '''plot data'''
    # figure = plt.figure()
    # plt.plot(x_avg, y_avg, 'r')
    # upper = np.add(y_avg, y_std)
    # lower = np.subtract(y_avg, y_std)
    # plt.fill_between(x_avg, upper, lower, color='0.8')  # plots with shaded error bars
    # plt.show()

    return (x_avg, y_avg, y_std)


def plot_avgRNAi_control_allSubplots(embd, strain='GLS', saveFlag=False, showFlag=False):
    """
    Takes in EMBD number & plots a figure consisting of subplots showing all measured features (avg) vs. control
    :param embd: EMBD number (i.e. 0423)
    :return:
    """

    r = RNAiClass(embd)
    r.getEmbryos()
    r.asignUse()

    if strain == 'GLS':
        folders = r.GLSUseEmbs  # retrieves a list of embryo objects that are used (considered by automated analysis)
        subplots = ['spotsG', 'spotsR', 'spotsY', 'MoI0G', 'MoI1G', 'MoI0R', 'MoI1R', 'MoI0Y', 'MoI1Y', 'CoM0R',
                    'CoM1R', 'CoM0Y', 'CoM1Y']  # defines curves to plot
        labels = ['Number of Endoderm/Pharynx Nuclei (green)', 'Number of Ectoderm Nuclei (red)',
                  'Number of Mesoderm Nuclei (yellow)', 'Endo/Pharynx Shape- Short Axis (MoI0G)',
                  'Endo/Pharynx Shape- Long Axis (MoI1G)', 'Ectoderm Shape- Short Axis (MoI0R)',
                  'Ectoderm Shape- Long Axis (MoI1R)',
                  'Mesoderm Shape- Short Axis (MoI0Y)', 'Mesoderm Shape- Long Axis (MoI1Y)',
                  'Ectoderm Position - Short Axis (Com0R)',
                  'Ectoderm Position- Long Axis (Com1R)', 'Mesoderm Position - Short Axis (Com0Y)',
                  'Mesoderm Position- Long Axis (Com1Y)']

    elif strain == 'MS':
        folders = r.MSUseEmbs
        subplots = ['tIntR', 'tIntG', 'headInt', 'lengthR', 'MoI0G', 'MoI1G', 'MoI0R', 'MoI1R', 'CoM0G', 'CoM1G',
                    'CoM0R', 'CoM1R']  # defines curves to plot
        labels = ['Total Intensity Neurons (red)', 'Total Intensity Epidermis (green)', 'Head Intensity',
                  'Neuronal Length', 'Epidermal Shape- Short Axis (MoI0G)',
                  'Epidermal Shape- Long Axis (MoI1G)', 'Neuronal Shape- Short Axis (MoI0R)',
                  'Neuronal Shape- Long Axis (MoI1R)',
                  'Epidermal Position- Short Axis (Com0G)', 'Epidermal Position- Long Axis (Com1G)',
                  'Neuronal Position- Short Axis (Com0R)',
                  'Neuronal Position- Long Axis (Com1R)']

    for emb_obj in folders:
        curves = emb_obj.getCurveNames()

    xtick_labels = ['-200', '-100', '0', '100', '200', '300', '400', '500', '600', '700']
    columns = 3

    if len(subplots) % columns > 0:  # rounds up to include sufficient rows to accommodate all plots
        rows = len(subplots) / columns + 1
    else:
        rows = len(subplots) / columns

    if plt.get_fignums():
        fig = plt.figure(num=1)
    else:
        fig = plt.figure(figsize=(12, 12))  # creates figure
    name = getGeneCommonName(embd)
    fig.suptitle("{g}(RNAi)\n EMBD{n:04} Average {s}".format(g=name, n=embd, s=strain), fontsize=14, style='italic',
                 y=.97)  # adds title to the top of all subplots

    for i in range(len(subplots)):  # toggles through subplot index and plots as it goes
        plt.subplot(rows, columns, i + 1)  # defines the subplot to draw on
        plt.title(labels[i], fontsize=12, y=1.05)  # makes subplot title
        plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
        ax = plt.subplot(rows, columns, i + 1)
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        ax.yaxis.get_offset_text().set_fontsize(8)
        x, y, err = get_adjusted_average_curve(embd, subplots[i],
                                               strain)  # finds the average curve from all considered embryos
        x = 20 * x  # adjusts x axis to time rather than frame
        plt.plot(x, y, 'r', label=labels[i])  # plots data in red
        upper_RNAi = np.add(y, err)
        lower_RNAi = np.subtract(y, err)
        plt.fill_between(x, upper_RNAi, lower_RNAi, color='darksalmon', alpha=0.5)  # plots with shaded error bars

        if subplots[i] in emb_obj.avgCurves:
            xR, yR, yerrR = emb_obj.avgCurves[subplots[i]]
            xR = xR * 20
            upper = np.add(yR, yerrR)
            lower = np.subtract(yR, yerrR)
        plt.plot(xR, yR, 'k')
        plt.fill_between(xR, upper, lower, color='0.8', alpha=0.5)  # plots with shaded error bars
        # plt.errorbar(xR[::10], yR[::10], yerrR[::10], color='k')  # plots with normal looking error bars
        plt.subplots_adjust(wspace=0.35, hspace=0.5)
    if saveFlag:
        if strain == 'GLS':
            plt.savefig(
                'Z:/EMBD{n:04}_avg_{s}.svg'.format(n=embd, s=strain))
        elif strain == 'MS':
            plt.savefig(
                'Z:/EMBD{n:04}_avg_{s}.svg'.format(n=embd, s=strain))

    if showFlag:
        plt.show()


def plot_3color_spots(EMBD):
    '''
    plots all spots data for a given RNAi condition on a single plot
    :param EMBD: EMBD number (i.e.0423)
    :return: plot
    '''

    showfit = True
    curves = ['spotsG', 'spotsR', 'spotsY']

    if EMBD == 0:
        r = RNAiClass(1)  # this is a work around- just need an embryo object to call to access the control avg curves
        r.getEmbryos()
        r.asignUse()
        c = ControlClass()
        folders = r.GLSUseEmbs  # retrieves a list of embryo objects that are used (considered by automated analysis)
        for folder in folders:
            emb_obj = folder

        fig = plt.figure()
        plt.title("Control\n EMBD{n:04} Nuclear Count".format(n=EMBD), fontsize=12, style='italic',
                  y=1.02)  # adds title to the top of all subplots
        plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
        ax = plt.subplot()
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        ax.yaxis.get_offset_text().set_fontsize(8)

        for curve in curves:
            if curve in emb_obj.avgCurves:
                x, y, err = emb_obj.avgCurves[curve]
                x = x * 20
                if curve == 'spotsG':
                    plt.plot(x, y, 'g', label='Endoderm/Pharynx')  # plots data in red
                    plt.legend(loc=2, prop={'size': 8})
                    upper_RNAi = np.add(y, err)
                    lower_RNAi = np.subtract(y, err)
                    plt.fill_between(x, upper_RNAi, lower_RNAi, color='palegreen',
                                     alpha=0.5)  # plots with shaded error bars
                    if showfit:
                        plt.plot(x, sigmoidGLS(x, c.paramsGLS['aG'], c.paramsGLS['mG'] * 20, c.paramsGLS['sG'] * 20, 1),
                                 color='k')

                elif curve == 'spotsR':
                    plt.plot(x, y, 'r', label='Ectoderm')  # plots data in red
                    plt.legend(loc=2, prop={'size': 8})
                    upper_RNAi = np.add(y, err)
                    lower_RNAi = np.subtract(y, err)
                    plt.fill_between(x, upper_RNAi, lower_RNAi, color='darksalmon',
                                     alpha=0.5)  # plots with shaded error bars
                    if showfit:
                        plt.plot(x, sigmoidGLS(x, c.paramsGLS['aR'], c.paramsGLS['mR'] * 20, c.paramsGLS['sR'] * 20, 1),
                                 color='k')

                elif curve == 'spotsY':
                    plt.plot(x, y, 'gold', label='Mesoderm')  # plots data in red
                    plt.legend(loc=2, prop={'size': 8})
                    upper_RNAi = np.add(y, err)
                    lower_RNAi = np.subtract(y, err)
                    plt.fill_between(x, upper_RNAi, lower_RNAi, color='khaki',
                                     alpha=0.5)  # plots with shaded error bars
                    if showfit:
                        plt.plot(x, sigmoidGLS(x, c.paramsGLS['aY'], c.paramsGLS['mY'] * 20, c.paramsGLS['sY'] * 20, 1),
                                 color='k')

        plt.show()


    else:
        r = RNAiClass(EMBD)
        r.getEmbryos()
        r.asignUse()
        fig = plt.figure()
        name = getGeneCommonName(EMBD)
        plt.title("{g}(RNAi)\n EMBD{n:04} Nuclear Count".format(g=name, n=EMBD), fontsize=12, style='italic',
                  y=1.02)  # adds title to the top of all subplots
        plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
        ax = plt.subplot()
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        ax.yaxis.get_offset_text().set_fontsize(8)

        for curve in curves:
            x, y, err = get_adjusted_average_curve(EMBD, curve,
                                                   'GLS')  # finds the average curve from all considered embryos
            x = 20 * x  # adjusts x axis to time rather than frame

            if curve == 'spotsG':
                plt.plot(x, y, 'g', label='Endoderm/Pharynx')  # plots data in red
                plt.legend(loc=2, prop={'size': 8})
                upper_RNAi = np.add(y, err)
                lower_RNAi = np.subtract(y, err)
                plt.fill_between(x, upper_RNAi, lower_RNAi, color='palegreen',
                                 alpha=0.5)  # plots with shaded error bars
                if showfit:
                    plt.plot(x, sigmoidGLS(x, r.paramsGLS['aG'], r.paramsGLS['mG'] * 20, r.paramsGLS['sG'] * 20, 1),
                             color='k')
            elif curve == 'spotsR':
                plt.plot(x, y, 'r', label='Ectoderm')  # plots data in red
                plt.legend(loc=2, prop={'size': 8})
                upper_RNAi = np.add(y, err)
                lower_RNAi = np.subtract(y, err)
                plt.fill_between(x, upper_RNAi, lower_RNAi, color='darksalmon',
                                 alpha=0.5)  # plots with shaded error bars
                if showfit:
                    plt.plot(x, sigmoidGLS(x, r.paramsGLS['aR'], r.paramsGLS['mR'] * 20, r.paramsGLS['sR'] * 20, 1),
                             color='k')

            elif curve == 'spotsY':
                plt.plot(x, y, 'gold', label='Mesoderm')  # plots data in red
                plt.legend(loc=2, prop={'size': 8})
                upper_RNAi = np.add(y, err)
                lower_RNAi = np.subtract(y, err)
                plt.fill_between(x, upper_RNAi, lower_RNAi, color='khaki',
                                 alpha=0.5)  # plots with shaded error bars
                if showfit:
                    plt.plot(x, sigmoidGLS(x, r.paramsGLS['aY'], r.paramsGLS['mY'] * 20, r.paramsGLS['sY'] * 20, 1),
                             color='k')

        plt.show()


def plot_2color_tInt(EMBD):
    '''
    plots green and red tInt avg data together for a given RNAi condition on a single plot
    :param EMBD: EMBD number (i.e.0423)
    :return: plot
    '''

    #note that to plot without y-scale go to get_adjusted_average_curve and change no_yscale to True

    showfit = True  # FIXME
    curves = ['tIntR', 'tIntG']

    if EMBD == 0:
        r = RNAiClass(1)  # this is a work around- just need an embryo object to call to access the control avg curves
        r.getEmbryos()
        r.asignUse()
        c = ControlClass()
        folders = r.MSUseEmbs  # retrieves a list of embryo objects that are used (considered by automated analysis)
        for folder in folders:
            emb_obj = folder

        fig = plt.figure()
        plt.title("Control\n EMBD{n:04} Total Intensity".format(n=EMBD), fontsize=12, style='italic',
                  y=1.02)  # adds title to the top of all subplots
        plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
        ax = plt.subplot()
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        ax.yaxis.get_offset_text().set_fontsize(8)

        for curve in curves:
            if curve in emb_obj.avgCurves:
                x, y, err = emb_obj.avgCurves[curve]
                x = x * 20
                if curve == 'tIntG':
                    plt.plot(x, y, 'g', label='Epidermis')  # plots data in red
                    plt.legend(loc=2, prop={'size': 8})
                    upper_RNAi = np.add(y, err)
                    lower_RNAi = np.subtract(y, err)
                    plt.fill_between(x, upper_RNAi, lower_RNAi, color='palegreen',
                                     alpha=0.5)  # plots with shaded error bars
                    if showfit:
                        plt.plot(x, sigmoidMS(x, c.paramsMS['aG'], c.paramsMS['bG'], c.paramsGLS['mG'] * 20,
                                              c.paramsMS['sG'] * 20,
                                              1), color='k')
                        # x / self.tScale, a, b, m / self.tScale, s, r) *self.scale[curveName],
                elif curve == 'tIntR':
                    plt.plot(x, y, 'r', label='Neurons')  # plots data in red
                    plt.legend(loc=2, prop={'size': 8})
                    upper_RNAi = np.add(y, err)
                    lower_RNAi = np.subtract(y, err)
                    plt.fill_between(x, upper_RNAi, lower_RNAi, color='darksalmon',
                                     alpha=0.5)  # plots with shaded error bars
                    if showfit:
                        plt.plot(x, sigmoidMS(x, c.paramsMS['aR'], c.paramsMS['bR'], c.paramsGLS['mR'] * 20,
                                              c.paramsMS['sR'] * 20,
                                              1), color='k')

        plt.show()


    else:
        r = RNAiClass(EMBD)
        r.getEmbryos()
        r.asignUse()
        fig = plt.figure()
        name = getGeneCommonName(EMBD)
        plt.title("{g}(RNAi)\n EMBD{n:04} Total Intensity".format(g=name, n=EMBD), fontsize=12, style='italic',
                  y=1.02)  # adds title to the top of all subplots
        plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
        ax = plt.subplot()
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        ax.yaxis.get_offset_text().set_fontsize(8)

        for curve in curves:
            x, y, err = get_adjusted_average_curve(EMBD, curve,
                                                   'MS')  # finds the average curve from all considered embryos
            x = 20 * x  # adjusts x axis to time rather than frame

            if curve == 'tIntG':
                plt.plot(x, y, 'g', label='Epidermis')  # plots data in red
                plt.legend(loc=2, prop={'size': 8})
                upper_RNAi = np.add(y, err)
                lower_RNAi = np.subtract(y, err)
                plt.fill_between(x, upper_RNAi, lower_RNAi, color='palegreen',
                                 alpha=0.5)  # plots with shaded error bars
                if showfit:
                    plt.plot(x, sigmoidMS(x, r.paramsMS['aG'], r.paramsMS['bG'], r.paramsGLS['mG'] * 20,
                                          r.paramsMS['sG'] * 20,
                                          1), color='k')
                    # plt.plot(x, sigmoidMS(x / r.tScale, r.paramsMS['aG'], r.paramsMS['bG'], r.paramsGLS['mG'] / r.tScale *20, r.paramsMS['sG']*20, 1) * r.scale[curve],color='k')
            elif curve == 'tIntR':
                plt.plot(x, y, 'r', label='Neurons')  # plots data in red
                plt.legend(loc=2, prop={'size': 8})
                upper_RNAi = np.add(y, err)
                lower_RNAi = np.subtract(y, err)
                plt.fill_between(x, upper_RNAi, lower_RNAi, color='darksalmon',
                                 alpha=0.5)  # plots with shaded error bars
                if showfit:
                    plt.plot(x, sigmoidMS(x, r.paramsMS['aR'], r.paramsMS['bR'], r.paramsGLS['mR'], r.paramsMS['sR'],
                                          1), color='k')
                    # plt.plot(x, sigmoidMS(x / r.tScale, r.paramsMS['aR'], r.paramsMS['bR'],
                    #                       r.paramsGLS['mR'] / r.tScale * 20, r.paramsMS['sR'] * 20, 1) * r.scale[curve],
                    #          color='k')

        plt.show()


def plot_color_MoI(EMBD, strain):
    '''
    plots all moment of inertia data for a given RNAi condition on plot
    :param EMBD: EMBD number (i.e.0423)
    :param strain: either MS or GLS
    :return: makes 2 plots (one for each axis)
    '''

    GLScurves = [['MoI0G', 'MoI0R', 'MoI0Y'], ['MoI1G', 'MoI1R', 'MoI1Y']]
    MScurves = [['MoI0G', 'MoI0R'], ['MoI1G', 'MoI1R']]
    labels = ['MoI- Short Axis', 'MoI- Long Axis']

    if EMBD == 0:
        r = RNAiClass(1)  # this is a work around- just need an embryo object to call to access the control avg curves
    elif EMBD > 0:
        r = RNAiClass(EMBD)
    r.getEmbryos()
    r.asignUse()
    GLSfolders = r.GLSUseEmbs  # retrieves a list of embryo objects that are used (considered by automated analysis)
    MSfolders = r.MSUseEmbs

    if strain == 'GLS':
        folders = GLSfolders
        curves = GLScurves
    elif strain == 'MS':
        folders = MSfolders
        curves = MScurves

    for folder in folders:
        emb_obj = folder

    fig = plt.figure()
    name = getGeneCommonName(EMBD)
    if EMBD == 0:
        fig.suptitle("Control\n EMBD{n:04} {s} Moment of Inertia".format(n=EMBD, s=strain), fontsize=12, style='italic',
                     y=.98)  # adds title to the top of all subplots
    elif EMBD > 0:
        fig.suptitle("{g}(RNAi)\n EMBD{n:04} {s} Moment of Inertia".format(n=EMBD, g=name, s=strain), fontsize=12,
                     style='italic',
                     y=.98)  # adds title to the top of all subplots
    plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
    ax = plt.subplot()
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.yaxis.get_offset_text().set_fontsize(8)

    for i in range(2):  # toggles through subplot index and plots as it goes
        plt.subplot(1, 2, i + 1)  # defines the subplot to draw on
        plt.title(labels[i], fontsize=12, y=.99)  # makes subplot title
        plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
        ax = plt.subplot(1, 2, i + 1)
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        ax.yaxis.get_offset_text().set_fontsize(8)
        for j in range(len(curves[i])):
            if EMBD == 0:
                x, y, err = emb_obj.avgCurves[curves[i][j]]
                x = 20 * x
                # plt.plot(x,y,'r', label = labels[i])
            elif EMBD > 0:
                x, y, err = get_adjusted_average_curve(EMBD, curves[i][j],
                                                       strain)  # finds the average curve from all considered embryos
                x = 20 * x  # adjusts x axis to time rather than frame
                # plt.plot(x,y,'r', label = labels[i])
            if curves[i][j][-1:] == 'G':
                if strain == 'GLS':
                    plt.plot(x, y, 'g', label='Endoderm/Pharynx')  # plots data in red
                elif strain == 'MS':
                    plt.plot(x, y, 'g', label='Epidermis')  # plots data in green
                plt.legend(loc=2, prop={'size': 8})
                upper_RNAi = np.add(y, err)
                lower_RNAi = np.subtract(y, err)
                plt.fill_between(x, upper_RNAi, lower_RNAi, color='palegreen',
                                 alpha=0.5)  # plots with shaded error bars
            elif curves[i][j][-1:] == 'R':
                if strain == 'GLS':
                    plt.plot(x, y, 'r', label='Ectoderm')  # plots data in red
                elif strain == 'MS':
                    plt.plot(x, y, 'r', label='Neurons')  # plots data in red
                plt.legend(loc=2, prop={'size': 8})
                upper_RNAi = np.add(y, err)
                lower_RNAi = np.subtract(y, err)
                plt.fill_between(x, upper_RNAi, lower_RNAi, color='darksalmon',
                                 alpha=0.5)  # plots with shaded error bars

            elif curves[i][j][-1:] == 'Y':
                plt.plot(x, y, 'gold', label='Mesoderm')  # plots data in yellow
                plt.legend(loc=2, prop={'size': 8})
                upper_RNAi = np.add(y, err)
                lower_RNAi = np.subtract(y, err)
                plt.fill_between(x, upper_RNAi, lower_RNAi, color='khaki',
                                 alpha=0.5)  # plots with shaded error bars
    plt.show()


def plot_color_CoM(EMBD, strain):
    '''
    plots all center of mass data for a given RNAi condition on a single plot
    :param EMBD: EMBD number (i.e.0423)
    :param strain: either MS or GLS
    :return: plot
    '''

    GLScurves = [['CoM0R'], ['CoM1R']]  # note that the CoMG curves are not used, and thus are not included
    MScurves = [['CoM0G', 'CoM0R'], ['CoM1G', 'CoM1R']]
    labels = ['CoM- Short Axis', 'CoM- Long Axis']

    if EMBD == 0:
        r = RNAiClass(1)  # this is a work around- just need an embryo object to call to access the control avg curves
    elif EMBD > 0:
        r = RNAiClass(EMBD)
    r.getEmbryos()
    r.asignUse()
    GLSfolders = r.GLSUseEmbs  # retrieves a list of embryo objects that are used (considered by automated analysis)
    MSfolders = r.MSUseEmbs

    if strain == 'GLS':
        folders = GLSfolders
        curves = GLScurves
    elif strain == 'MS':
        folders = MSfolders
        curves = MScurves

    for folder in folders:
        emb_obj = folder

    fig = plt.figure()
    name = getGeneCommonName(EMBD)
    if EMBD == 0:
        fig.suptitle("Control\n EMBD{n:04} {s} Center of Mass".format(n=EMBD, s=strain), fontsize=12, style='italic',
                     y=.98)  # adds title to the top of all subplots
    elif EMBD > 0:
        fig.suptitle("({g})RNAi\n EMBD{n:04} {s} Center of Mass".format(n=EMBD, g=name, s=strain), fontsize=12,
                     style='italic',
                     y=.98)  # adds title to the top of all subplots
    plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
    ax = plt.subplot()
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.yaxis.get_offset_text().set_fontsize(8)

    for i in range(2):  # toggles through subplot index and plots as it goes
        plt.subplot(1, 2, i + 1)  # defines the subplot to draw on
        plt.title(labels[i], fontsize=12, y=.99)  # makes subplot title
        plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
        ax = plt.subplot(1, 2, i + 1)
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        ax.yaxis.get_offset_text().set_fontsize(8)
        for j in range(len(curves[i])):
            if EMBD == 0:
                x, y, err = emb_obj.avgCurves[curves[i][j]]
                x = 20 * x
                # plt.plot(x,y,'r', label = labels[i])
            elif EMBD > 0:
                x, y, err = get_adjusted_average_curve(EMBD, curves[i][j],
                                                       strain)  # finds the average curve from all considered embryos
                x = 20 * x  # adjusts x axis to time rather than frame
                # plt.plot(x,y,'r', label = labels[i])
            if curves[i][j][-1:] == 'G':
                if strain == 'GLS':
                    plt.plot(x, y, 'g', label='Endoderm/Pharynx')  # plots data in red
                elif strain == 'MS':
                    plt.plot(x, y, 'g', label='Epidermis')  # plots data in green
                plt.legend(loc=2, prop={'size': 8})
                upper_RNAi = np.add(y, err)
                lower_RNAi = np.subtract(y, err)
                plt.fill_between(x, upper_RNAi, lower_RNAi, color='palegreen',
                                 alpha=0.5)  # plots with shaded error bars
            elif curves[i][j][-1:] == 'R':
                if strain == 'GLS':
                    plt.plot(x, y, 'r', label='Ectoderm')  # plots data in red
                elif strain == 'MS':
                    plt.plot(x, y, 'r', label='Neurons')  # plots data in red
                plt.legend(loc=2, prop={'size': 8})
                upper_RNAi = np.add(y, err)
                lower_RNAi = np.subtract(y, err)
                plt.fill_between(x, upper_RNAi, lower_RNAi, color='darksalmon',
                                 alpha=0.5)  # plots with shaded error bars

            elif curves[i][j][-1:] == 'Y':
                plt.plot(x, y, 'gold', label='Mesoderm')  # plots data in yellow
                plt.legend(loc=2, prop={'size': 8})
                upper_RNAi = np.add(y, err)
                lower_RNAi = np.subtract(y, err)
                plt.fill_between(x, upper_RNAi, lower_RNAi, color='khaki',
                                 alpha=0.5)  # plots with shaded error bars
    plt.show()


def plot_pretty_graphs(x_label, x_points, cv1_label, y_curve1, std_curve1, cv2_label=None, y_curve2=None,
                       std_curve2=None):
    '''
    plots python style graphs with shaded error bars, given input data
    :param x_label: a string i.e. 'Time(min)'
    :param x_points: a list of X-axis values [0,1,2,3]
    :param cv1_label: curve 1 label-- a string i.e. 'Neuronal Intensity'
    :param curve1: a list of numbers i.e. [0,1,2,4,7]
    :param std1: a list of numbers i.e. [0,1,2,4,7]
    :param cv2_label: curve 2 label-- a string i.e. 'Neuronal Intensity'
    :param curve2: a list of numbers i.e. [0,1,2,4,7]
    :param std2: a list of numbers i.e. [0,1,2,4,7]
    :return:
    '''

    fig = plt.figure()
    fig.suptitle("{0} vs.{1}".format(cv1_label, x_label), fontsize=14, style='italic',
                 y=.97)  # adds title to the top of all subplots

    plt.xlabel(x_label, fontsize=10)  # labels x-axis
    ax = plt.subplot()
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.yaxis.get_offset_text().set_fontsize(8)
    x, y, err = x_points, y_curve1, std_curve1
    plt.plot(x, y, '.g-', label=cv1_label)  # plots data in red
    upper = np.add(y, err)
    lower = np.subtract(y, err)
    plt.fill_between(x, upper, lower, color='palegreen',
                     alpha=0.5)  # plots with shaded error bars

    x2, y2, err2 = x_points, y_curve2, std_curve2
    plt.plot(x2, y2, '.r-', label=cv2_label)  # plots data in red
    upper = np.add(y2, err2)
    lower = np.subtract(y2, err2)
    plt.fill_between(x2, upper, lower, color='darksalmon',
                     alpha=0.5)  # plots with shaded error bars
    plt.show()


def show_500gene_heatmap():
    '''
    plots a heatmap readout of the defects scored in the 500 gene pilot data set from csv
    :return:
    '''
    import pandas as pd
    import seaborn as sns
    from string import ascii_letters
    import matplotlib.pyplot as plt

    # data = pd.read_csv(
    #     'Z:/EMBD500genepilotScoring012319_sorted.csv')
    data = pd.read_csv(
        'Z:/EMBD500genepilot_high_prev.csv')
    # data = pd.read_csv(
    #     'Z:/pilot_for_heatmapAllSens.csv')

    print(data.head(5))
    labels = data[data.columns[2]].astype(str)  # retrieves gene names
    data = data[data.columns[9:]].astype(float)
    f, ax = plt.subplots(figsize=(12, 12))
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns.set(style="white", font_scale=0.2)
    sns.heatmap(data, mask=False, cmap="RdBu", robust=False, vmax=3, center=0,
                square=False, cbar_kws={"shrink": .5}, xticklabels='auto',
                yticklabels=labels)  # cmap="PRGn" (Green), cmap="RdBu" (Blue), cmap=cmap (Red), linewidths=.5
    # sns.set(font_scale=0.4)
    plt.show()

def show_manual_heatmap_by_gene():
    '''
    plots a heatmap readout of the defects manually scored from csv
    IN PROGRESS
    :return:
    '''
    import pandas as pd
    import seaborn as sns
    from string import ascii_letters
    import matplotlib.pyplot as plt
    from fileLookup import FILENAME_WEB_CONTENT


    file = FOLDER_IN + 'Gene_lists/manualScoringForSeabornB.csv'

    data = pd.read_csv(file)
    # print(data.head(5))

    labels = data[data.columns[0]].astype(str)  # retrieves embd number
    data = data[data.columns[1:35]].astype(float)  # skips first column (labels)
    for index, row in data.iterrows():
        print(index)
        lab = labels[index]  # retrieves EMBD number
        label = ["${}(RNAi)$".format(getGeneCommonName(lab))] # converts to genename and makes label italics
        row = pd.DataFrame(row)  # converts series back to dataframe so seaborn can properly see the 1D data

        fig = plt.figure(figsize=(5, 15))  # generates matplotlib fig
        ax = plt.subplot(1,2,2)  # plots to the second subplot (did this to avoid text cutoff during saving step)
        # f, ax = plt.subplots(figsize=(0.1, 10))
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        sns.set(style="dark", font_scale=0.8)

        sns.heatmap(row, mask=False, cmap="RdBu", robust=False, vmax=3, center=0,
                    square=True, cbar=False, cbar_kws={"shrink": 0.2,}, xticklabels= label,
                    yticklabels= data.columns, linewidths=0.1, linecolor="dimgrey" )  # cmap="PRGn" (Green), cmap="RdBu" (Blue), cmap=cmap (Red), linewidths=.5
        ax.tick_params(left=False, bottom=False, labelbottom=False, labeltop=True, labelsize=12)
        # plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
        fig.tight_layout()  # optimizes figure layout before saving
        fig.savefig(FILENAME_WEB_CONTENT + 'EMBD{n:04}/EMBD{n:04}_manual.svg'.format(n=int(lab), bbox_inches="tight", pad_inches=1))  # bbox_inches="tight", pad_inches=0.5
        # plt.show()
        plt.close()

def show_manual_heatmap_by_gene_germline():
    '''
    plots a heatmap readout of the defects manually scored from csv
    IN PROGRESS
    :return:
    '''
    import pandas as pd
    import seaborn as sns
    from string import ascii_letters
    import matplotlib.pyplot as plt
    from fileLookup import FILENAME_WEB_CONTENT_GERMLINE
    from GermlineClass import GermlineClass
    import os

    file = FOLDER_IN + 'Gene_lists/manual_germline_for_seaborn.csv'
    saveFlag = True

    data = pd.read_csv(file)
    # print(data.head(5))
    labels = data[data.columns[0]].astype(str)  # retrieves embd number
    data = data[data.columns[1:95]].astype(float)  # skips first column (labels)
    for index, row in data.iterrows():
        print(index)
        int_lab = int(float(labels[index]))  # retrieves SS number as integer
        g = GermlineClass(int_lab)
        # lab = [g.get_genename_from_germ_id(int_lab)]
        lab = ['${}(RNAi)$'.format(g.get_genename_from_germ_id(int_lab))]
        # label = ["${}(RNAi)$".format(getGeneCommonName(lab))] # converts to genename and makes label italics
        row = pd.DataFrame(row)  # converts series back to dataframe so seaborn can properly see the 1D data

        fig = plt.figure(figsize=(8, 20))  # generates matplotlib fig
        ax = plt.subplot(1,2,2)  # plots to the second subplot (did this to avoid text cutoff during saving step)
        # f, ax = plt.subplots(figsize=(0.1, 10))
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        sns.set(style="dark", font_scale=0.5)

        sns.heatmap(row, mask=False, cmap="RdBu", robust=False, vmax=1, center=0,
                    square=True, cbar=False, cbar_kws={"shrink": 0.2,}, xticklabels= lab,
                    yticklabels= data.columns, linewidths=0.1, linecolor="dimgrey" )  # cmap="PRGn" (Green), cmap="RdBu" (Blue), cmap=cmap (Red), linewidths=.5
        ax.tick_params(left=False, bottom=False, labelbottom=False, labeltop=True, labelsize=12)
        plt.xticks(rotation=0)
        # plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
        # plt.tight_layout()  # optimizes figure layout before saving
        if saveFlag:
            folder = FILENAME_WEB_CONTENT_GERMLINE + 'SS{n:03}'.format(n=int_lab)
            if not os.path.exists(folder):
                os.makedirs(folder)
            plt.savefig(FILENAME_WEB_CONTENT_GERMLINE + 'SS{n:03}/SS{n:03}_manual.svg'.format(n=int_lab),
                        bbox_inches="tight", pad_inches=1)  # bbox_inches="tight", pad_inches=0.5
        # fig.savefig(FILENAME_WEB_CONTENT_GERMLINE + 'SS{n:03}/SS{n:03}_manual.svg'.format(n=int_lab), bbox_inches="tight", pad_inches=1)  # bbox_inches="tight", pad_inches=0.5
        # plt.show()
        plt.close()

def plot_tIntVsPtArrest():
    '''
    Contains tInt data comparison with manually scored point of arrest (see mySQL_man_auto_comparison020719.xls)
    :return:
    '''
    cv1_label = 'Green'
    curve1 = [70320054.77, 212760998.2, 304947413.2, 385072889.8, 373130221.6, 389941444.6, 426171857.2]
    std1 = [139071541.3, 203854968.3, 167659088.9, 111694381.3, 122814475.5, 77271653.01, 112482686.1]
    cv2_label = 'Red'
    curve2 = [166867576.6, 337398232.3, 498832440, 535070746.7, 586401961.5, 583489611.2, 603453252.9]
    std2 = [195261491.7, 196334983.4, 181263849.4, 124229347.1, 127678347.9, 100582800.8, 115528329.6]
    x_label = 'Point of arrest'
    # x_points= ['no markers','early','comma','1.5 fold','2 fold','3 fold', 'control' ]
    x_points = [0, 1, 2, 3, 4, 5, 6]
    n = [84, 150, 179, 46, 72, 231, 447]
    return (x_label, x_points, cv1_label, curve1, std1, cv2_label, curve2, std2)


def plot_maxIntGLSVsPtArrest():
    '''
    Contains max Int data GLS comparison with manually scored point of arrest (see mySQL_man_auto_comparison020719.xls)
    :return:
    '''
    cv1_label = 'Max Int Green'
    curve1 = [466118896.8, 921905991.1, 1060149594, 1172954259, 1196549050, 1239860549, 1349671774]
    std1 = [364250316.2, 414110712.7, 202625254.4, 212474503.4, 149360301.5, 188736449.5, 184547703.4]
    cv2_label = 'Max Int Red'
    curve2 = [530685396, 979253209.7, 1294025414, 1371486514, 1415590433, 1545887783, 1678901654]
    std2 = [422512273.2, 445843009.2, 228715909, 186119848.5, 162804850.1, 188734228.1, 180873951.8]
    x_label = 'Point of arrest'
    # x_points= ['no markers','early','comma','1.5 fold','2 fold','3 fold', 'control' ]
    x_points = [0, 1, 2, 3, 4, 5, 6]
    n = [50, 235, 330, 115, 115, 360, 413]  # number of embryos
    return (x_label, x_points, cv1_label, curve1, std1, cv2_label, curve2, std2)


def sigmoidGLS(x, a, m, s, r):
    return a * (1. - (1. + np.exp((x - m) / s)) ** (-r))


def sigmoidMS(x, a, b, m, s, r):
    # return a - b * (1. + np.exp((x - m) / s)) ** (-r)
    return a - a * (1. + np.exp((x - m) / s)) ** (-r)


# fig.plot(x, sigmoidal(x / self.tScale, a, b, m / self.tScale, s, r) * self.scale[curveName], color='k')
# fig.plot(((self.tMove - self.t0) * self.tScale, (self.tMove - self.t0) * self.tScale), (0, a),
#          color='b')

def make_CoM_tif_series(embd, embryo, strain='MS', blackback=False):
    '''
    Takes in EMBD# (i.e. 0403), embryo number (i.e. 2) and strain and generates a fig with two subplots. The first is a
    max projection with CoM dots that overlay the position of tissue within the image (2) a graph that shows the control
     traces in grey and adds in the data for the specified embryo by timepoint to create a real-time tif series that can
     be made into a movie. Shows (if ShowFlag) and saves (if saveFlag) tif series. Note if GLS is selected, it shows the center
     of mass, not the distance between the two tissues' center of mass (as measured)-- this is less useful to display.
    :param embd: embd number
    :param embryo: embryo number
    :param strain: 'GLS' or 'MS'
    :param blackback: allows you to generate the tifs with a black background for presentation purposes
    :return:shows and/or saves an image or series of images +plots for CoM data for either GLS or MS
    '''

    from Embryos import getCenterOfMass, MSEmbryo
    import matplotlib

    folders = [FOLDER_IN + 'cropped/EMBD{n:04}/MS/Emb{x}/'.format(n=embd, x=embryo)]
    for folder in folders[:]:
        emb = MSEmbryo(folder, check_version=False)

    saveFlag = True
    showFlag = False
    # subplots = ['max project', 'CoM0G', 'CoM1G','CoM0R','CoM1R']
    rows = 2
    columns = 1

    if blackback:
        text_col = 'white'
    else:
        text_col = "black"

    '''this gives access to control ave curves'''
    r = RNAiClass(embd)  # this is a work around- just need an embryo object to call to access the control avg curves
    r.getEmbryos()
    r.asignUse()
    c = ControlClass()

    if strain == 'MS':
        folders = r.MSUseEmbs  # retrieves a list of embryo objects that are used (considered by automated analysis)
    elif strain == 'GLS':
        folders = r.GLSUseEmbs  # retrieves a list of embryo objects that are used (considered by automated analysis)

    for folder in folders:
        emb_obj = folder

    '''this accesses the embryo'''
    if strain == 'MS':
        folders = [FOLDER_IN + 'cropped/EMBD{n:04}/MS/Emb{x}/'.format(n=embd, x=embryo)]
        for folder in folders[:]:
            emb = MSEmbryo(folder, check_version=False)
    elif strain == 'GLS':
        folders = [FOLDER_IN + 'cropped/EMBD{n:04}/GLS/Emb{x}/'.format(n=embd, x=embryo)]
        for folder in folders[:]:
            emb = GSEmbryo(folder, check_version=False)

    for i in range(1, 32):  # steps through the timepoints to plot
        '''generate the figure'''
        fig = plt.figure(figsize=(5, 8))  # creates figure
        if blackback:
            plt.rcParams.update({
                "lines.color": "white",
                "patch.edgecolor": "white",
                "text.color": "black",
                "axes.facecolor": "white",
                "axes.edgecolor": "lightgray",
                "axes.labelcolor": "white",
                "xtick.color": "white",
                "ytick.color": "white",
                "grid.color": "lightgray",
                "figure.facecolor": "black",
                "figure.edgecolor": "black",
                "savefig.facecolor": "black",
                "savefig.edgecolor": "black"})
        name = getGeneCommonName(embd)
        fig.suptitle("{g}(RNAi) embryo{e}\n EMBD{n:04} {s}".format(g=name, n=embd, e=embryo, s=strain), fontsize=14,
                     style='italic', y=.97, color=text_col)  # adds title to the top of all subplots

        if strain == 'MS':
            time_thres = 500

            '''plot max project with CoM dots in the first subplot'''
            plt.subplot(rows, columns, 1)  # defines the subplot to draw on
            plt.title('Maximum Projection', fontsize=12, y=1.0, color=text_col)  # makes subplot title
            plt.markerSize = 300
            t = emb.tScale * (i - emb.t0)  # this scales time to be aligned with plots

            xe, ye, ze = emb.getEmbCenter()
            plt.scatter([xe], [ye], color='gray',edgecolors= 'white', s=150) # plot embryo center

            im, (xG, yG), (xR, yR) = emb.getMaxProj(t)

            'testing out a fix below'
            # xR = emb.CMPos[1, :, 0][i]
            # yR = emb.CMPos[1, :, 1][i]
            # xG = emb.CMPos[0, :, 0][i]
            # yG = emb.CMPos[0, :, 1][i]

            # if bi:
            #     im[0] = im[0] > 0
            #     im[1] = im[1] > 0
            plt.imshow(im)
            plt.scatter([xG], [yG], color='darkseagreen', edgecolors= 'white', s=150)
            plt.scatter([xR], [yR], color='salmon', edgecolors= 'white', s=150)
            ax = plt.subplot(rows, columns, 1)
            ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.axis('off')
            ax.set_axis_bgcolor('k')



            '''plot CoM values in the second subplot'''
            plt.subplot(rows, columns, 2)  # defines the subplot to draw on
            plt.subplots_adjust(hspace=0.1, bottom=0.1)
            plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
            ax = plt.subplot(rows, columns, 2)
            ax.tick_params(axis='x', bottom=True, top=False, left=False, right=False, labelsize=8)
            ax.tick_params(axis='y', bottom=False, top=False, left=True, right=False, labelsize=8)
            xtick_labels = ['-200', '-100', '0', '100', '200', '300', '400', '500', '600', '700']
            ax.yaxis.get_offset_text().set_fontsize(8)

            plt.title('Distance to Center of Mass- D Long, t ={0}'.format(t), fontsize=12, y=1.0,
                      color=text_col)  # makes subplot title

            # plot red curve
            x, y = emb.getCurve('CoM1R')
            # print(x,y)
            x = 20 * x  # adjusts x axis to time rather than frame
            x = x[np.where(x<time_thres)]
            y = y[np.where(x<time_thres)]
            plt.plot(x[0:i + 1], y[0:i + 1], 'r', label='Neurons')  # plots data in red

            # plot green curve
            x1, y1 = emb.getCurve('CoM1G')
            x1 = 20 * x1  # adjusts x axis to time rather than frame
            x1 = x1[np.where(x1<time_thres)]
            y1 = y1[np.where(x1<time_thres)]
            plt.plot(x1[0:i + 1], y1[0:i + 1], 'g',
                     label='Epidermis')  # plots data in green

            '''plot control aves on second subplot'''
            if 'CoM1R' in emb_obj.avgCurves:
                xR, yR, yerrR = emb_obj.avgCurves['CoM1R']
                xR = xR * 20
                xR = xR[np.where(xR<time_thres)]
                yR = yR[np.where(xR<time_thres)]
                yerrR = yerrR[np.where(xR<time_thres)]
                upper = np.add(yR, yerrR)
                lower = np.subtract(yR, yerrR)
                plt.plot(xR, yR, 'darksalmon')
                plt.fill_between(xR, upper, lower, color='0.8', alpha=0.5)  # plots with shaded error bars

            if 'CoM1G' in emb_obj.avgCurves:
                xG, yG, yerrG = emb_obj.avgCurves['CoM1G']
                xG = xG * 20
                xG = xG[np.where(xG<time_thres)]
                yG = yG[np.where(xG<time_thres)]
                yerrG = yerrG[np.where(xG<time_thres)]
                upper = np.add(yG, yerrG)
                lower = np.subtract(yG, yerrG)
                plt.plot(xG, yG, 'yellowgreen')
                plt.fill_between(xG, upper, lower, color='0.8', alpha=0.5)  # plots with shaded error bars
                plt.subplots_adjust(wspace=0.35, hspace=0.5)

        elif strain == 'GLS':
            '''plot max project with CoM dots in the first subplot'''

            plt.subplot(rows, columns, 1)  # defines the subplot to draw on
            plt.title('Maximum Projection', fontsize=12, y=1.0, color=text_col)  # makes subplot title
            plt.markerSize = 300
            t = emb.tScale * (i - emb.t0)  # this scales time to be aligned with plots
            im, (xG, yG), (xR, yR) = emb.getMaxProj(t)
            # if bi:
            #     im[0] = im[0] > 0
            #     im[1] = im[1] > 0
            plt.imshow(im)

            emb.loadAllSpots()
            rounded_t= np.round(i)
            if emb.getSpotsPositions(0,rounded_t) is not None:
                xG, yG, zG = getCenterOfMass(emb.getSpotsPositions(0, rounded_t))
                plt.scatter([xG], [yG], color='darkseagreen', edgecolors= "white", s=150)
            if emb.getSpotsPositions(1, rounded_t) is not None:
                xR, yR, zR = getCenterOfMass(emb.getSpotsPositions(1,rounded_t))
                plt.scatter([xR], [yR], color='salmon', edgecolors= "white", s=150)
            if emb.getSpotsPositions(2, rounded_t) is not None:
                xY, yY, zY = getCenterOfMass(emb.getSpotsPositions(2, rounded_t))
                plt.scatter([xY], [yY], color='gold', edgecolors= "white", s=150)
            # plt.scatter([xG], [yG], color='darkseagreen', s=150)
            # plt.scatter([xR], [yR], color='salmon', s=150)
            # plt.scatter([xY], [yY], color='gold', s=150)

            ax = plt.subplot(rows, columns, 1)
            ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.axis('off')
            ax.set_axis_bgcolor('k')

            '''plot CoM values in the second subplot'''
            plt.subplot(rows, columns, 2)  # defines the subplot to draw on
            plt.subplots_adjust(hspace=0.1, bottom=0.1)
            plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
            ax = plt.subplot(rows, columns, 2)
            ax.tick_params(axis='x', bottom=True, top=False, left=False, right=False, labelsize=8)
            ax.tick_params(axis='y', bottom=False, top=False, left=True, right=False, labelsize=8)
            xtick_labels = ['-200', '-100', '0', '100', '200', '300', '400', '500', '600', '700']
            ax.yaxis.get_offset_text().set_fontsize(8)
            plt.title('Distance Between Centers of Mass- D Long', fontsize=12, y=1.0,
                      color=text_col)  # makes subplot title
            time_thres = 500 #set time cutoff for plotting
            # plot red curve
            x, y = emb.getCurve('CoM1R')
            x = 20 * x  # adjusts x axis to time rather than frame
            x = x[np.where(x<time_thres)]
            y = y[np.where(x<time_thres)]
            plt.plot(x[0:i + 1], y[0:i + 1], 'r',
                     label='Ectoderm to Endoderm')  # plots data in red

            # plot green curve
            x1, y1 = emb.getCurve('CoM1Y')
            x1 = 20 * x1  # adjusts x axis to time rather than frame
            x1 = x1[np.where(x1<time_thres)]
            y1 = y1[np.where(x1<time_thres)]
            plt.plot(x1[0:i + 1], y1[0:i + 1], 'gold',
                     label='Mesoderm to Endoderm')  # plots data in yellow

            '''plot control aves on second subplot'''
            if 'CoM1R' in emb_obj.avgCurves:
                xR, yR, yerrR = emb_obj.avgCurves['CoM1R']
                xR = xR * 20
                xR = xR[np.where(xR<time_thres)] #limits the plot to the time range specified in time_thres
                yR = yR[np.where(xR<time_thres)]
                yerrR = yerrR[np.where(xR<time_thres)]
                upper = np.add(yR, yerrR)
                lower = np.subtract(yR, yerrR)
                plt.plot(xR, yR, 'darksalmon')
                plt.fill_between(xR, upper, lower, color='0.8', alpha=0.5)  # plots with shaded error bars

            if 'CoM1Y' in emb_obj.avgCurves:
                xY, yY, yerrY = emb_obj.avgCurves['CoM1Y']
                xY = xY * 20
                xY = xY[np.where(xY<time_thres)] #limits the plot to the time range specified in time_thres
                yY = yY[np.where(xY<time_thres)]
                yerrY = yerrY[np.where(xY<time_thres)]
                upper = np.add(yY, yerrY)
                lower = np.subtract(yY, yerrY)
                plt.plot(xY, yY, 'khaki')
                plt.fill_between(xY, upper, lower, color='0.8', alpha=0.5)  # plots with shaded error bars
                plt.subplots_adjust(wspace=0.35, hspace=0.5)
        plt.legend(loc=1, prop={'size': 8}, frameon=False)  # loc=1 is upper right, 2=upper left

        if saveFlag:
            if strain == 'GLS':
                plt.savefig(
                    'Z:/EMBD{n:04}_{s}_{t}.tif'.format(n=embd,
                                                                                                                s=strain,
                                                                                                                t=i))
            elif strain == 'MS':
                plt.savefig(
                    'Z:/EMBD{n:04}_{s}_{t}.tif'.format(n=embd,
                                                                                                               s=strain,
                                                                                                               t=i))

        if showFlag:
            plt.show()


def show_head_endon(embd, embryo, strain='MS', blackback=False):
    '''
    Takes in EMBD# (i.e. 0403), embryo number (i.e. 2) and generates a fig with two subplots. The first is an end-on
    max projection of the head region (2) a graph of the number of pixels above threshold from projected on the end-on view of 1/3 of the embryo (normalized
    by the area of the embryo). Shows (if ShowFlag) and saves (if saveFlag) tif series. Can save with a black background for presentation purposes
    :param embd: embd number
    :param embryo: embryo number
    :param strain: 'MS'
    :param blackback: if True, generates a black background version of the figure
    :return:shows and/or saves an image or series of images +plots for CoM data for MS
    '''
    from scipy.ndimage import zoom
    from myFigure import *
    from cv2 import imread

    saveFlag = True
    showFlag = False
    no_yscale = False
    rows = 2
    columns = 1

    if blackback:
        text_col = 'white'
    else:
        text_col = "black"

    '''this gives access to control ave curves'''
    r = RNAiClass(embd)  # this is a work around- just need an embryo object to call to access the control avg curves
    r.getEmbryos()
    r.asignUse()
    # c = ControlClass()

    folders = r.MSUseEmbs  # retrieves a list of embryo objects that are used (considered by automated analysis)
    for folder in folders:
        emb_obj = folder

    '''this accesses the embryo'''
    if strain == 'MS':
        folders = [FOLDER_IN + 'cropped/EMBD{n:04}/MS/Emb{x}/'.format(n=embd, x=embryo)]
        for folder in folders[:]:
            emb = MSEmbryo(folder, check_version=False)

    if emb.image is None:
        emb.loadImages()

    '''plot max project of head region in binary mask'''
    # for t in emb.time:
    for t in emb.time:
        ims = emb.image.images[int(t), :,
              0]  # retrieves all z planes for timepoint t and channel 0, shape is (18,122,214)
        length = ims[0].shape[1]
        ims = ims[:, :, :int(length / 3.)]  # redefines ims as all z, all y, and 1/3 of x
        ims = (ims > emb.thresh[0])
        # imProj = np.max(ims, axis=2)
        # imProj = (imProj > emb.thresh[0])

        mainfig = plt.figure(figsize=(5, 8))  # creates figure
        if blackback:
            plt.rcParams.update({
                "lines.color": "white",
                "patch.edgecolor": "white",
                "text.color": "black",
                "axes.facecolor": "white",
                "axes.edgecolor": "lightgray",
                "axes.labelcolor": "white",
                "xtick.color": "white",
                "ytick.color": "white",
                "grid.color": "lightgray",
                "figure.facecolor": "black",
                "figure.edgecolor": "black",
                "savefig.facecolor": "black",
                "savefig.edgecolor": "black"})
        name = getGeneCommonName(embd)
        mainfig.suptitle("{g}(RNAi) embryo{e}\n EMBD{n:04} {s}".format(g=name, n=embd, e=embryo, s=strain), fontsize=14,
                         style='italic', y=.97, color=text_col)  # adds title to the top of all subplots
        plt.subplot(2, 1, 1)  # defines the first subplot to draw on
        ax = plt.subplot(2, 1, 1)
        ax.set_axis_bgcolor('k')
        ax.spines['bottom'].set_color('white')
        ax.spines['left'].set_color('white')
        ax.tick_params(axis='x', colors='white')
        ax.tick_params(axis='y', colors='white')
        ax.yaxis.label.set_color('white')
        ax.xaxis.label.set_color('white')
        ax.set_aspect('equal')
        ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
        ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())
        ax.axis('off')
        # fig1 = myFigure()  # makes the side view projection figure
        # fig1.makeBlack()
        # fig1.set_axes_equal()
        # fig1.noAxis()

        imProj = np.max(ims,
                        axis=2)  # makes a projection on x-axis to give a z-y display (shape is 18,2,128, which is Z, channels, y)
        imProjG = zoom(imProj, (emb.dz, 1), order=0)
        plt.imshow(imProjG)
        ax.imshow(imProjG, cmap='Greys_r')
        plt.legend(loc=None)

        i = int(t)  # lazy fix-- steps through the timepoints to plot
        '''plot CoM values in the second subplot'''
        plt.subplot(2, 1, 0)  # defines the subplot to draw on
        plt.subplots_adjust(hspace=0.1, bottom=0.1)
        plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
        ax = plt.subplot(2, 1, 0)
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        xtick_labels = ['-200', '-100', '0', '100', '200', '300', '400', '500', '600', '700']
        ax.yaxis.get_offset_text().set_fontsize(8)

        plt.title('Head Intensity', fontsize=12, y=1.0, color=text_col)  # makes subplot title

        # plot green curve
        x1, y1 = emb.getCurve('headInt')
        x1 = 20 * x1  # adjusts x axis to time rather than frame
        if no_yscale:
            y1= y1/emb.scale['headInt']
        inds = np.invert(np.isnan(y1))
        first_true = np.where(inds == True)[0][0]
        # plt.plot(x1[first_true:i], y1[first_true:i], 'g', label='Head Intensity (headInt)')  # plots data in green

        if inds[i]:
            plt.plot(x1[first_true:i + 1], y1[first_true:i + 1], 'g',
                     label='Head Intensity (headInt)')  # plots data in green

        # plt.plot(x1[i], y1[i], 'g', label='Head Intensity (headInt)')  # plots data in green
        #     print('True')
        # else:
        #     print('false')

        # plt.plot(x1[0:i], y1[0:i], 'g', label='Head Intensity (headInt)')  # plots data in green

        '''plot control aves on second subplot'''

        if 'headInt' in emb_obj.avgCurves:
            xG, yG, yerrG = emb_obj.avgCurves['headInt']
            xG = xG * 20
            upper = np.add(yG, yerrG)
            lower = np.subtract(yG, yerrG)
            plt.plot(xG, yG, 'yellowgreen')
            plt.fill_between(xG, upper, lower, color='0.8', alpha=0.5)  # plots with shaded error bars
            plt.subplots_adjust(wspace=0.35, hspace=0.5)
        if showFlag:
            plt.show()
            mainfig.show()
        if saveFlag:
            mainfig.savefig(
                'Z:/EMBD{n:04}_{s}_{t}.tif'.format(
                    n=embd, s=strain,
                    t=i))


def make_CoM_tif_series_endon(embd, embryo, strain='MS', blackback=False):
    '''
    Takes in EMBD# (i.e. 0403), embryo number (i.e. 2) and generates a fig with two subplots. The first is an end-on
    max projection with CoM dots that overlay the position of tissue within the image (2) a graph that shows the control
     traces in grey and adds in the data for the specified embryo by timepoint to create a real-time tif series that can
     be made into a movie. Shows (if ShowFlag) and saves (if saveFlag) tif series. Not config to show GLS  (it shows the center
     of mass, not the distance between the two tissues' center of mass (as measured)-- this is less useful to display).
    :param embd: embd number
    :param embryo: embryo number
    :param strain: 'MS'
    :param blackback: allows you to make a figure with a black background for presentation purposes
    :return:shows and/or saves an image or series of images +plots for CoM data for MS
    '''
    from scipy.ndimage import zoom
    from myFigure import *
    from cv2 import imread
    from Embryos import getCenterOfMass

    saveFlag = True
    showFlag = False
    rows = 2
    columns = 1

    if blackback:
        text_col = 'white'
    else:
        text_col = "black"

    '''this gives access to control ave curves'''
    r = RNAiClass(embd)  # this is a work around- just need an embryo object to call to access the control avg curves
    r.getEmbryos()
    r.asignUse()
    # c = ControlClass()

    if strain == 'MS':
        folders = r.MSUseEmbs  # retrieves a list of embryo objects that are used (considered by automated analysis)
    elif strain == 'GLS':
        folders = r.GLSUseEmbs  # retrieves a list of embryo objects that are used (considered by automated analysis)
    for folder in folders:
        emb_obj = folder

    '''this accesses the embryo'''
    if strain == 'MS':
        folders = [FOLDER_IN + 'cropped/EMBD{n:04}/MS/Emb{x}/'.format(n=embd, x=embryo)]
        for folder in folders[:]:
            emb = MSEmbryo(folder, check_version=False)
    elif strain == 'GLS':
        folders = [FOLDER_IN + 'cropped/EMBD{n:04}/GLS/Emb{x}/'.format(n=embd, x=embryo)]
        for folder in folders[:]:
            emb = GSEmbryo(folder, check_version=False)


    if emb.image is None:
        emb.loadImages()

    '''plot max project with CoM dots in the first subplot'''

    # for t in emb.time:
    for i in range(0, 31):
        t = emb.tScale * (i - emb.t0)
        if int(np.round(t)) >= 0:  # waits until R, distR = emb.getCurve('CoM0R') at t 0 before starting to plot
            mainfig = plt.figure(figsize=(5, 8))  # creates figure
            if blackback:
                plt.rcParams.update({
                    "lines.color": "white",
                    "patch.edgecolor": "white",
                    "text.color": "black",
                    "axes.facecolor": "white",
                    "axes.edgecolor": "lightgray",
                    "axes.labelcolor": "white",
                    "xtick.color": "white",
                    "ytick.color": "white",
                    "grid.color": "lightgray",
                    "figure.facecolor": "black",
                    "figure.edgecolor": "black",
                    "savefig.facecolor": "black",
                    "savefig.edgecolor": "black"})
            name = getGeneCommonName(embd)
            mainfig.suptitle("{g}(RNAi) embryo{e}\n EMBD{n:04} {s}".format(g=name, n=embd, e=embryo, s=strain), fontsize=14,
                             style='italic', y=.97, color=text_col)  # adds title to the top of all subplots
            plt.subplot(2, 1, 1)  # defines the first subplot to draw on
            ax = plt.subplot(2, 1, 1)
            ax.set_axis_bgcolor('k')
            ax.spines['bottom'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.tick_params(axis='x', colors='white')
            ax.tick_params(axis='y', colors='white')
            ax.yaxis.label.set_color('white')
            ax.xaxis.label.set_color('white')
            ax.set_aspect('equal')
            ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.axis('off')
            plt.title('Maximum Projection', fontsize=12, y=1.0, color=text_col)  # makes subplot title


            if emb.tScale>1.2: tSc = 1.2
            elif emb.tScale <0.8: tSc = 0.8
            else: tSc = emb.tScale

            slide = int(np.round(i-2))  # frame to show or plot

            # slide = int(np.round(i / tSc))  # frame to show or plot
            if slide >30: slide = 30
            # slide= int(np.round(i))

            ims = emb.image.images[slide, :,
                  0:2]  # retrieves all z planes for frame i and channel 0, shape is (18,2,128,206)... this is (frame, z, channels,y,x)


            imProj = np.max(ims,
                            axis=3)  # makes a projection on x-axis to give a z-y display (shape is 18,2,128, which is Z, channels, y)
            imProjG = zoom(imProj[:, 0, :], (emb.dz, 1), order=0)
            imProjR = zoom(imProj[:, 1, :], (emb.dz, 1), order=0)
            imProjG = (imProjG / 256).astype('uint8')  # converts to 8-bit
            imProjR = (imProjR / 256).astype('uint8')
            rg_image = np.stack([imProjR, imProjG, np.zeros_like(imProjR)],
                                axis=2)  # makes an R,G,B stack from the grayscale images

            plt.imshow(rg_image)
            plt.legend(loc=None)


            # i = int(np.round(t))  # lazy fix-- steps through the timepoints to plot

            if strain == "MS":
                #note that control emb 400_1, seems to need to be shifted by 2. i.e. slide = slide-2. I do not understand why
                time_thres = 500

                xe, ye, ze = emb.getEmbCenter()

                ze = 144-ze # 144 is the total Z range. Subtracting ze here is plotting from bottom up rather than top down
                ye = 128-ye
                plt.scatter(ye, ze, color='gray', edgecolors='white', s=150)  # plot embryo center

                # xR, yR, zR = emb.image.getIntCenter(slide, 1, 0)  # note that this puts the spot in closer to center of embryo than emb.CMPos
                # xG, yG, zG = emb.image.getIntCenter(slide, 0, 0)
#
                'FIXME- testing out a fix below-- note that getIntCenter and CMPos should be the same, but they are not. CMPos values seem more similar to expectations ' \
                'than getIntCenter values unsure which to use'
                G2_z = emb.CMPos[0, :, 2]
                G2_y = emb.CMPos[0, :, 1]
                R2_z = emb.CMPos[1, :, 2]
                R2_y = emb.CMPos[1, :, 1]
                zR = R2_z[i]
                yR = R2_y[i]
                zG = G2_z[i]
                yG = G2_y[i]



#               # fig1 = myFigure()  # makes the side view projection figure
                # fig1.makeBlack()
                # fig1.set_axes_equal()
                # fig1.noAxis()

                # plt.scatter(yG, zG * emb.dz, color='darkseagreen', edgecolors='white', s=150)
                # plt.scatter(yR, zR * emb.dz, color='salmon', edgecolors='white', s=150)

                plt.scatter(yG, zG, color='darkseagreen', edgecolors='white', s=150)
                plt.scatter(yR, zR, color='salmon', edgecolors='white', s=150)

                '''plot CoM values in the second subplot'''
                plt.subplot(2, 1, 0)  # defines the subplot to draw on
                plt.subplots_adjust(hspace=0.1, bottom=0.1)
                plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
                ax = plt.subplot(2, 1, 0)
                ax.tick_params(axis='x', bottom=True, top=False, left=False, right=False, labelsize=8)
                ax.tick_params(axis='y', bottom=False, top=False, left=True, right=False, labelsize=8)
                xtick_labels = ['-200', '-100', '0', '100', '200', '300', '400', '500', '600', '700']
                ax.yaxis.get_offset_text().set_fontsize(8)
                plt.xlabel('time(mins)', fontsize=10)  # labels x-axis

                plt.title('Distance to Center of Mass- D Short', fontsize=12, y=1.0, color=text_col)  # makes subplot title
                # plot red curve
                x, y = emb.getCurve('CoM0R')  # already time shifted and tScale adjusted, i think
                x = 20 * (x)  # adjusts x axis to time rather than frame
                x = x[np.where(x < time_thres)]
                y = y[np.where(x < time_thres)]


                # plt.plot(x[0:i], y[0:i], 'r', label='Neuronal Position- Short Axis (Com0R)')  # plots data in red
                plt.plot(x[0:np.round(i)], y[0:np.round(i)], 'r', label='Neurons')  # plots data in red

                # plot green curve
                x1, y1 = emb.getCurve('CoM0G')
                x1 = 20 * (x1)  # adjusts x axis to time rather than frame

                x1 = x1[np.where(x1 < time_thres)] #limits plot to time threshold
                y1 = y1[np.where(x1 < time_thres)]

                # plt.plot(x1[0:i], y1[0:i], 'g', label='Epidermal Position- Short Axis (Com0G)')  # plots data in green
                plt.plot(x1[0:np.round(i)], y1[0:np.round(i)], 'g', label='Epidermis')  # plots data in green

                '''plot control aves on second subplot'''
                if 'CoM0R' in emb_obj.avgCurves:
                    xR, yR, yerrR = emb_obj.avgCurves['CoM0R']
                    xR = xR * 20
                    xR = xR[np.where(xR < time_thres)]
                    yR = yR[np.where(xR < time_thres)]
                    yerrR = yerrR[np.where(xR < time_thres)]

                    upper = np.add(yR, yerrR)
                    lower = np.subtract(yR, yerrR)
                    plt.plot(xR, yR, 'darksalmon')
                    plt.fill_between(xR, upper, lower, color='0.8', alpha=0.5)  # plots with shaded error bars

                if 'CoM0G' in emb_obj.avgCurves:
                    xG, yG, yerrG = emb_obj.avgCurves['CoM0G']
                    xG = xG * 20
                    xG = xG[np.where(xG < time_thres)]
                    yG = yG[np.where(xG < time_thres)]
                    yerrG = yerrR[np.where(xG < time_thres)]

                    upper = np.add(yG, yerrG)
                    lower = np.subtract(yG, yerrG)
                    plt.plot(xG, yG, 'yellowgreen')
                    plt.fill_between(xG, upper, lower, color='0.8', alpha=0.5)  # plots with shaded error bars
                    plt.subplots_adjust(wspace=0.35, hspace=0.5)
                    plt.legend(loc=1, prop={'size': 8}, frameon=False) # loc=1 is upper right, 2=upper left

                if showFlag:
                    plt.show()
                    mainfig.show()
                if saveFlag:
                    mainfig.savefig(
                        'Z:/EMBD{n:04}_{s}_{t}.tif'.format(n=embd,
                                                                                                                         s=strain, t=i))

            '''in progress below!'''
            if strain == "GLS":
                emb.loadAllSpots()
                rounded_t = np.round(i)
                if emb.getSpotsPositions(0, rounded_t) is not None:
                    xG, yG, zG = getCenterOfMass(emb.getSpotsPositions(0, rounded_t))
                    plt.scatter([yG], [zG], color='darkseagreen', edgecolors="white", s=150)
                if emb.getSpotsPositions(1, rounded_t) is not None:
                    xR, yR, zR = getCenterOfMass(emb.getSpotsPositions(1, rounded_t))
                    plt.scatter([yR], [zR], color='salmon', edgecolors="white", s=150)
                if emb.getSpotsPositions(2, rounded_t) is not None:
                    xY, yY, zY = getCenterOfMass(emb.getSpotsPositions(2, rounded_t))
                    plt.scatter([yY], [zY], color='gold', edgecolors="white", s=150)

                '''plot CoM values in the second subplot'''
                plt.subplot(2, 1, 0)  # defines the subplot to draw on
                plt.subplots_adjust(hspace=0.1, bottom=0.1)
                plt.xlabel('time(mins)', fontsize=10)  # labels x-axis
                ax = plt.subplot(rows, columns, 2)
                ax.tick_params(axis='x', bottom=True, top=False, left=False, right=False, labelsize=8)
                ax.tick_params(axis='y', bottom=False, top=False, left=True, right=False, labelsize=8)
                xtick_labels = ['-200', '-100', '0', '100', '200', '300', '400', '500']
                # xtick_labels = ['-200', '-100', '0', '100', '200', '300', '400', '500', '600', '700']

                ax.yaxis.get_offset_text().set_fontsize(8)
                plt.title('Distance Between Centers of Mass- D Short', fontsize=12, y=1.0,
                          color=text_col)  # makes subplot title
                plt.xlabel('time(mins)', fontsize=10)  # labels x-axis

                # plot red curve
                time_thres = 500  # allows you to truncate the plot to not include so much movement
                x, y = emb.getCurve('CoM0R')
                x = 20 * (x)  # adjusts x axis to time rather than frame
                x = x[np.where(x<time_thres)]
                y = y[np.where(x<time_thres)]
                plt.plot(x[0:i + 1], y[0:i + 1], 'r',
                         label='Ectoderm to Endoderm')  # plots data in red
                        # label = 'Ectoderm to Endoderm Position- Short Axis (Com0R)')  # plots data in red


                # plot green curve
                x1, y1 = emb.getCurve('CoM0Y')
                x1 = 20 * (x1)  # adjusts x axis to time rather than frame
                x1 = x1[np.where(x1<time_thres)]
                y1 = y1[np.where(x1<time_thres)]
                plt.plot(x1[0:i + 1], y1[0:i + 1], 'gold',
                         label='Mesoderm to Endoderm')  # plots data in yellow

                '''plot control aves on second subplot'''

                if 'CoM0R' in emb_obj.avgCurves:
                    xR, yR, yerrR = emb_obj.avgCurves['CoM0R']
                    xR = xR * 20
                    upper = np.add(yR[np.where(xR<time_thres)], yerrR[np.where(xR<time_thres)])
                    lower = np.subtract(yR[np.where(xR<time_thres)], yerrR[np.where(xR<time_thres)])
                    plt.plot(xR[np.where(xR<time_thres)], yR[np.where(xR<500)], 'darksalmon')
                    plt.fill_between(xR[np.where(xR<time_thres)], upper, lower, color='0.8', alpha=0.5)  # plots with shaded error bars

                if 'CoM0Y' in emb_obj.avgCurves:
                    xY, yY, yerrY = emb_obj.avgCurves['CoM0Y']
                    xY = xY * 20
                    upper = np.add(yY[np.where(xY<500)], yerrY[np.where(xY<500)])
                    lower = np.subtract(yY[np.where(xY<500)], yerrY[np.where(xY<500)])
                    plt.plot(xY[np.where(xY<500)], yY[np.where(xY<500)], 'khaki')
                    plt.fill_between(xY[np.where(xY<500)], upper, lower, color='0.8', alpha=0.5)  # plots with shaded error bars
                    plt.subplots_adjust(wspace=0.35, hspace=0.5)
                plt.legend(loc=1, prop={'size': 8}, frameon=False)  # loc=1 is upper right, 2=upper left

                # if 'CoM0Y' in emb_obj.avgCurves:
                #     xY, yY, yerrY = emb_obj.avgCurves['CoM0Y']
                #     xY = xY * 20
                #     upper = np.add(yY, yerrY)
                #     lower = np.subtract(yY, yerrY)
                #     plt.plot(xY, yY, 'khaki')
                #     plt.fill_between(xY, upper, lower, color='0.8', alpha=0.5)  # plots with shaded error bars
                #     plt.subplots_adjust(wspace=0.35, hspace=0.5)

                if showFlag:
                    plt.show()
                    mainfig.show()
                if saveFlag:
                    mainfig.savefig(
                        'Z:/EMBD{n:04}_{s}_{t}.tif'.format(n=embd,
                                                                                                                      s=strain,
                                                                                                                         t=i))

def make_allplot_tif_series(embd, embryo, strain='MS', blackback=False):
    '''FIXEME-- IN PROGRESS
    # Takes in EMBD# (i.e. 0403), embryo number (i.e. 2) and strain and generates a fig with  a series of graphs that shows the control
    #  traces in grey and adds in the data for the specified embryo by timepoint to create a real-time tif series that can
    #  be made into a movie. Shows each step (if ShowFlag) and saves (if saveFlag) tif series. Works with both strains
    # :param embd: embd number (5)
    # :param embryo: embryo number (1)
    # :param strain: 'GLS' or 'MS'
    # :return:shows and/or saves an image or series of plots for  either GLS or MS for a single embryo
    '''

    from Embryos import getCenterOfMass, MSEmbryo
    import matplotlib
    from matplotlib import pyplot as plt
    import matplotlib.gridspec as gridspec
    #
    saveFlag = True
    showFlag = False

    no_yscale = True  # allows you to plot without y-scale (if true, plots without yscale-- if false, plots with y scale)
    if no_yscale:
        print("----Plotting WITHOUT Y-Scale----")

    if strain =='MS':
        folders = [FOLDER_IN + 'cropped/EMBD{n:04}/MS/Emb{x}/'.format(n=embd, x=embryo)]
        for folder in folders[:]:
            emb = MSEmbryo(folder, check_version=False)  # loads embryo object

        #
        # for i in range(1, 32):  # steps through the timepoints to plot
        '''generate the figure'''
        # fig = plt.figure(figsize=(5, 8))  # creates figure

        subplots = ['tIntR', 'tIntG', 'headInt', 'lengthR', 'MoI0G', 'MoI1G', 'MoI0R', 'MoI1R', 'CoM0G', 'CoM1G',
                    'CoM0R',
                    'CoM1R']  # defines curves to plot
        labels = [r'Total Intensity Neurons $(I^R_{tot}$)', r'Total Intensity Epidermis $(I^G_{tot}$)',
                  r'Area of Head Epidermis $(\bar A^G_{head})$', r'Neuronal Length $(L^R)$',
                  r'Epidermal Shape- Short $(\bar K^G_{short})$',
                  r'Epidermal Shape- Long $(\bar K^G_{long})$', r'Neuronal Shape- Short $(\bar K^R_{short})$',
                  r'Neuronal Shape- Long $(\bar K^R_{long})$',
                  r'Epidermal Position- Short $(\bar D^G_{short})$',
                  r'Epidermal Position- Long $(\bar D^G_{long})$',
                  r'Neuronal Position- Short $(\bar D^R_{short})$',
                  r'Neuronal Position- Long $(\bar D^R_{long})$']

        xtick_labels = ['-200', '-100', '0', '100', '200', '300', '400', '500', '600']
        columns = 3
        if len(subplots) % columns > 0:  # rounds up to include sufficient rows to accomodate all plots
            rows = len(subplots) / columns + 1
        else:
            rows = len(subplots) / columns
    #
        if plt.get_fignums():  # checks to see if figure exists already, so multiple embs can be plotted at once
            fig = plt.figure(num=1)
        else:
            fig = plt.figure(figsize=(12, 12))  # creates figure
    #
        name = getGeneCommonName(embd)
        fig.suptitle("{g}(RNAi)\n EMBD{n:04} Emb{x}  Morphogenesis Strain".format(g=name, n=embd, x=embryo),
                     fontsize=12, style='italic',
                     y=.97)  # adds title to the top of all subplots

        time_step = [w for w in range(32)]

        for w in time_step:
            for i in range(len(subplots)):  # toggles through subplot index and plots as it goes
                time_thres = 450
                plt.subplot(rows, columns, i + 1)  # defines the subplot to draw on
                plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
                plt.title(labels[i], fontsize=12, y=1.05)  # makes subplot title
                plt.xlabel('Time (minutes)', fontsize=11)  # labels x-axis
                ax = plt.subplot(rows, columns, i + 1)
                ax.tick_params(axis='x', bottom=True, top=False, left=False, right=False,
                               labelsize=8)  # get rid of ticks on top and right
                ax.tick_params(axis='y', bottom=False, top=False, left=True, right=False, labelsize=8)

                ax.yaxis.get_offset_text().set_fontsize(8)
                x, y = emb.getCurve(subplots[i])
                x = 20 * x
                if no_yscale:
                    y = y / emb.scale[subplots[i]]

                x = x[np.where(x < time_thres)]
                y = y[np.where(x < time_thres)]

                plt.plot(x, y, 'w', label=labels[i], alpha=0.5) # plots in white to set the axis for the final plot
                if w < len(x) and not np.isnan(y[w]): #plots in red frame by frame
                    z=0
                    plt.plot(x[z:time_step[w+1]],y[z:time_step[w+1]], 'r', label=labels[i])  # plots data in red
                    # plt.plot([x[time_step[w - 1]], x[time_step[w]]], [y[time_step[w - 1]], y[time_step[w]]], 'r',
                    #      label=labels[i])  # plots data in red
                elif w == len(x):
                    z = 0
                    plt.plot(x[z:time_step[w]], y[z:time_step[w]], 'r', label=labels[i])  # plots data in red


                #             plt.plot(x[:emb.cutoff[subplots[i]]], y[:emb.cutoff[subplots[i]]], 'k--',
    #                      dashes=(3, 2))  # plots fit-region with dashed line
    #
                if subplots[i] in emb.avgCurves:
                    xR, yR, yerrR = emb.avgCurves[subplots[i]]
                    xR = xR * 20
                    xR = xR[np.where(xR < time_thres)]  # makes the plots without late timepoints
                    yR = yR[np.where(xR < time_thres)]
                    yerrR = yerrR[np.where(xR < time_thres)]
                    upper = np.add(yR, yerrR)
                    lower = np.subtract(yR, yerrR)
                # if w <= len(xR):
                plt.plot(xR, yR, 'k')
                plt.fill_between(xR, upper, lower, color='0.8')  # plots with shaded error bars
                # plt.errorbar(xR[::10], yR[::10], yerrR[::10], color='k')  # plots with normal looking error bars
                plt.subplots_adjust(wspace=0.35, hspace=0.5)

            print(w)
            if showFlag:
                plt.show()
            if saveFlag:
                plt.savefig('Z:/EMBD{n:04}_{s}_{t}.tif'.format(n=embd,s=strain,t=w))
        #
    elif strain == 'GLS':
        folders = [FOLDER_IN + 'cropped/EMBD{n:04}/GLS/Emb{x}/'.format(n=embd, x=embryo)]
        for folder in folders[:]:
            emb = GSEmbryo(folder, check_version=False)  # loads embryo object
        curves = emb.getCurveNames()
        print(curves)

        subplots = ['spotsG', 'spotsR', 'spotsY', 'MoI0G', 'MoI1G', 'MoI0R', 'MoI1R', 'MoI0Y', 'MoI1Y', 'CoM0R',
                    'CoM1R',
                    'CoM0Y', 'CoM1Y']  # defines curves to plot
        labels = [r'Nuclear Count- Endoderm/Pharynx $(N^G_{nuc})$', r'Nuclear Count- Ectoderm $(N^R_{nuc})$',
                  r'Nuclear Count- Mesoderm $ (N^Y_{nuc})$', r'Endoderm/Pharynx Shape- Short $(\bar K^G_{short})$',
                  r'Endoderm/Pharynx Shape- Long $(\bar K^G_{long})$', r'Ectoderm Shape- Short $(\bar K^R_{short})$',
                  r'Ectoderm Shape- Long $(\bar K^R_{long})$',
                  r'Mesoderm Shape- Short $(\bar K^Y_{short})$', r'Mesoderm Shape- Long $(\bar K^Y_{long})$',
                  r'Ectoderm Position - Short $(\bar D^R_{short})$',
                  r'Ectoderm Position- Long $(\bar D^R_{long})$', r'Mesoderm Position - Short $(\bar D^Y_{short})$',
                  r'Mesoderm Position- Long $(\bar D^Y_{long})$']

        # labels = ['Number of Endoderm/Pharynx Nuclei (green)', 'Number of Ectoderm Nuclei (red)',
        #           'Number of Mesoderm Nuclei (yellow)', 'Endo/Pharynx Shape- Short Axis (MoI0G)',
        #           'Endo/Pharynx Shape- Long Axis (MoI1G)', 'Ectoderm Shape- Short Axis (MoI0R)',
        #           'Ectoderm Shape- Long Axis (MoI1R)',
        #           'Mesoderm Shape- Short Axis (MoI0Y)', 'Mesoderm Shape- Long Axis (MoI1Y)',
        #           'Ectoderm Position - Short Axis (Com0R)',
        #           'Ectoderm Position- Long Axis (Com1R)', 'Mesoderm Position - Short Axis (Com0Y)',
        #           'Mesoderm Position- Long Axis (Com1Y)']

        xtick_labels = ['-200', '-100', '0', '100', '200', '300', '400', '500', '600', '700']
        columns = 3

        if len(subplots) % columns > 0:  # rounds up to include sufficient rows to accommodate all plots
            rows = len(subplots) / columns + 1
        else:
            rows = len(subplots) / columns

        if plt.get_fignums():
            fig = plt.figure(num=1)
        else:
            fig = plt.figure(figsize=(12, 12))  # creates figure
        name = getGeneCommonName(embd)
        fig.suptitle("{g}(RNAi)\n EMBD{n:04} Emb{x}  Germ Layer Strain".format(g=name, n=embd, x=embryo), fontsize=12,
                     style='italic',
                     y=.99)  # adds title to the top of all subplots

        Gcurves = ['MoI0G',
                   'MoI1G']  # this block is in place to allow for dropping skewed values before markers turn on.
        Rcurves = ['MoI0R', 'MoI1R', 'CoM0R', 'CoM1R']
        Ycurves = ['MoI0Y', 'MoI1Y', 'CoM0Y', 'CoM1Y']
        fixGcurves = False
        fixGinds = 0
        fixRcurves = False
        fixRinds = 0
        fixYcurves = False
        fixYinds = 0

        time_step = [w for w in range(32)]

        for w in time_step:
            for i in range(len(subplots)):  # toggles through subplot index and plots as it goes
                time_thres = 450  # switch back to 550
                plt.subplot(rows, columns, i + 1)  # defines the subplot to draw on
                plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)
                plt.title(labels[i], fontsize=12, y=1.00)  # makes subplot title

                plt.xlabel('Time (minutes)', fontsize=11)  # labels x-axis
                ax = plt.subplot(rows, columns, i + 1)
                # ax.tick_params(axis='x', labelsize=8)
                # ax.tick_params(axis='y', labelsize=8)
                ax.tick_params(axis='x', bottom=True, top=False, left=False, right=False,
                               labelsize=8)  # get rid of ticks on top and right
                ax.tick_params(axis='y', bottom=False, top=False, left=True, right=False, labelsize=8)
                ax.yaxis.get_offset_text().set_fontsize(8)

                x, y = emb.getCurve(subplots[i])
                x = 20 * x
                if no_yscale:
                    y = y / emb.scale[subplots[i]]
                x = x[np.where(x < time_thres)]
                y = y[np.where(x < time_thres)]

                if subplots[i] in ['spotsG', 'spotsR',
                                   'spotsY']:  # this is a fix put in place to deal with many plots where the early timepoints had massively skewed spots counts (i.e 1000)
                    if np.any(y[0:10] > 400):
                        inds = np.where(y[0:10] > 400)
                        y = pd.DataFrame(y)
                        y.iloc[[inds]] = np.nan
                        if subplots[i] == 'spotsG':
                            fixGcurves = True
                            fixGinds = inds
                        elif subplots[i] == 'spotsR':
                            fixRcurves = True
                            fixRinds = inds
                        elif subplots[i] == 'spotsY':
                            fixYcurves = True
                            fixYinds = inds
                if subplots[i] in Gcurves and fixGcurves:
                    y = pd.DataFrame(y)
                    y.iloc[fixGinds] = np.nan
                if subplots[i] in Rcurves and fixRcurves:
                    y = pd.DataFrame(y)
                    y.iloc[fixRinds] = np.nan
                if subplots[i] in Ycurves and fixYcurves:
                    y = pd.DataFrame(y)
                    y.iloc[fixYinds] = np.nan
                plt.plot(x, y, 'w', label=labels[i], alpha=0.5) # plots in white to set the axis for the final plot
                if w < len(x) and not np.isnan(y[w]): #plots in red frame by frame
                    z=0
                    plt.plot(x[z:time_step[w+1]],y[z:time_step[w+1]], 'r', label=labels[i])  # plots data in red
                    # plt.plot([x[time_step[w - 1]], x[time_step[w]]], [y[time_step[w - 1]], y[time_step[w]]], 'r',
                    #      label=labels[i])  # plots data in red
                elif w == len(x):
                    z = 0
                    plt.plot(x[z:time_step[w]], y[z:time_step[w]], 'r', label=labels[i])  # plots data in red
                    plt.plot(x, y, 'r', label=labels[i])  # plots data in red

                if subplots[i] in emb.avgCurves:
                    xR, yR, yerrR = emb.avgCurves[subplots[i]]
                    xR = xR * 20
                    xR = xR[np.where(xR < time_thres)]  # makes the plots without late timepoints
                    yR = yR[np.where(xR < time_thres)]
                    yerrR = yerrR[np.where(xR < time_thres)]
                    upper = np.add(yR, yerrR)
                    lower = np.subtract(yR, yerrR)

                plt.plot(xR, yR, 'k')
                plt.fill_between(xR, upper, lower, color='0.8')  # plots with shaded error bars
                # plt.errorbar(xR[::10], yR[::10], yerrR[::10], color='k')  # plots with normal looking error bars
                plt.subplots_adjust(wspace=0.35, hspace=0.5)
                print(w)
            if saveFlag:
                plt.savefig('Z:/EMBD{n:04}_GLS_{s}_{t}.tif'.format(n=embd,
                                                                                                                    s=strain,
                                                                                                                   t=w))
            if showFlag:
                plt.show()

def make_long_and_short_axis_tif_series(embd, embryo, strain='MS', blackback=False):
    '''
    Takes in EMBD# (i.e. 0403), embryo number (i.e. 2) and strain and generates a tif series with two subplots. The first is a
    max projection along the long axis and (2) max projection along the short axis.
    :param embd: embd number
    :param embryo: embryo number
    :param strain: 'GLS' or 'MS'
    :param blackback: allows you to generate the tifs with a black background for presentation purposes
    :return:shows and/or saves an image or series of images +plots for CoM data for either GLS or MS
    '''

    from Embryos import MSEmbryo, GSEmbryo
    import matplotlib
    from scipy.ndimage import zoom
    from myFigure import *


    saveFlag = True
    showFlag = False

    rows = 2
    columns = 1

    if blackback:
        text_col = 'white'
    else:
        text_col = "black"

    '''this accesses the embryo'''
    if strain == 'MS':
        folders = [FOLDER_IN + 'cropped/EMBD{n:04}/MS/Emb{x}/'.format(n=embd, x=embryo)]
        for folder in folders[:]:
            emb = MSEmbryo(folder, check_version=False)
    elif strain == 'GLS':
        folders = [FOLDER_IN + 'cropped/EMBD{n:04}/GLS/Emb{x}/'.format(n=embd, x=embryo)]
        for folder in folders[:]:
            emb = GSEmbryo(folder, check_version=False)


    for i in range(1, 32):  # steps through the timepoints to plot
        '''generate the figure'''
        fig = plt.figure(figsize=(5, 8))  # creates figure
        if blackback:
            plt.rcParams.update({
                "lines.color": "white",
                "patch.edgecolor": "white",
                "text.color": "black",
                "axes.facecolor": "white",
                "axes.edgecolor": "lightgray",
                "axes.labelcolor": "white",
                "xtick.color": "white",
                "ytick.color": "white",
                "grid.color": "lightgray",
                "figure.facecolor": "black",
                "figure.edgecolor": "black",
                "savefig.facecolor": "black",
                "savefig.edgecolor": "black"})
        name = getGeneCommonName(embd)
        fig.suptitle("{g}(RNAi) embryo{e}\n EMBD{n:04} {s}".format(g=name, n=embd, e=embryo, s=strain), fontsize=14,
                     style='italic', y=.97, color=text_col)  # adds title to the top of all subplots

        if strain == 'MS':
            '''plot long axis max project in the first subplot'''
            plt.subplot(rows, columns, 1)  # defines the subplot to draw on
            plt.title('Lateral view', fontsize=16, y=1.0, color=text_col)  # makes subplot title
            plt.markerSize = 300
            t = emb.tScale * (i - emb.t0)  # this scales time to be aligned with plots


            # xe, ye, ze = emb.getEmbCenter()
            # plt.scatter([xe], [ye], color='gray',edgecolors= 'white', s=150) # plot embryo center

            im, (xG, yG), (xR, yR) = emb.getMaxProj(t)
            plt.imshow(im)
            ax = plt.subplot(rows, columns, 1)
            ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.axis('off')
            ax.set_axis_bgcolor('k')

            '''plot short axis max project in the second subplot'''
            plt.subplot(rows, columns, 2)  # defines the first subplot to draw on
            ax = plt.subplot(rows, columns, 2)
            ax.set_axis_bgcolor('k')
            ax.spines['bottom'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.tick_params(axis='x', colors='white')
            ax.tick_params(axis='y', colors='white')
            ax.yaxis.label.set_color('white')
            ax.xaxis.label.set_color('white')
            ax.set_aspect('equal')
            ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.axis('off')
            plt.title('End-on view', fontsize=16, y=1.0, color=text_col)  # makes subplot title


            # if emb.tScale>1.2: tSc = 1.2
            # elif emb.tScale <0.8: tSc = 0.8
            # else: tSc = emb.tScale

            # slide = int(np.round(t / emb.tScale + emb.t0))
            # slide = int(np.round(t) * emb.tScale + emb.t0)

            slide = int(np.round(i-1))  # frame to show or plot


            # slide = int(np.round(i / tSc))  # frame to show or plot
            if slide >31: slide = 31
            # slide= int(np.round(i))

            ims = emb.image.images[slide, :,
                  0:2]  # retrieves all z planes for frame i and channel 0, shape is (18,2,128,206)... this is (frame, z, channels,y,x)


            imProj = np.max(ims,
                            axis=3)  # makes a projection on x-axis to give a z-y display (shape is 18,2,128, which is Z, channels, y)
            imProjG = zoom(imProj[:, 0, :], (emb.dz, 1), order=0)
            imProjR = zoom(imProj[:, 1, :], (emb.dz, 1), order=0)
            imProjG = (imProjG / 256).astype('uint8')  # converts to 8-bit
            imProjR = (imProjR / 256).astype('uint8')
            rg_image = np.stack([imProjR, imProjG, np.zeros_like(imProjR)],
                                axis=2)  # makes an R,G,B stack from the grayscale images

            plt.imshow(rg_image)
            plt.legend(loc=None)
            if saveFlag:
                plt.savefig(
                    'Z:/long_short_tif_series/EMBD{n:04}_MS_{s}_{t}.tif'.format(
                        n=embd,
                        s=strain,
                        t=i))
            if showFlag:
                plt.show()

        elif strain == 'GLS':
            plt.subplot(rows, columns, 1)  # defines the subplot to draw on
            plt.title('Lateral view', fontsize=16, y=1.0, color=text_col)  # makes subplot title
            plt.markerSize = 300
            t = emb.tScale * (i - emb.t0)  # this scales time to be aligned with plots
            im, (xG, yG), (xR, yR) = emb.getMaxProj(t)
            # if bi:
            #     im[0] = im[0] > 0
            #     im[1] = im[1] > 0
            plt.imshow(im)
            ax = plt.subplot(rows, columns, 1)
            ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.axis('off')
            ax.set_axis_bgcolor('k')

            '''plot short axis max project in the second subplot'''
            plt.subplot(rows, columns, 2)  # defines the first subplot to draw on
            ax = plt.subplot(rows, columns, 2)
            ax.set_axis_bgcolor('k')
            ax.spines['bottom'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.tick_params(axis='x', colors='white')
            ax.tick_params(axis='y', colors='white')
            ax.yaxis.label.set_color('white')
            ax.xaxis.label.set_color('white')
            ax.set_aspect('equal')
            ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())
            ax.axis('off')
            plt.title('End-on view', fontsize=16, y=1.0, color=text_col)  # makes subplot title
            slide = int(np.round(i-1))  # frame to show or plot
            if slide >31: slide = 31
            ims = emb.image.images[slide, :,
                  0:2]  # retrieves all z planes for frame i and channel 0, shape is (18,2,128,206)... this is (frame, z, channels,y,x)
            imProj = np.max(ims,
                            axis=3)  # makes a projection on x-axis to give a z-y display (shape is 18,2,128, which is Z, channels, y)
            imProjG = zoom(imProj[:, 0, :], (emb.dz, 1), order=0)
            imProjR = zoom(imProj[:, 1, :], (emb.dz, 1), order=0)
            imProjG = (imProjG / 256).astype('uint8')  # converts to 8-bit
            imProjR = (imProjR / 256).astype('uint8')
            rg_image = np.stack([imProjR, imProjG, np.zeros_like(imProjR)],
                                axis=2)  # makes an R,G,B stack from the grayscale images

            plt.imshow(rg_image)
            plt.legend(loc=None)
            if saveFlag:
                plt.savefig(
                    'Z:/EMBD{n:04}_GLS_{s}_{t}.tif'.format(
                        n=embd,
                        s=strain,
                        t=i))
            if showFlag:
                plt.show()



def make_reciprocal_match_movies (seed):
    '''
    Takes in seed (EMBD number) and identifies reciprocal matches, converts output to format compatable
    with the show_best_match_movies format, and runs and saves from that function. IN PROGRESS
    :param seed:
    :return: makes a tif file of top 5 recip matches
    '''
    from db_utils_embryos import initialize_db
    from db_utils_genes import get_reciprocal_matches
    from PADassist import show_best_match_movies

    conn, curs = initialize_db()
    limit = 5
    recip_matches = get_reciprocal_matches(seed, 10, curs, thr=0.55)  #retrieves a list of recip matches [245,12,32,4, 7]

    seed_string = '{r:04}'.format(r=seed)  # converts to '0304' format for movie function
    recip_matches_string = [seed_string]
    PAD_list =[]
    for item in recip_matches:
        hit_str = '{r:04}'.format(r=item)  # converts each hit to 'EMBD0304' format for movie function
        recip_matches_string.append(hit_str)
        sql = 'SELECT id, rnai1, rnai2, PAD FROM pair_lookup WHERE rnai1 = {0} AND rnai2 = {1};'.format(seed, item)
        curs.execute(sql)
        id, rnai1, rnai2, PAD = curs.fetchall()[0]
        PAD_list.append(PAD)
    PAD_list_str = map(str,PAD_list)

    recip_matches_string = recip_matches_string[0:limit+1]
    PAD_list_str = PAD_list_str[0:limit+1]

    show_best_match_movies(recip_matches_string,seed_string,PAD_list_str, saveFlag=False, recip_saveFlag=True)


    print(recip_matches_string,seed_string,PAD_list_str)

def make_ALL_reciprocal_match_movies():
    '''
    make reciprocal match movies for the first 500 genes
    :return: saves reciprocal match movies to Database folder
    '''
    for i in range(401,504):
        make_reciprocal_match_movies(i)

def initialize():
    print("starting")
    # make_MS_fig_with_subplots(15, 9, saveFlag=True, showFlag=True)
    # make_MS_fig_with_subplots(661, 6, showFlag=True, saveFlag=True)
    # make_GLS_fig_with_subplots(661, 6, saveFlag=True, showFlag=True)
    # plot_avgRNAi_control_allSubplots(264, strain='MS', saveFlag=False, showFlag=True)
    # plot_allEmbs_allSubplots(1905, strain='GLS', showFlag=False, saveAllFlag=True)
    # for i in range(5,20):
    #     make_GLS_fig_with_subplots(441,i, saveFlag=True, showFlag=True)
    # for i in range(1,10):
    #     make_MS_fig_with_subplots(5,i, saveFlag=True, showFlag=True)

    # neuro_list = [5,6,7,16,32,45,68,78,91,115,119,238,239,255,264,268,288,291,320,356,357,359,364,386,387,388,398,408,422,435,44,441,513,518]
    # antRupt_list = [463, 117, 118, 9, 501, 133, 34] #239,440,281,219,255,330,95,386,268,264,422,222,122,32,369,386,455,242,343,141
    #
    # for embd in antRupt_list:
    #     plot_allEmbs_allSubplots(embd, strain='MS', saveAllFlag=True)

    # from fileLookup import redo_batch1, redo_batch2
    # for embd in redo_batch2:
    #     plot_allEmbs_allSubplots(embd, strain='GLS', saveAllFlag=True)

    '''plot and save all subplots, both strains'''
    for i in [1906,1907,1908,1909]:
        from fileLookup import SKIPPED

        if 'EMBD{0}'.format(i) not in SKIPPED:
            for strain in ['GLS', 'MS']:  #
                plot_allEmbs_allSubplots(i, strain=strain,
                                     saveAllFlag=True, showFlag=False)  # note this is not working for controls, but does work for RNAi conditions

        # get_adjusted_average_curve(10, 'spotsR', 'GLS')
        # list = [274]
        # for l in list:
        #         plot_avgRNAi_control_allSubplots(i, strain=strain, saveFlag=True, showFlag=True)

    # for i in range(3020, 3027):
    # plot_3color_spots(281)
    # plot_2color_tInt(359)
        # plot_color_MoI(0, 'MS')
        #
        # plot_color_CoM(0, 'GLS')

    ''' make reciprocal hit movies'''
    # seed = 212
    # make_reciprocal_match_movies(seed)

    # make_ALL_reciprocal_match_movies()
    #
    '''make heatmap'''
    # show_500gene_heatmap()



        # x_label, x_points, cv1_label, curve1, std1, cv2_label, curve2, std2 = plot_maxIntGLSVsPtArrest()
        # # x_label, x_points, cv1_label, curve1, std1, cv2_label, curve2, std2= plot_tIntVsPtArrest()
        # plot_pretty_graphs(x_label, x_points, cv1_label, curve1,std1,cv2_label, curve2, std2)
    '''make tif series'''
    # make_CoM_tif_series(5,4, strain='GLS', blackback=False) #448,1 is a good wt-like emb gls
    # make_CoM_tif_series_endon(88, 2, strain='MS', blackback=False)


        # show_head_endon(6, 1)
        # CoM_endon(6,1)

    # make_CoM_tif_series_endon(88, 1, strain='MS', blackback=False)
    # show_head_endon(95, 2, strain='MS', blackback=False)  #400, 1 for control images

    # show_manual_heatmap_by_gene()

    # show_manual_heatmap_by_gene_germline()

    '''make timeseries_all plots'''
    # make_allplot_tif_series(5, 4, strain='GLS', blackback=False)
    # make_long_and_short_axis_tif_series(5, 1, strain='MS', blackback=False)


if __name__ == '__main__':
    initialize()
