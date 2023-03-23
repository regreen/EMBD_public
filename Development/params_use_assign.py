"""
Library that assigns PARAMS_NAMES_USE_GLS and PARAMS_NAMES_USE_MS
"""

import numpy as np


def get_params_use():
    """
    assigns names for GLS & MS parameters
    :return: names of parameters
    """
    return _group_rank_opt() #*use this one!
    # return _lda_group()
    # return _all()
    # return _get_pars_rank_all()
    # return _lda_RNAi()
    # return _gls_params()
    # return _ms_params()

def _all():
    print("!!!!!!!!!USING ALL PARAMETERS!!!!!!!!!")
    return ['maxR', 'maxG', 'aG', 'sG', 'Gtail', 'aR', 'sR', 'Rtail', 'aY', 'sY', 'Ytail', 'mRmG', 'mYmG',
            'tMove'] + \
           ['fracR', 'fracG', 'fracY', 'tScale'] + \
           ['CoM0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['CoM1R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['CoM0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['CoM1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI0G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI1G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI1R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']], \
           ['maxG', 'aG', 'mG', 'sG', 'maxR', 'aR', 'mR', 'sR', 'tMove', 'mRmG', 'tScale'] + \
           ['maxHead', 'tailHead', 'scaleHead', 'devHead'] + \
           ['scaleLength', 'devLength', 'tailLength'] + \
           ['CoM0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['CoM1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['CoM0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['CoM1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]

def _gls_params():
    print("USING ONLY GLS-------------- PARAMETERS!!!!!!!!!")
    return ['maxR', 'maxG', 'aG', 'sG', 'Gtail', 'aR', 'sR', 'Rtail', 'aY', 'sY', 'Ytail', 'mRmG', 'mYmG',
            'tMove'] + \
           ['fracR', 'fracG', 'fracY', 'tScale'] + \
           ['CoM0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['CoM1R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['CoM0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['CoM1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI0G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI1G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI1R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']], []

def _ms_params():
    print("USING ONLY MS------------ PARAMETERS!!!!!!!!!")
    return [],['maxG', 'aG', 'mG', 'sG', 'maxR', 'aR', 'mR', 'sR', 'tMove', 'mRmG', 'tScale'] + \
           ['maxHead', 'tailHead', 'scaleHead', 'devHead'] + \
           ['scaleLength', 'devLength', 'tailLength'] + \
           ['CoM0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['CoM1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['CoM0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['CoM1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
           ['MoI1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]

def _lda_RNAi():
    '''optimal params based on LDA on RNAi vs. Control'''
    print("!!!!!!!!!USING LDA RNAi v Control PARAMETERS!!!!!!!!!")
    return ['Rtail'], ['CoM1RstdTail_MS', 'CoM0RstdTail_MS','MoI1RstdTail_MS','devHead']

def _lda_group():
    # """parameters are obtained after LDA on groups using th cutoff on its eigenvector values 1/24/18"""
    # th = 0.088
    # p_gls = ['MoI1YstdTail_GLS', 'MoI1GstdTail_GLS', 'MoI0YavgRes_GLS', 'MoI0YstdTail_GLS', 'maxR', 'sR', 'aR',
    #          'tScale',
    #          'MoI0GstdTail_GLS', 'MoI0Rtd_GLS', 'aG', 'Rtail', 'CoM0Ytd_GLS', 'MoI1RavgRes_GLS', 'CoM1Rtd_GLS',
    #          'CoM1Ytd_GLS', 'maxG', 'CoM0YstdTail_GLS', 'CoM0RavgRes_GLS', 'MoI0Gtd_GLS', 'MoI1Ytd_GLS', 'Gtail',
    #          'CoM1RstdTail_GLS', 'CoM1Yscale_GLS', 'tMove', 'CoM1RavgRes_GLS', 'MoI0RavgRes_GLS', 'fracG', 'Ytail',
    #          'CoM0RstdTail_GLS', 'MoI1RstdTail_GLS', 'MoI0Ytd_GLS', 'MoI1GavgRes_GLS', 'MoI0GavgRes_GLS', 'MoI1Rtd_GLS',
    #          'MoI1YavgRes_GLS', 'MoI1Gtd_GLS', 'CoM1YstdTail_GLS', 'sY', 'MoI0RstdTail_GLS', 'fracY', 'CoM0Rtd_GLS',
    #          'MoI1Rscale_GLS', 'MoI0Rscale_GLS', 'MoI1Gscale_GLS', 'CoM0YavgRes_GLS', 'mRmG', 'MoI0Yscale_GLS', 'aY',
    #          'mYmG', 'sG', 'CoM0Yscale_GLS', 'MoI1Yscale_GLS', 'CoM1Rscale_GLS', 'CoM0Rscale_GLS', 'fracR',
    #          'CoM1YavgRes_GLS', 'MoI0Gscale_GLS']
    # v_gls = [0.22862684523308058, 0.16792232185004394, 0.15491250946584145, 0.14548945688672074, 0.14359544826767479,
    #          0.14336317788850098, 0.12958008388895495, 0.12230792364119915, 0.11621145099508468, 0.11179407281168961,
    #          0.10972065164958157, 0.1006074713093692, 0.097651692123804118, 0.097456626083825465, 0.092653765996182885,
    #          0.091202393674338245, 0.091169699556468864, 0.089878385996220153, 0.086763608529237105,
    #          0.086294925783954771,
    #          0.07756664213073404, 0.074759306672032594, 0.072148183336436883, 0.068512845209481282,
    #          0.066766783064704044,
    #          0.062994254032786329, 0.062253062859336412, 0.062044660175684364, 0.061788859388522965,
    #          0.061569434307922041,
    #          0.06111827423440009, 0.060647456663547941, 0.060231614750425988, 0.057858807051929417,
    #          0.055781012010882758,
    #          0.054924833584544926, 0.053695280490711657, 0.052545075470618138, 0.050073539977837513,
    #          0.049001979708468472,
    #          0.04465647571068114, 0.041548348992208828, 0.041308362089235995, 0.041300076634846562,
    #          0.036862454921166417,
    #          0.036272572304336004, 0.035004286925127664, 0.034042030370988741, 0.033334652698454446,
    #          0.033097238173637772,
    #          0.033078699519458118, 0.032318718492014301, 0.030734302456456912, 0.03072720430120475,
    #          0.029242263239322559,
    #          0.027077645164877116, 0.020966809851794749, 0.016378894534519278]
    # p_ms = ['CoM1RstdTail_MS', 'CoM0GstdTail_MS', 'MoI0Gscale_MS', 'MoI0RstdTail_MS', 'maxHead', 'aR',
    #         'CoM0RstdTail_MS',
    #         'MoI1RstdTail_MS', 'devHead', 'CoM1GstdTail_MS', 'CoM0Gtd_MS', 'MoI1GstdTail_MS', 'CoM1GavgRes_MS',
    #         'tScale',
    #         'CoM1Rtd_MS', 'MoI1Gtd_MS', 'MoI1Gscale_MS', 'scaleHead', 'devLength', 'tailHead', 'mR', 'sG', 'MoI0Rtd_MS',
    #         'MoI0Gtd_MS', 'CoM0Gscale_MS', 'MoI0GstdTail_MS', 'MoI0RavgRes_MS', 'MoI1Rscale_MS', 'MoI0GavgRes_MS',
    #         'CoM1RavgRes_MS', 'sR', 'MoI1Rtd_MS', 'CoM0RavgRes_MS', 'CoM1Gscale_MS', 'MoI1GavgRes_MS', 'CoM1Gtd_MS',
    #         'maxR',
    #         'CoM0GavgRes_MS', 'scaleLength', 'aG', 'maxG', 'CoM0Rscale_MS', 'CoM1Rscale_MS', 'CoM0Rtd_MS',
    #         'MoI0Rscale_MS',
    #         'MoI1RavgRes_MS', 'tailLength', 'mRmG', 'tMove', 'mG']
    # v_ms = [0.29343773228483733, 0.13150641407637054, 0.12163220067099569, 0.11368628547713333, 0.1125497312283425,
    #         0.10866716323077871, 0.10442258641338353, 0.10087630670927278, 0.091964698863814989, 0.091690105830659327,
    #         0.079130294052409705, 0.071290584320351749, 0.066791964943058502, 0.066025882291725652,
    #         0.065665491167468393,
    #         0.065278983147626965, 0.063591505622943703, 0.063460221664974142, 0.061200762803162854,
    #         0.060487028766152327,
    #         0.058063261350302715, 0.056939623372920904, 0.056140127617482688, 0.055308967165431926,
    #         0.052038858711532847,
    #         0.050869053738418477, 0.047754917481341017, 0.046661548072726919, 0.046604526682922642,
    #         0.040627261285743489,
    #         0.038281604185199096, 0.037933253808961724, 0.037541801205703745, 0.037529533271163221,
    #         0.037517147168959514,
    #         0.037118234181638426, 0.033438682995446271, 0.032562692947378694, 0.032280850564926591,
    #         0.031931767053915434,
    #         0.031906044105975015, 0.030801330837760682, 0.030447626882883182, 0.029002712710494299,
    #         0.028547816315750332,
    #         0.0282282671265481, 0.023823682988020494, 0.015312550279593152, 0.015227107733229623, 0.015106561915269063]
    """parameters are obtained after LDA on groups using th cutoff on its eigenvector values 3/2/18"""
    # th = 0.0
    # p_gls = ['MoI0YstdTail_GLS', 'CoM0YavgRes_GLS', 'MoI0YavgRes_GLS', 'MoI1YstdTail_GLS', 'CoM0Ytd_GLS',
    #          'MoI1RavgRes_GLS', 'CoM0Rtd_GLS', 'MoI1Gtd_GLS', 'CoM1Ytd_GLS', 'maxR', 'MoI1GavgRes_GLS', 'tScale',
    #          'CoM1Rscale_GLS', 'CoM1Rtd_GLS', 'Ytail', 'MoI1GstdTail_GLS', 'CoM0YstdTail_GLS', 'MoI0GavgRes_GLS',
    #          'Gtail', 'MoI0Rscale_GLS', 'aG', 'MoI0Yscale_GLS', 'CoM0Rscale_GLS', 'MoI0Gtd_GLS', 'CoM1RavgRes_GLS',
    #          'MoI1RstdTail_GLS', 'fracG', 'Rtail', 'MoI0Gscale_GLS', 'MoI0RstdTail_GLS', 'MoI0RavgRes_GLS',
    #          'CoM1Yscale_GLS', 'maxG', 'MoI0Rtd_GLS', 'MoI0GstdTail_GLS', 'CoM1YavgRes_GLS', 'tMove', 'MoI1YavgRes_GLS',
    #          'aY', 'MoI1Yscale_GLS', 'CoM1RstdTail_GLS', 'sG', 'sY', 'fracY', 'CoM0Yscale_GLS', 'sR', 'MoI0Ytd_GLS',
    #          'CoM0RstdTail_GLS', 'MoI1Rtd_GLS', 'MoI1Gscale_GLS', 'MoI1Rscale_GLS', 'aR', 'mYmG', 'CoM0RavgRes_GLS',
    #          'CoM1YstdTail_GLS', 'MoI1Ytd_GLS', 'mRmG', 'fracR']
    # v_gls = [0.12358848976024053, 0.11694291941551142, 0.10067037338665243, 0.094230546308524885, 0.08638603407482863,
    #          0.080238860654955363, 0.080127718463656919, 0.07953831318697259, 0.078520125838972193,
    #          0.077421376779478565, 0.076700809490016175, 0.073270738335379448, 0.071765025776379623,
    #          0.071741824163678858, 0.069154189433268926, 0.067244937628271634, 0.06474628395337971,
    #          0.060287026455415386, 0.056149152061122788, 0.056017022686304922, 0.055678756663772207,
    #          0.051426138445763075, 0.049673519104348557, 0.0463128550374357, 0.046237389513429594, 0.046005804272494816,
    #          0.041773502407515366, 0.040038002499950265, 0.039800546707354702, 0.039300800371132254,
    #          0.039166519929804453, 0.038455173177286334, 0.036763088899249896, 0.035696798471742307,
    #          0.035613439743327133, 0.032414022277088318, 0.031472426156418688, 0.030898772764649075,
    #          0.029431030425967609, 0.029231896403001743, 0.025679909226421341, 0.025552236253075133,
    #          0.025128987426089478, 0.023765896964903586, 0.023326571978520904, 0.022771421852835958,
    #          0.022239273996535563, 0.021507615529642633, 0.021073430699338438, 0.020992560928320674,
    #          0.017559834411274376, 0.017358581090659537, 0.016703629729689438, 0.016487193694780092,
    #          0.013682861880314796, 0.013442028428013253, 0.010097854945524419, 0.0060501643718536887]
    # p_ms = ['devHead', 'CoM0GstdTail_MS', 'tailHead', 'CoM1RavgRes_MS', 'MoI1RavgRes_MS', 'MoI0RavgRes_MS',
    #         'MoI0Gscale_MS', 'CoM1RstdTail_MS', 'CoM0RstdTail_MS', 'CoM0Rtd_MS', 'sR', 'CoM1GstdTail_MS',
    #         'CoM0Gscale_MS', 'CoM0Gtd_MS', 'MoI1GstdTail_MS', 'CoM0GavgRes_MS', 'MoI1RstdTail_MS', 'CoM1Gscale_MS',
    #         'aR', 'MoI0Rtd_MS', 'MoI0RstdTail_MS', 'CoM0RavgRes_MS', 'MoI1GavgRes_MS', 'MoI1Rscale_MS', 'CoM1Rtd_MS',
    #         'maxG', 'maxR', 'MoI0GstdTail_MS', 'sG', 'CoM0Rscale_MS', 'devLength', 'MoI0Gtd_MS', 'maxHead',
    #         'MoI0GavgRes_MS', 'CoM1Rscale_MS', 'tailLength', 'scaleLength', 'tScale', 'MoI1Rtd_MS', 'MoI0Rscale_MS',
    #         'scaleHead', 'CoM1Gtd_MS', 'mR', 'aG', 'mRmG', 'CoM1GavgRes_MS', 'MoI1Gtd_MS', 'MoI1Gscale_MS', 'mG',
    #         'tMove']
    # v_ms = [0.37504524885992252, 0.28519311956762144, 0.21534070126065169, 0.14807218210003073, 0.13959157137988351,
    #         0.1248970539736497, 0.1126304652706252, 0.11095189409749832, 0.10964276609739677, 0.10655638100042966,
    #         0.10149554440001199, 0.097519178404212079, 0.092539304004950015, 0.08423680414257019, 0.084197267362921199,
    #         0.075519119051671413, 0.073578961388327771, 0.068315925459514287, 0.065292435954675654,
    #         0.063454659536040348, 0.060879712936452217, 0.060512479043650536, 0.057671485703592429,
    #         0.057422462085320654, 0.05546119625686903, 0.054774521396825558, 0.049474407409567593, 0.048862103933614821,
    #         0.043659571428217768, 0.040399435697980365, 0.039298492766910463, 0.03599133162247526, 0.032871771019355636,
    #         0.031516897499219786, 0.030669528850844101, 0.030415258502000702, 0.028772961179846599,
    #         0.028127529873993275, 0.025233846114175087, 0.021247105432917419, 0.020901331870789078, 0.01993529298039224,
    #         0.019111790013377452, 0.017065483152847619, 0.016149494989640186, 0.016049449981868767,
    #         0.015037097181178275, 0.013414739436723, 0.012961700268639485, 0.0060619865825997247]

    """parameters obtained after LDA on groups using th cutoff on its eigenvector values 3/8/18- origin reset"""
    th = 0.055
    p_gls = ['MoI0YavgRes_GLS', 'MoI0YstdTail_GLS', 'MoI1YavgRes_GLS', 'Ytail', 'MoI0Gtd_GLS', 'maxG', 'Rtail',
             'CoM0Ytd_GLS', 'sR', 'MoI1GavgRes_GLS', 'CoM0YavgRes_GLS', 'MoI1Gtd_GLS', 'aR', 'CoM0Rtd_GLS',
             'MoI1YstdTail_GLS', 'sY', 'fracG', 'CoM0YstdTail_GLS', 'MoI0GstdTail_GLS', 'MoI1Ytd_GLS', 'aY',
             'MoI0Gscale_GLS', 'maxR', 'CoM0RavgRes_GLS', 'Gtail', 'MoI0Yscale_GLS', 'MoI1GstdTail_GLS', 'tScale',
             'MoI0RstdTail_GLS', 'MoI1Rscale_GLS', 'MoI0Rtd_GLS', 'CoM0RstdTail_GLS', 'MoI0RavgRes_GLS', 'CoM1Ytd_GLS',
             'fracY', 'CoM0Yscale_GLS', 'mYmG', 'MoI1Rtd_GLS', 'MoI1Yscale_GLS', 'CoM0Rscale_GLS', 'mRmG',
             'MoI1RavgRes_GLS', 'CoM1Rscale_GLS', 'CoM1RstdTail_GLS', 'tMove', 'CoM1Rtd_GLS', 'MoI1Gscale_GLS',
             'CoM1YavgRes_GLS', 'MoI0GavgRes_GLS', 'MoI0Ytd_GLS', 'MoI0Rscale_GLS', 'CoM1YstdTail_GLS',
             'MoI1RstdTail_GLS', 'aG', 'CoM1RavgRes_GLS', 'sG', 'fracR', 'CoM1Yscale_GLS']
    v_gls = [0.29823035392553759, 0.1421347891815633, 0.14072971290052466, 0.11693545424339889, 0.11336956363945554,
           0.10022759027895517, 0.095885509446955536, 0.080459419744971325, 0.073756808047627506, 0.073186491330521866,
           0.070812291238252667, 0.068424131366836768, 0.067849260829824717, 0.067800556343059495, 0.067686774258687535,
           0.062910770208758071, 0.056995943884258403, 0.052449449252897115, 0.052156420586540554, 0.050981592119708471,
           0.044371803892961199, 0.042539497789123457, 0.04161250115932693, 0.040023182247361511, 0.039564684728712607,
           0.038917057889446348, 0.03842221759119413, 0.036733932669208956, 0.035832309258601203, 0.035766366017486458,
           0.034604611727751207, 0.033785920708759784, 0.03334844174495788, 0.032499549137162534, 0.032444118723116984,
           0.031550788649843516, 0.030138579892527578, 0.029645777611712169, 0.029644382818629748, 0.028964185398265813,
           0.02871686787212771, 0.028691856241008343, 0.028324684710344546, 0.02776899542838776, 0.027106628702398591,
           0.02584509371790783, 0.024937711263752178, 0.024908645267453129, 0.02430631263938051, 0.02283468626997897,
           0.019181768970258153, 0.016171222340859336, 0.015838190531326676, 0.014566300790568735, 0.014374103793500062,
           0.01201228355636238, 0.011833571081746284, 0.0087066852690607949]

    p_ms = ['devHead', 'CoM0GstdTail_MS', 'CoM1GstdTail_MS', 'tailHead', 'CoM1RstdTail_MS', 'CoM0RstdTail_MS',
             'CoM0RavgRes_MS', 'MoI1RstdTail_MS', 'aG', 'CoM0GavgRes_MS', 'MoI0GstdTail_MS', 'MoI0RstdTail_MS', 'maxR',
             'CoM1RavgRes_MS', 'MoI1Gscale_MS', 'aR', 'CoM0Gscale_MS', 'MoI1RavgRes_MS', 'MoI1GavgRes_MS',
             'CoM0Rscale_MS', 'MoI1GstdTail_MS', 'MoI0RavgRes_MS', 'sR', 'maxG', 'tScale', 'MoI1Rtd_MS', 'CoM0Rtd_MS',
             'MoI0Gscale_MS', 'tailLength', 'MoI0GavgRes_MS', 'MoI1Rscale_MS', 'CoM1Gscale_MS', 'mRmG', 'devLength',
             'maxHead', 'scaleLength', 'CoM1Rtd_MS', 'scaleHead', 'MoI0Rtd_MS', 'mR', 'CoM1GavgRes_MS', 'sG', 'mG',
             'CoM1Rscale_MS', 'MoI0Gtd_MS', 'MoI1Gtd_MS', 'MoI0Rscale_MS', 'CoM1Gtd_MS', 'tMove', 'CoM0Gtd_MS']

    v_ms = [0.37883067205066695, 0.21597978899855869, 0.19248219171214037, 0.18049988014320054, 0.12482024270868899,
           0.12316507084605793, 0.09858341330038245, 0.097536539185328272, 0.090079348959183944, 0.083301800459458591,
           0.081093849407204405, 0.077906866237797856, 0.072768202916958527, 0.070622253274024591, 0.067473327166272079,
           0.064863143351408151, 0.062355223074682734, 0.050035551176754435, 0.049537706913244654, 0.048096778155381985,
           0.046572673967095461, 0.046460220115730401, 0.045932327261501274, 0.044422629895322058, 0.039863082058108054,
           0.038237042537036928, 0.036065808665949294, 0.03538945928895023, 0.032316338601718127, 0.029551408456978165,
           0.027973015280028216, 0.026357170995809576, 0.026029739524535273, 0.02586035833673132, 0.02464108807357698,
           0.024231757482836874, 0.023121205015504631, 0.022764697938818829, 0.020834070320024226, 0.017690312700166862,
           0.017218532231207, 0.016854213109864447, 0.01657347902235522, 0.013118782460133774, 0.011767328691949071,
           0.010764960988070414, 0.010498745889696535, 0.010352413232187939, 0.0064979743917378376,
           0.0056051950771095109]


    print(len(list(np.array(p_gls)[np.array(v_gls) > th]))), (len(list(np.array(p_ms)[np.array(v_ms) > th])))
    return list(np.array(p_gls)[np.array(v_gls) > th]), list(np.array(p_ms)[np.array(v_ms) > th])




def _group_rank_opt():
    """ date = 3/7/18
    optimized parameters (using PAD) for average group rank normalized by number morph and cfs groups number
    Note: used training set (split rnaiList 70/30 train/test) to generate and required at least 3 members of each group
    for training-- origin was reset prior to this run.
    train list = [[21, 364, 408, 130], [359, 386, 422, 264, 118, 388], [417, 115, 64], [63, 19, 501, 77], [181, 182, 357], [363, 108, 154, 45, 117], [34, 16, 31, 264], [45, 52, 398, 38], [95, 277, 4], [503, 403, 404, 67, 90], [261, 453, 379, 420, 327, 225, 289], [10, 177, 217, 28, 18, 15], [110, 101, 98, 142], [321, 396, 498, 495], [387, 385, 422, 31], [320, 58, 125, 288], [426, 197, 1, 76]]
    """

    th = 4  # NOTE DIRECTION OF > < in PARAMS USE LINES

    ordered_names = ['CoM0RstdTail_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'CoM0RavgRes_GLS_GLS',
           'CoM1Ytd_GLS_GLS', 'mRmG_GLS', 'mR_MS', 'Ytail_GLS', 'aY_GLS',
           'MoI1YstdTail_GLS_GLS', 'mRmG_MS', 'tScale_GLS',
           'MoI1GavgRes_MS_MS', 'CoM1RstdTail_GLS_GLS', 'MoI0GavgRes_MS_MS',
           'CoM1YstdTail_GLS_GLS', 'sG_GLS', 'scaleLength_MS',
           'CoM0Yscale_GLS_GLS', 'scaleHead_MS', 'tMove_MS', 'mG_MS',
           'CoM0Rscale_GLS_GLS', 'CoM1Rscale_MS_MS', 'MoI1RstdTail_GLS_GLS',
           'CoM1Rtd_GLS_GLS', 'Gtail_GLS', 'MoI0Gtd_GLS_GLS',
           'MoI0Yscale_GLS_GLS', 'tScale_MS', 'MoI1GavgRes_GLS_GLS',
           'Rtail_GLS', 'CoM1Gtd_MS_MS', 'CoM0YstdTail_GLS_GLS',
           'MoI0Gtd_MS_MS', 'CoM0YavgRes_GLS_GLS', 'MoI0Rtd_MS_MS', 'mYmG_GLS',
           'MoI0Rscale_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'MoI0RstdTail_GLS_GLS',
           'MoI1Ytd_GLS_GLS', 'maxG_MS', 'CoM1RavgRes_GLS_GLS',
           'MoI0GstdTail_GLS_GLS', 'maxR_GLS', 'CoM0Ytd_GLS_GLS',
           'CoM0RstdTail_MS_MS', 'aR_GLS', 'sR_GLS', 'CoM1GstdTail_MS_MS',
           'CoM1YavgRes_GLS_GLS', 'CoM0GavgRes_MS_MS', 'MoI0Rscale_MS_MS',
           'MoI1Gtd_MS_MS', 'MoI0GstdTail_MS_MS', 'CoM0GstdTail_MS_MS',
           'fracY_GLS', 'MoI1Rtd_MS_MS', 'MoI0YstdTail_GLS_GLS',
           'CoM1Yscale_GLS_GLS', 'MoI0RavgRes_MS_MS', 'MoI1RavgRes_GLS_GLS',
           'tMove_GLS', 'CoM0Rtd_MS_MS', 'MoI0RstdTail_MS_MS', 'aG_MS',
           'MoI1Gtd_GLS_GLS', 'devHead_MS', 'tailLength_MS',
           'MoI0RavgRes_GLS_GLS', 'MoI0GavgRes_GLS_GLS', 'MoI1Yscale_GLS_GLS',
           'MoI1Rscale_MS_MS', 'MoI0Ytd_GLS_GLS', 'CoM1RstdTail_MS_MS',
           'CoM0Gscale_MS_MS', 'CoM1Rtd_MS_MS', 'MoI1GstdTail_MS_MS',
           'maxG_GLS', 'MoI0Rtd_GLS_GLS', 'sG_MS', 'fracG_GLS',
           'CoM0Rtd_GLS_GLS', 'tailHead_MS', 'MoI1YavgRes_GLS_GLS',
           'maxHead_MS', 'CoM0Rscale_MS_MS', 'CoM1GavgRes_MS_MS',
           'MoI1Gscale_GLS_GLS', 'CoM0RavgRes_MS_MS']

    ordered_vals = [18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
           18, 18, 18, 18, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15,
           15, 15, 15, 14, 13, 13, 13, 12, 11, 11, 11, 11, 10, 10, 10, 10, 10,
           9, 9, 9, 9, 8, 8, 7, 7, 6, 6, 6, 6, 6, 6, 6, 5, 4,
           3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1]


    """ date = 3/4/18
    optimized parameters for average group rank normalized by number morph and cfs groups number
    Note: used training set (split rnaiList 70/30 train/test) to generate and required at least 3 members of each group
    for training
    """
    # th = 4  ###NOTE DIRECTION OF > < in PARAMS USE LINES
    # ordered_names = ['mG_MS', 'MoI0GavgRes_MS_MS', 'tScale_MS', 'tScale_GLS',
    #                  'MoI1YstdTail_GLS_GLS', 'CoM1YstdTail_GLS_GLS',
    #                  'CoM1RstdTail_GLS_GLS', 'sG_GLS', 'mR_MS', 'MoI0Rtd_MS_MS',
    #                  'scaleHead_MS', 'tMove_MS', 'CoM0RstdTail_GLS_GLS',
    #                  'MoI1Rscale_GLS_GLS', 'CoM1Gtd_MS_MS', 'mRmG_MS', 'mRmG_GLS',
    #                  'aY_GLS', 'CoM0Yscale_GLS_GLS', 'MoI1GavgRes_MS_MS',
    #                  'MoI0Yscale_GLS_GLS', 'devHead_MS', 'Ytail_GLS', 'CoM1Ytd_GLS_GLS',
    #                  'MoI1RstdTail_GLS_GLS', 'scaleLength_MS', 'Gtail_GLS',
    #                  'MoI0Gtd_MS_MS', 'mYmG_GLS', 'Rtail_GLS', 'CoM1Rtd_GLS_GLS',
    #                  'MoI0Gtd_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'CoM0Rtd_MS_MS',
    #                  'tailHead_MS', 'MoI1GavgRes_GLS_GLS', 'MoI1Gtd_MS_MS',
    #                  'CoM0YavgRes_GLS_GLS', 'CoM1Rscale_MS_MS', 'CoM0RavgRes_GLS_GLS',
    #                  'MoI0Gscale_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'MoI0Rscale_GLS_GLS',
    #                  'tMove_GLS', 'fracY_GLS', 'MoI0RstdTail_GLS_GLS', 'MoI1Ytd_GLS_GLS',
    #                  'CoM0YstdTail_GLS_GLS', 'sR_GLS', 'MoI0Rscale_MS_MS',
    #                  'MoI0GstdTail_GLS_GLS', 'MoI1Rtd_MS_MS', 'CoM1YavgRes_GLS_GLS',
    #                  'MoI0RavgRes_MS_MS', 'aR_GLS', 'CoM0Ytd_GLS_GLS',
    #                  'MoI0YstdTail_GLS_GLS', 'CoM0GstdTail_MS_MS', 'MoI1Rscale_MS_MS',
    #                  'maxG_MS', 'maxR_GLS', 'CoM1GstdTail_MS_MS', 'CoM1Rtd_MS_MS',
    #                  'MoI0GstdTail_MS_MS', 'MoI1Yscale_GLS_GLS', 'MoI1RavgRes_GLS_GLS',
    #                  'MoI0RstdTail_MS_MS', 'CoM0GavgRes_MS_MS', 'MoI0RavgRes_GLS_GLS',
    #                  'CoM1Gscale_MS_MS', 'tailLength_MS', 'MoI0Ytd_GLS_GLS', 'aG_MS',
    #                  'MoI1Gtd_GLS_GLS', 'CoM1Yscale_GLS_GLS', 'sG_MS',
    #                  'CoM0RstdTail_MS_MS', 'CoM0Gscale_MS_MS', 'MoI1GstdTail_GLS_GLS',
    #                  'MoI0Rtd_GLS_GLS', 'CoM0RavgRes_MS_MS', 'MoI0GavgRes_GLS_GLS',
    #                  'CoM0Rscale_MS_MS', 'MoI1RstdTail_MS_MS', 'CoM0Gtd_MS_MS',
    #                  'CoM1Rscale_GLS_GLS', 'CoM1GavgRes_MS_MS', 'CoM0Rtd_GLS_GLS',
    #                  'MoI1Gscale_GLS_GLS', 'MoI1YavgRes_GLS_GLS']
    #
    # ordered_vals = [18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 17,
    #                 17, 17, 17, 17, 17, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 15, 15,
    #                 15, 14, 14, 13, 13, 13, 12, 12, 12, 11, 11, 11, 11, 11, 11, 10, 10,
    #                 10, 10, 10, 10, 9, 9, 8, 8, 8, 8, 8, 8, 6, 6, 6, 5, 5,
    #                 5, 4, 4, 4, 4, 4, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 1,
    #                 1, 1, 1, 1, 1]

    """ date = 3/2/18
    optimized parameters for average group rank normalized by number morph and cfs groups number
    Note: used training set (split rnaiList 50/50) to generate
    """

    # th = 17 ###NOTE DIRECTION OF > < in PARAMS USE LINES
    #
    # ordered_vals = [17, 17, 17, 17, 17, 17, 17, 17, 16, 16, 16, 16, 16, 15, 15, 15, 15,
    #    15, 15, 15, 15, 15, 14, 14, 14, 14, 14, 14, 14, 13, 13, 13, 13, 12,
    #    12, 12, 11, 11, 10, 10,  9,  9,  8,  8,  8,  8,  8,  8,  8,  7,  7,
    #     7,  7,  6,  6,  6,  6,  5,  5,  4,  4,  4,  4,  4,  4,  4,  3,  3,
    #     3,  3,  3,  3,  3,  3,  2,  2,  2,  2,  2,  2,  1,  1,  1,  1,  1,
    #     1,  1,  1,  1,  1]
    #
    # ordered_names = ['tScale_GLS', 'CoM1YstdTail_GLS_GLS', 'aY_GLS',
    #    'CoM0RstdTail_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'Ytail_GLS',
    #    'MoI0Gtd_MS_MS', 'mYmG_GLS', 'mG_MS', 'MoI0GavgRes_MS_MS',
    #    'fracY_GLS', 'mRmG_GLS', 'Rtail_GLS', 'MoI0Rtd_GLS_GLS', 'mR_MS',
    #    'MoI1GavgRes_MS_MS', 'aR_GLS', 'CoM1YavgRes_GLS_GLS',
    #    'MoI0Rtd_MS_MS', 'mRmG_MS', 'CoM1Gtd_MS_MS', 'MoI1Gtd_GLS_GLS',
    #    'MoI1Rscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS',
    #    'CoM0YstdTail_GLS_GLS', 'MoI1YstdTail_GLS_GLS', 'CoM0Ytd_GLS_GLS',
    #    'maxR_GLS', 'MoI0GstdTail_GLS_GLS', 'MoI1RstdTail_GLS_GLS',
    #    'CoM1Rscale_MS_MS', 'MoI1Rscale_MS_MS', 'sG_GLS',
    #    'CoM0Rscale_GLS_GLS', 'MoI1Gtd_MS_MS', 'MoI1YavgRes_GLS_GLS',
    #    'maxG_GLS', 'Gtail_GLS', 'MoI1Ytd_GLS_GLS', 'CoM0RstdTail_MS_MS',
    #    'MoI1GstdTail_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'MoI0Gtd_GLS_GLS',
    #    'tailLength_MS', 'CoM1Ytd_GLS_GLS', 'CoM1RstdTail_MS_MS',
    #    'MoI1GavgRes_GLS_GLS', 'MoI1Rtd_GLS_GLS', 'scaleHead_MS', 'aG_GLS',
    #    'scaleLength_MS', 'tMove_MS', 'tScale_MS', 'MoI1Yscale_GLS_GLS',
    #    'CoM1GstdTail_MS_MS', 'CoM0GstdTail_MS_MS', 'CoM1Rtd_GLS_GLS',
    #    'CoM1Rtd_MS_MS', 'MoI0GstdTail_MS_MS', 'CoM1Gscale_MS_MS',
    #    'MoI0Gscale_GLS_GLS', 'aG_MS', 'tMove_GLS', 'devHead_MS',
    #    'CoM0Rtd_GLS_GLS', 'MoI1Rtd_MS_MS', 'sY_GLS',
    #    'MoI0RstdTail_GLS_GLS', 'sR_GLS', 'CoM1GavgRes_MS_MS',
    #    'MoI0YavgRes_GLS_GLS', 'MoI0RavgRes_GLS_GLS', 'tailHead_MS',
    #    'fracR_GLS', 'MoI0Rscale_GLS_GLS', 'MoI0RstdTail_MS_MS',
    #    'fracG_GLS', 'devLength_MS', 'CoM0Gtd_MS_MS', 'CoM0YavgRes_GLS_GLS',
    #    'MoI0YstdTail_GLS_GLS', 'CoM0Rscale_MS_MS', 'MoI0RavgRes_MS_MS',
    #    'CoM1Yscale_GLS_GLS', 'MoI1Gscale_GLS_GLS', 'MoI1RavgRes_MS_MS',
    #    'maxHead_MS', 'MoI1RavgRes_GLS_GLS', 'MoI1Gscale_MS_MS',
    #    'MoI0Gscale_MS_MS']

    # th = 1 ###NOTE DIRECTION OF > < in PARAMS USE LINES
    # """date = 1/7/18 optimized parameters for average group rank normalized by number morph and cfs groups number
    # Note: added in 3 new groups from previous optimization method"""
    # ordered_vals = [17,17,17,17,17,17,17,17,17,17,17,16,16,16,16,15,15,15,15,15,15,15,14,14,14,14,14,14,13,13,13,12,12,11,11,11,11,11,10,9,9,8,8,8,8,7,6,6,6,5,5,5,5,4,4,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    # ordered_names = ['tScale_MS','MoI0GavgRes_MS_MS','mRmG_MS','tMove_MS','CoM0Ytd_GLS_GLS','MoI0Rscale_MS_MS','mG_MS','CoM0YavgRes_GLS_GLS','CoM1YstdTail_GLS_GLS','CoM1RstdTail_GLS_GLS','MoI1Ytd_GLS_GLS','MoI1Rscale_MS_MS','mRmG_GLS','MoI1RstdTail_GLS_GLS','MoI1Gtd_MS_MS','scaleLength_MS','CoM0RstdTail_GLS_GLS','CoM1Rscale_MS_MS','CoM0Gtd_MS_MS','MoI0Yscale_GLS_GLS','tMove_GLS','MoI0Rtd_MS_MS','CoM0Rscale_GLS_GLS','CoM0GavgRes_MS_MS','CoM1Rscale_GLS_GLS','CoM1Rtd_MS_MS','MoI1Gscale_GLS_GLS','MoI1Rtd_MS_MS','devLength_MS','CoM0RavgRes_GLS_GLS','MoI0GstdTail_MS_MS','MoI1Rscale_GLS_GLS','MoI1GstdTail_MS_MS','MoI0Ytd_GLS_GLS','sR_GLS','MoI0GstdTail_GLS_GLS','MoI1Rtd_GLS_GLS','tScale_GLS','CoM0RavgRes_MS_MS','MoI1Gtd_GLS_GLS','MoI1YstdTail_GLS_GLS','CoM1Ytd_GLS_GLS','CoM1GavgRes_MS_MS','CoM0GstdTail_MS_MS','MoI0RstdTail_MS_MS','MoI0YstdTail_GLS_GLS','CoM1Gscale_MS_MS','sG_GLS','fracY_GLS','sG_MS','mR_MS','MoI1RavgRes_GLS_GLS','MoI0RavgRes_MS_MS','CoM0YstdTail_GLS_GLS','CoM1GstdTail_MS_MS','CoM0Rtd_GLS_GLS','tailHead_MS','CoM1RavgRes_GLS_GLS','aY_GLS','CoM1Rtd_GLS_GLS','CoM0Yscale_GLS_GLS','mYmG_GLS','MoI0RstdTail_GLS_GLS','fracG_GLS','MoI1Yscale_GLS_GLS','MoI0Rtd_GLS_GLS','scaleHead_MS','MoI0RavgRes_GLS_GLS','MoI0Gtd_GLS_GLS','CoM0Gscale_MS_MS','MoI1GavgRes_MS_MS','MoI0Gscale_GLS_GLS','MoI0Rscale_GLS_GLS','CoM1RavgRes_MS_MS','MoI1GstdTail_GLS_GLS','sY_GLS','CoM1YavgRes_GLS_GLS','CoM0Rtd_MS_MS','MoI1GavgRes_GLS_GLS','Ytail_GLS','CoM1Yscale_GLS_GLS','tailLength_MS','aG_MS','maxR_GLS']
    gls_ordered_names = []
    gls_ordered_vals = []
    ms_ordered_names = []
    ms_ordered_vals = []

    for i in range(len(ordered_names)):
        name = ordered_names[i]
        if name[-2:] == 'LS':
            gls_ordered_names.append(name[:-4])
            gls_ordered_vals.append(ordered_vals[i])
        else:
            ms_ordered_names.append(name[:-3])
            ms_ordered_vals.append(ordered_vals[i])

    PARAM_REMOVE_GLS = list(np.array(gls_ordered_names)[np.array(gls_ordered_vals) > th])
    PARAM_REMOVE_MS = list(np.array(ms_ordered_names)[np.array(ms_ordered_vals) > th])
    PARAM_NAMES_CONSIDER_GLS = ['maxR', 'maxG', 'aG', 'sG', 'Gtail', 'aR', 'sR', 'Rtail', 'aY', 'sY', 'Ytail', 'mRmG',
                                'mYmG', 'tMove'] + \
                               ['fracR', 'fracG', 'fracY', 'tScale'] + \
                               ['CoM0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['CoM1R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['CoM0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['CoM1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['MoI0G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['MoI1G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['MoI0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['MoI1R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['MoI0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['MoI1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]  # +\

    PARAM_NAMES_CONSIDER_MS = ['maxG', 'aG', 'mG', 'sG', \
                               'maxR', 'aR', 'mR', 'sR', \
                               'tMove', 'mRmG', 'tScale'] + \
                              ['maxHead', 'tailHead', 'scaleHead', 'devHead'] + \
                              ['scaleLength', 'devLength', 'tailLength'] + \
                              ['CoM0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['CoM1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['CoM0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['CoM1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['MoI0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['MoI1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['MoI0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['MoI1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]  # +\

    PARAM_NAMES_USE_GLS = []
    PARAM_NAMES_USE_MS = []

    for i in PARAM_NAMES_CONSIDER_GLS:
        if i not in PARAM_REMOVE_GLS:
            PARAM_NAMES_USE_GLS.append(i)

    for i in PARAM_NAMES_CONSIDER_MS:
        if i not in PARAM_REMOVE_MS:
            PARAM_NAMES_USE_MS.append(i)
    print(PARAM_NAMES_USE_GLS)
    print(PARAM_NAMES_USE_MS)
    print (len(PARAM_NAMES_USE_GLS), len(PARAM_NAMES_USE_MS))
    return PARAM_NAMES_USE_GLS, PARAM_NAMES_USE_MS


def _get_pars_rank_all():
    from embdFunc import split_params
    from varLookup import PARAM_NAMES
    p_pop = ['mG_MS', 'CoM0RstdTail_GLS_GLS', 'scaleHead_MS', 'MoI1GavgRes_MS_MS', 'mRmG_MS', 'MoI1Rscale_GLS_GLS',
             'MoI0Yscale_GLS_GLS', 'tMove_MS', 'CoM1Gtd_MS_MS', 'CoM0Yscale_GLS_GLS', 'MoI0GavgRes_MS_MS', 'mR_MS',
             'sG_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'devHead_MS', 'CoM0Ytd_GLS_GLS',
             'CoM1YstdTail_GLS_GLS', 'CoM0Rtd_MS_MS', 'tMove_GLS', 'MoI1YstdTail_GLS_GLS', 'mRmG_GLS', 'aY_GLS',
             'MoI1GavgRes_GLS_GLS', 'CoM0GavgRes_MS_MS', 'tScale_GLS', 'tScale_MS', 'CoM1Ytd_GLS_GLS',
             'MoI1RstdTail_GLS_GLS', 'scaleLength_MS', 'CoM0YavgRes_GLS_GLS', 'Gtail_GLS', 'MoI0Gtd_MS_MS',
             'CoM1Rscale_MS_MS', 'MoI0RstdTail_GLS_GLS', 'MoI0Rscale_MS_MS', 'mYmG_GLS', 'Rtail_GLS', 'tailHead_MS',
             'maxR_GLS', 'MoI0Rtd_MS_MS', 'CoM0RavgRes_GLS_GLS', 'MoI0Gscale_GLS_GLS', 'CoM1YavgRes_GLS_GLS',
             'CoM1Rtd_GLS_GLS', 'Ytail_GLS', 'aR_GLS', 'MoI0Rscale_GLS_GLS', 'MoI1Gtd_MS_MS', 'sR_GLS',
             'MoI0RstdTail_MS_MS', 'CoM1Yscale_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'CoM0YstdTail_GLS_GLS',
             'CoM1Gscale_MS_MS', 'MoI0GstdTail_GLS_GLS', 'MoI1RavgRes_GLS_GLS', 'maxG_MS', 'CoM0RstdTail_MS_MS',
             'MoI0GstdTail_MS_MS', 'MoI0Gtd_GLS_GLS', 'CoM1GstdTail_MS_MS', 'MoI0RavgRes_MS_MS', 'CoM0GstdTail_MS_MS',
             'CoM1Rtd_MS_MS', 'MoI1Ytd_GLS_GLS']
    PARAM_REMOVE_GLS, PARAM_REMOVE_MS = [], []
    for pn in p_pop:
        if pn[-3:] == '_MS':
            PARAM_REMOVE_MS.append(pn[:-3])
        else:
            PARAM_REMOVE_GLS.append(pn[:-4])
    PARAM_NAMES_CONSIDER_GLS = ['maxR', 'maxG', 'aG', 'sG', 'Gtail', 'aR', 'sR', 'Rtail', 'aY', 'sY', 'Ytail', 'mRmG',
                                'mYmG', 'tMove'] + \
                               ['fracR', 'fracG', 'fracY', 'tScale'] + \
                               ['CoM0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['CoM1R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['CoM0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['CoM1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['MoI0G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['MoI1G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['MoI0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['MoI1R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['MoI0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                               ['MoI1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]  # +\

    PARAM_NAMES_CONSIDER_MS = ['maxG', 'aG', 'mG', 'sG', \
                               'maxR', 'aR', 'mR', 'sR', \
                               'tMove', 'mRmG', 'tScale'] + \
                              ['maxHead', 'tailHead', 'scaleHead', 'devHead'] + \
                              ['scaleLength', 'devLength', 'tailLength'] + \
                              ['CoM0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['CoM1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['CoM0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['CoM1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['MoI0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['MoI1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['MoI0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']] + \
                              ['MoI1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]  # +\

    PARAM_NAMES_USE_GLS = []
    PARAM_NAMES_USE_MS = []

    for i in PARAM_NAMES_CONSIDER_GLS:
        if i not in PARAM_REMOVE_GLS:
            PARAM_NAMES_USE_GLS.append(i)

    for i in PARAM_NAMES_CONSIDER_MS:
        if i not in PARAM_REMOVE_MS:
            PARAM_NAMES_USE_MS.append(i)
    return PARAM_NAMES_USE_GLS, PARAM_NAMES_USE_MS


""" Other parameter assignments that are not used anymore """
# # PARAM_NAMES_USE_GLS=[]
# PARAM_NAMES_USE_MS=[]

# PARAM_NAMES_USE_GLS=['aG', 'sG', 'Gtail', 'aR', 'Rtail', 'sY', 'Ytail', 'tMove']+\
#                     ['fracR','fracG']+\
#                     ['CoM0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['CoM1R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes']]+\
#                     ['CoM0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['CoM1Y{0}_GLS'.format(par) for par in ['avgRes']]+\
#                     ['MoI0G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes']]+\
#                     ['MoI1G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes']]+\
#                     ['MoI0R{0}_GLS'.format(par) for par in ['td', 'avgRes', 'stdTail']]+\
#                     ['MoI1R{0}_GLS'.format(par) for par in ['td', 'avgRes']]+\
#                     ['MoI0Y{0}_GLS'.format(par) for par in ['td', 'scale']]+\
#                     ['MoI1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes']]#+\
#
# PARAM_NAMES_USE_MS=['maxG', 'aG', 'mG',\
#                     'maxR', 'aR', 'mR','sR',\
#                     'tMove','mRmG', 'tScale']+\
#                    ['aSigHead', 'mSigHead', 'sSigHead', 'mGaussHead', 'maxHead', 'tailHead', 'devHead']+\
#                    ['scaleLength', 'devLength', 'tailLength']+\
#                    ['CoM0G{0}_MS'.format(par) for par in ['scale', 'avgRes', 'stdTail']]+\
#                    ['CoM1G{0}_MS'.format(par) for par in ['td', 'scale', 'stdTail']]+\
#                    ['CoM0R{0}_MS'.format(par) for par in ['td', 'scale', 'stdTail']]+\
#                    ['CoM1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes']]+\
#                    ['MoI0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes']]#+\
#
# # '''removed params after optimization 7/20/17'''
# PARAM_NAMES_USE_GLS=['aG', 'sG', 'Gtail', 'aR', 'sR', 'Rtail', 'sY', 'tMove']+\
#                     ['fracR','fracY']+\
#                     ['CoM0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes']]+\
#                     ['CoM1R{0}_GLS'.format(par) for par in ['scale']]+\
#                     ['CoM0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes']]+\
#                     ['CoM1Y{0}_GLS'.format(par) for par in ['scale', 'avgRes']]+\
#                     ['MoI0G{0}_GLS'.format(par) for par in ['scale', 'avgRes']]+\
#                     ['MoI1G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes']]+\
#                     ['MoI0R{0}_GLS'.format(par) for par in ['td', 'avgRes']]+\
#                     ['MoI1R{0}_GLS'.format(par) for par in ['td', 'scale']]+\
#                     ['MoI0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes']]+\
#                     ['MoI1Y{0}_GLS'.format(par) for par in ['td', 'avgRes']]#+\
# #                     ['{0}eval_{t}_{s}'.format(cn,t=t,s='GLS') for t in TIME_POINTS_EVAL_GLS for cn in ['CoM{0}{1}'.format(ax, c) for c in ['R', 'Y']]]+\
# #                     ['{0}eval_{t}_{s}'.format(cn,t=t,s='GLS') for t in TIME_POINTS_EVAL_GLS for cn in ['MoI{0}{1}'.format(ax, c) for c in ['G', 'R', 'Y']]]
#
# PARAM_NAMES_USE_MS=['maxG', 'aG', 'sG',\
#                     'maxR', 'aR', 'sR',\
#                     'tMove','mRmG', 'tScale']+\
#                    ['aSigHead', 'mSigHead', 'sSigHead', 'aGaussHead', 'mGaussHead', 'sGaussHead', 'maxHead', 'tailHead', 'devHead']+\
#                    ['scaleLength', 'devLength', 'tailLength']+\
#                    ['CoM0G{0}_MS'.format(par) for par in ['scale', 'avgRes', 'stdTail']]+\
#                    ['CoM1G{0}_MS'.format(par) for par in ['scale', 'avgRes']]+\
#                    ['CoM0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM1R{0}_MS'.format(par) for par in ['td', 'avgRes', 'stdTail']]+\
#                    ['MoI0G{0}_MS'.format(par) for par in ['scale', 'stdTail']]+\
#                    ['MoI1G{0}_MS'.format(par) for par in ['scale', 'avgRes', 'stdTail']]+\
#                    ['MoI0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]#+\
# #                    ['{0}eval_{t}_{s}'.format(cn,t=t,s='MS') for t in TIME_POINTS_EVAL_MS for cn in ['CoM{0}{1}'.format(ax, c) for c in ['G', 'R']]]+\
# #                    ['{0}eval_{t}_{s}'.format(cn,t=t,s='MS') for t in TIME_POINTS_EVAL_MS for cn in ['MoI{0}{1}'.format(ax, c) for c in ['G', 'R']]]

'''removed params after PAD optimization 7/22/17'''
# PARAM_NAMES_USE_GLS=['aG', 'Gtail', 'aR', 'sR', 'Rtail', 'sY', 'Ytail', 'mRmG', 'tMove']+\
#                     ['fracR','fracG','fracY', 'tScale']+\
#                     ['CoM0R{0}_GLS'.format(par) for par in ['td', 'avgRes', 'stdTail']]+\
#                     ['CoM1R{0}_GLS'.format(par) for par in ['scale', 'stdTail']]+\
#                     ['CoM0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes']]+\
#                     ['CoM1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes']]+\
#                     ['MoI0G{0}_GLS'.format(par) for par in ['scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1R{0}_GLS'.format(par) for par in ['avgRes', 'stdTail']]+\
#                     ['MoI0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1Y{0}_GLS'.format(par) for par in ['td', 'avgRes', 'stdTail']]#+\
#
# PARAM_NAMES_USE_MS=['maxG', 'aG', 'mG','sG',\
#                     'maxR', 'aR', 'sR',\
#                     'mRmG', 'tScale']+\
#                    ['maxHead', 'tailHead', 'scaleHead', 'devHead']+\
#                    ['scaleLength', 'devLength', 'tailLength']+\
#                    ['CoM0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM0R{0}_MS'.format(par) for par in ['td', 'scale', 'stdTail']]+\
#                    ['CoM1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI1G{0}_MS'.format(par) for par in ['td', 'scale']]+\
#                    ['MoI0R{0}_MS'.format(par) for par in ['td', 'scale', 'stdTail']]+\
#                    ['MoI1R{0}_MS'.format(par) for par in ['scale', 'avgRes', 'stdTail']]#+\

# '''removed params after PAD optimization 7/27/17'''
# PARAM_NAMES_USE_GLS=['aG', 'sG', 'Gtail', 'aR', 'sR', 'Rtail', 'aY', 'sY', 'Ytail', 'mRmG', 'tMove']+\
#                     ['fracR','fracG','fracY', 'tScale']+\
#                     ['CoM0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['CoM1R{0}_GLS'.format(par) for par in ['td', 'scale', 'stdTail']]+\
#                     ['CoM0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['CoM1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'stdTail']]+\
#                     ['MoI0G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1R{0}_GLS'.format(par) for par in ['td', 'avgRes', 'stdTail']]+\
#                     ['MoI0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1Y{0}_GLS'.format(par) for par in ['avgRes', 'stdTail']]#+\
#
# PARAM_NAMES_USE_MS=['maxG', 'aG', 'mG','sG',\
#                     'maxR', 'aR', 'mR','sR',\
#                     'tMove','mRmG', 'tScale']+\
#                    ['maxHead', 'tailHead', 'scaleHead', 'devHead']+\
#                    ['scaleLength', 'devLength', 'tailLength']+\
#                    ['CoM0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI1R{0}_MS'.format(par) for par in ['scale', 'avgRes', 'stdTail']]

'''removed params after PAD optimization by maximizing gap 9/29/17'''
# PARAM_NAMES_USE_GLS=['aG', 'Gtail', 'aR', 'sY', 'mRmG', 'mYmG', 'tMove']+\
#                     ['fracR','fracG','fracY', 'tScale']+\
#                     ['CoM0R{0}_GLS'.format(par) for par in ['td']]+\
#                     ['CoM1R{0}_GLS'.format(par) for par in ['scale']]+\
#                     ['CoM0Y{0}_GLS'.format(par) for par in ['scale']]+\
#                     ['CoM1Y{0}_GLS'.format(par) for par in ['td']]+\
#                     ['MoI0G{0}_GLS'.format(par) for par in ['scale', 'stdTail']]+\
#                     ['MoI0R{0}_GLS'.format(par) for par in ['td', 'avgRes', 'stdTail']]+\
#                     ['MoI0Y{0}_GLS'.format(par) for par in ['td', 'stdTail']]+\
#                     ['MoI1Y{0}_GLS'.format(par) for par in ['avgRes', 'stdTail']]
#
# PARAM_NAMES_USE_MS=['maxG', 'aG', 'mG','sG',\
#                     'maxR', 'aR', 'mR',\
#                     'mRmG', 'tScale']+\
#                    ['maxHead', 'tailHead', 'scaleHead', 'devHead']+\
#                    ['devLength', 'tailLength']+\
#                    ['CoM0G{0}_MS'.format(par) for par in ['td', 'scale']]+\
#                    ['CoM1G{0}_MS'.format(par) for par in ['td']]+\
#                    ['CoM0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM1R{0}_MS'.format(par) for par in ['td', 'scale']]+\
#                    ['MoI0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes']]+\
#                    ['MoI1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes']]+\
#                    ['MoI0R{0}_MS'.format(par) for par in ['stdTail']]+\
#                    ['MoI1R{0}_MS'.format(par) for par in ['scale', 'avgRes', 'stdTail']]

'''removed params with PCA<0.6 9/29/17'''
# PARAM_NAMES_USE_GLS = ['maxR', 'maxG', 'aG', 'Gtail', 'aR', 'sR', 'Rtail', 'aY', 'sY', 'Ytail', 'mRmG', 'mYmG',
#                        'tMove'] + \
#                       ['fracR', 'fracG', 'fracY', 'tScale'] + \
#                       ['CoM0R{0}_GLS'.format(par) for par in ['td']] + \
#                       ['CoM1R{0}_GLS'.format(par) for par in ['td']] + \
#                       ['CoM0Y{0}_GLS'.format(par) for par in ['td', 'stdTail']] + \
#                       ['CoM1Y{0}_GLS'.format(par) for par in ['td']] + \
#                       ['MoI0G{0}_GLS'.format(par) for par in ['td', 'avgRes', 'stdTail']] + \
#                       ['MoI1G{0}_GLS'.format(par) for par in ['td', 'avgRes', 'stdTail']] + \
#                       ['MoI0R{0}_GLS'.format(par) for par in ['td', 'stdTail']] + \
#                       ['MoI1R{0}_GLS'.format(par) for par in ['td', 'avgRes', 'stdTail']] + \
#                       ['MoI0Y{0}_GLS'.format(par) for par in ['td', 'avgRes', 'stdTail']] + \
#                       ['MoI1Y{0}_GLS'.format(par) for par in ['td', 'avgRes', 'stdTail']]  # +\
# #                     ['{0}eval_{t}_{s}'.format(cn,t=t,s='GLS') for t in TIME_POINTS_EVAL_GLS for cn in ['CoM{0}{1}'.format(ax, c) for c in ['R', 'Y']]]+\
# #                     ['{0}eval_{t}_{s}'.format(cn,t=t,s='GLS') for t in TIME_POINTS_EVAL_GLS for cn in ['MoI{0}{1}'.format(ax, c) for c in ['G', 'R', 'Y']]]
#
# PARAM_NAMES_USE_MS = ['maxG', 'aG', \
#                       'maxR', 'aR', 'mR', \
#                        'mRmG'] + \
#                      ['maxHead', 'tailHead', 'scaleHead', 'devHead'] + \
#                      ['devLength', 'tailLength'] + \
#                      ['CoM0G{0}_MS'.format(par) for par in ['stdTail']] + \
#                      ['CoM1G{0}_MS'.format(par) for par in ['stdTail']] + \
#                      ['CoM0R{0}_MS'.format(par) for par in ['stdTail']] + \
#                      ['CoM1R{0}_MS'.format(par) for par in ['stdTail']] + \
#                      ['MoI0G{0}_MS'.format(par) for par in ['td', 'scale', 'stdTail']] + \
#                      ['MoI1G{0}_MS'.format(par) for par in ['td', 'scale', 'stdTail']] + \
#                      ['MoI0R{0}_MS'.format(par) for par in ['stdTail']] + \
#                      ['MoI1R{0}_MS'.format(par) for par in ['stdTail']]  # +\
# #                    ['aSigHead', 'mSigHead', 'sSigHead', 'aGaussHead', 'mGaussHead', 'sGaussHead']+\
# #                    ['{0}eval_{t}_{s}'.format(cn,t=t,s='MS') for t in TIME_POINTS_EVAL_MS for cn in ['CoM{0}{1}'.format(ax, c) for c in ['G', 'R']]]+\
# #                    ['{0}eval_{t}_{s}'.format(cn,t=t,s='MS') for t in TIME_POINTS_EVAL_MS for cn in ['MoI{0}{1}'.format(ax, c) for c in ['G', 'R']]]

'''PCA RNAi with th cutoff'''
# th = 0.1 ###NOTE DIRECTION OF > < in PARAMS USE LINES
# ordered_names = ['sG_GLS', 'CoM1Rscale_GLS_GLS', 'CoM0RavgRes_MS_MS', 'mR_MS', 'MoI1Yscale_GLS_GLS', 'CoM0Yscale_GLS_GLS', 'CoM1RavgRes_MS_MS', 'CoM0YavgRes_GLS_GLS', 'MoI0RavgRes_MS_MS', 'MoI1Gscale_GLS_GLS', 'sR_GLS', 'CoM1Yscale_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'mRmG_MS', 'sY_GLS', 'CoM1RavgRes_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'CoM0GavgRes_MS_MS', 'CoM1GavgRes_MS_MS', 'CoM1Gscale_MS_MS', 'fracY_GLS', 'mRmG_GLS', 'MoI0Gscale_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'tMove_GLS', 'CoM0GstdTail_MS_MS', 'CoM0RavgRes_GLS_GLS', 'CoM0Gscale_MS_MS', 'MoI1RstdTail_GLS_GLS', 'MoI0GstdTail_MS_MS', 'mYmG_GLS', 'scaleHead_MS', 'MoI0RstdTail_GLS_GLS', 'MoI0RavgRes_GLS_GLS', 'MoI0YstdTail_GLS_GLS', 'CoM0Rscale_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'fracR_GLS', 'sR_MS', 'CoM0Rscale_MS_MS', 'MoI1Rscale_MS_MS', 'CoM0Rtd_MS_MS', 'MoI1GstdTail_MS_MS', 'aR_MS', 'MoI0GavgRes_GLS_GLS', 'maxR_MS', 'MoI0Rscale_MS_MS', 'CoM1GstdTail_MS_MS', 'tailHead_MS', 'MoI1Rtd_MS_MS', 'MoI1GavgRes_GLS_GLS', 'MoI0GavgRes_MS_MS', 'MoI0GstdTail_GLS_GLS', 'devHead_MS', 'CoM1Rtd_MS_MS', 'maxG_MS', 'MoI1GavgRes_MS_MS', 'maxHead_MS', 'MoI0Gscale_MS_MS', 'MoI1RavgRes_MS_MS', 'CoM1YstdTail_GLS_GLS', 'sG_MS', 'CoM0Gtd_MS_MS', 'mG_MS', 'CoM1Rscale_MS_MS', 'MoI1GstdTail_GLS_GLS', 'aG_MS', 'MoI1YstdTail_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'CoM0RstdTail_GLS_GLS', 'Gtail_GLS', 'fracG_GLS', 'tScale_GLS', 'MoI1YavgRes_GLS_GLS', 'aY_GLS', 'devLength_MS', 'CoM0RstdTail_MS_MS', 'MoI0RstdTail_MS_MS', 'MoI0YavgRes_GLS_GLS', 'MoI1RstdTail_MS_MS', 'maxG_GLS', 'MoI1Gscale_MS_MS', 'MoI1RavgRes_GLS_GLS', 'CoM1Gtd_MS_MS', 'aR_GLS', 'CoM1RstdTail_MS_MS', 'MoI0Rtd_MS_MS', 'aG_GLS', 'MoI1Gtd_GLS_GLS', 'Rtail_GLS', 'MoI0Gtd_GLS_GLS', 'Ytail_GLS', 'tScale_MS', 'CoM1Rtd_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'CoM0Ytd_GLS_GLS', 'MoI0Gtd_MS_MS', 'maxR_GLS', 'CoM0Rtd_GLS_GLS', 'MoI1Gtd_MS_MS', 'MoI1Ytd_GLS_GLS', 'MoI1Rtd_GLS_GLS', 'MoI0Rtd_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'tailLength_MS', 'scaleLength_MS', 'tMove_MS']
# ordered_vals = [0.012039431709191618, 0.019749992937013246, 0.020108365850790962, 0.021445662170860625, 0.02448890301163046, 0.02952734079977215, 0.029946747808872687, 0.030608497592490375, 0.033313679171588007, 0.033580415241705318, 0.034105903218143691, 0.039643709215971171, 0.039676821846647461, 0.045452041459597699, 0.045731480837153608, 0.052062039095040646, 0.053822701831792476, 0.055511433751136763, 0.060602089064392885, 0.060852504358968829, 0.061159036078833225, 0.063774176534408145, 0.068443308003993661, 0.069718255176573929, 0.071224014222644322, 0.072924629976550803, 0.076272078268297405, 0.081360753460642338, 0.081729251111061968, 0.082888078621574726, 0.086973928620760504, 0.090939294280949795, 0.091417412978738316, 0.092018966417296288, 0.095539017314471231, 0.095608254156653666, 0.10076505481750203, 0.10403976558078926, 0.10481743801595673, 0.10532702400796974, 0.10743920552659991, 0.11412759204753554, 0.11539046561514059, 0.11682207636658437, 0.116860524001627, 0.11840694547995549, 0.1187907595752557, 0.12209170847728917, 0.12569384814040768, 0.12582079108985905, 0.12647734641995556, 0.1268497163538348, 0.12737815929861024, 0.1277137372445758, 0.12859914799221583, 0.12920194952058728, 0.12977094388811858, 0.1305809466726964, 0.13146854567826274, 0.13166214861697534, 0.13176327055141246, 0.13405203570141824, 0.1355984255957107, 0.13584466254373809, 0.13775504685360584, 0.13787475836293545, 0.1383053985227761, 0.14164118977343856, 0.14277610911777069, 0.15038586740367244, 0.15117872729031123, 0.1515653894690363, 0.159447089567366, 0.16397766382778306, 0.16426422782387967, 0.16452562130731432, 0.16901769206558037, 0.17022634054061236, 0.17337774784814269, 0.17588121510993837, 0.17898398065044629, 0.18231293309663274, 0.18650791011836471, 0.18726165954149632, 0.19103836190950058, 0.19772481358577876, 0.20059699878199821, 0.2014396200676406, 0.20203097070241183, 0.20740700358178332, 0.21106898594566087, 0.2125216315916848, 0.21292674288070798, 0.21297339324097447, 0.21321949733893109, 0.21570795127096104, 0.21598114007959329, 0.21831949848908636, 0.22475932256369788, 0.22755557951911698, 0.24436056901433673, 0.24803173704580792, 0.25230490585863896, 0.25316449878777203, 0.26182391460504539, 0.2908422949580407, 0.30996172020805718, 1.1792152955686281]
#
# gls_ordered_names = []
# gls_ordered_vals = []
# ms_ordered_names = []
# ms_ordered_vals = []
#
# for i in range(len(ordered_names)):
#     name = ordered_names[i]
#     if name[-2:] == 'LS':
#         gls_ordered_names.append(name[:-4])
#         gls_ordered_vals.append(ordered_vals[i])
#     else:
#         ms_ordered_names.append(name[:-3])
#         ms_ordered_vals.append(ordered_vals[i])
#
# PARAM_NAMES_USE_GLS = list(np.array(gls_ordered_names)[np.array(gls_ordered_vals) > th])
# PARAM_NAMES_USE_MS = list(np.array(ms_ordered_names)[np.array(ms_ordered_vals) > th])


'''LDA Control VS RNAi with th cutoff'''
# th = 0.06 ###NOTE DIRECTION OF > < in PARAMS USE LINES
#
# gls_ordered_names = ['fracG', 'aG', 'CoM0Rtd_GLS', 'Rtail', 'MoI1RstdTail_GLS', 'MoI0RstdTail_GLS', 'MoI0YstdTail_GLS', 'aR', 'tScale', 'MoI0YavgRes_GLS', 'MoI1Gtd_GLS', 'MoI0GstdTail_GLS', 'Ytail', 'MoI0Gtd_GLS', 'sG', 'MoI0Ytd_GLS', 'maxG', 'MoI1RavgRes_GLS', 'MoI0Yscale_GLS', 'MoI1GstdTail_GLS', 'CoM0RavgRes_GLS', 'tMove', 'fracR', 'CoM0RstdTail_GLS', 'MoI0Gscale_GLS', 'aY', 'MoI1YstdTail_GLS', 'fracY', 'CoM1YstdTail_GLS', 'CoM0Ytd_GLS', 'maxR', 'CoM1Rscale_GLS', 'MoI0RavgRes_GLS', 'CoM1Ytd_GLS', 'MoI1Gscale_GLS', 'mYmG', 'MoI0Rscale_GLS', 'MoI1Ytd_GLS', 'MoI1GavgRes_GLS', 'CoM1RavgRes_GLS', 'CoM0YavgRes_GLS', 'MoI1Rtd_GLS', 'MoI0GavgRes_GLS', 'CoM1YavgRes_GLS', 'CoM1Yscale_GLS', 'CoM1RstdTail_GLS', 'CoM0YstdTail_GLS', 'mRmG', 'MoI1Yscale_GLS', 'MoI1YavgRes_GLS', 'CoM0Yscale_GLS', 'sR', 'MoI1Rscale_GLS', 'CoM1Rtd_GLS', 'Gtail', 'CoM0Rscale_GLS', 'sY', 'MoI0Rtd_GLS']
# gls_ordered_vals = [0.48113002183964271, 0.31372124260847484, 0.26800419175383322, 0.26554057092918593, 0.25730235786651634, 0.25035349572662707, 0.24504077143726233, 0.17353345968960446, 0.1612681468539581, 0.15533381767560903, 0.15351343469514078, 0.13364008605135758, 0.12563629488868144, 0.11964059389615248, 0.11273165627372082, 0.10978158752233884, 0.10821964309159536, 0.10556910732668535, 0.098714886409082478, 0.096497504884584598, 0.09421896801641734, 0.093383087318861954, 0.091175411624697306, 0.089363096786373539, 0.088657776491762777, 0.088012244764361813, 0.081218780457695122, 0.078281900357824938, 0.076774060867646446, 0.07436061801682381, 0.06335585705887832, 0.063222037719830732, 0.061551545138731914, 0.057875712913133676, 0.057458737666241953, 0.057147956378518223, 0.055533933476098626, 0.049526459635091315, 0.04652357110346058, 0.045674979633002981, 0.044315771286400359, 0.043912778686117288, 0.043447899420101681, 0.043335579737369394, 0.041786506361917157, 0.040798113476917211, 0.038223197190260252, 0.029811718346664809, 0.018860258172206875, 0.017659509012673248, 0.016515273493303777, 0.014745469510940485, 0.014114827161759444, 0.013229058964782849, 0.012357322444825777, 0.0097696807924913028, 0.0051036510569498691, 0.0012317901648923198]
# PARAM_NAMES_USE_GLS = list(np.array(gls_ordered_names)[np.array(gls_ordered_vals) > th])
#
# ms_ordered_names = ['CoM0Gtd_MS', 'aG', 'devLength', 'MoI0Rtd_MS', 'CoM1RstdTail_MS', 'CoM0RstdTail_MS', 'MoI1RstdTail_MS', 'tailLength', 'scaleLength', 'maxHead', 'CoM1Gtd_MS', 'CoM0RavgRes_MS', 'MoI0RstdTail_MS', 'CoM1Rscale_MS', 'sG', 'MoI0Gtd_MS', 'MoI0GstdTail_MS', 'devHead', 'CoM0Gscale_MS', 'tScale', 'mR', 'sR', 'MoI0Gscale_MS', 'MoI1GstdTail_MS', 'MoI1Rscale_MS', 'tailHead', 'MoI0Rscale_MS', 'MoI0RavgRes_MS', 'CoM0Rscale_MS', 'maxR', 'CoM0Rtd_MS', 'maxG', 'scaleHead', 'CoM1Gscale_MS', 'CoM0GavgRes_MS', 'aR', 'MoI0GavgRes_MS', 'CoM1GstdTail_MS', 'MoI1GavgRes_MS', 'MoI1RavgRes_MS', 'MoI1Gtd_MS', 'MoI1Gscale_MS', 'mG', 'CoM1GavgRes_MS', 'CoM0GstdTail_MS', 'tMove', 'MoI1Rtd_MS', 'CoM1RavgRes_MS', 'mRmG', 'CoM1Rtd_MS']
# ms_ordered_vals = [0.35915290423651047, 0.30494068839544852, 0.30066021626832756, 0.24742876078895712, 0.24504258494698347, 0.23849751509205433, 0.2134108669003392, 0.1877179183520179, 0.18684426129763551, 0.17650225737604075, 0.17258792835616468, 0.16809444273231944, 0.15725927866881886, 0.15528715130244056, 0.15001482902556074, 0.14663407575531426, 0.12890216017887893, 0.12250791975025853, 0.1223151553408722, 0.12049850954779555, 0.11611334382274001, 0.113130480072474, 0.11222520033232142, 0.10232842049242719, 0.098575627415847925, 0.097044424105865112, 0.095212416895697621, 0.09501615918097879, 0.093741427327345583, 0.091304431672518477, 0.085960115020098107, 0.075939926083960302, 0.071447226657642762, 0.064667678886869501, 0.064238055208283976, 0.063810679671787746, 0.058261287317130442, 0.04572364323812636, 0.045208403989114604, 0.044763732638655677, 0.039440590821692073, 0.034836530689212483, 0.033824316233025301, 0.031325815197525163, 0.02999163253306544, 0.022709592117525039, 0.014986664416965174, 0.013737328214082693, 0.004596950717461387, 0.0039588301925902007]
# PARAM_NAMES_USE_MS = list(np.array(ms_ordered_names)[np.array(ms_ordered_vals) > th])


'''LDA group separation with th cutoff'''
# th = 0.152 ###NOTE DIRECTION OF > < in PARAMS USE LINES
#
# ordered_vals = [0.24881804484742129, 0.19252744183340248, 0.15716907685463088, 0.15487762049743811, 0.15460749149196362, 0.15151019175426497, 0.14264678937368394, 0.13185823862391741, 0.11956463271069107, 0.11640593899031337, 0.11405210828491567, 0.10929527591756413, 0.10469005473918411, 0.10372296256504947, 0.1030019763263773, 0.10095859611157891, 0.096603883188199097, 0.095789318599585421, 0.094323765743408325, 0.093655817993354526, 0.091529757344144022, 0.087159806542470492, 0.084308649231702809, 0.081600018601915847, 0.08153669054356795, 0.081032826481230294, 0.080575374484341045, 0.078428472349276454, 0.07736784481205701, 0.077336339504827475, 0.077312001592042393, 0.073607746637832455, 0.073499584057848566, 0.071474176516757651, 0.070913804591033525, 0.070694301603590845, 0.067890073644661986, 0.067408725060521313, 0.067179861185292039, 0.066642836658361959, 0.066295832259342991, 0.065769908219571807, 0.063375508740769484, 0.062708184349171872, 0.060417056926668856, 0.059618208998629907, 0.059574408989533825, 0.059364450982678882, 0.057918590249701134, 0.057022373328694552, 0.056360284349917758, 0.056273512617820404, 0.056183279133283207, 0.055770781982951108, 0.055256112415800661, 0.054958735835662458, 0.054769777783566113, 0.051027287539727559, 0.049381332047579059, 0.04866233667540959, 0.047890560621973534, 0.047671866344498354, 0.047359036635925236, 0.047140352192845861, 0.044196857882386621, 0.043244489842710343, 0.042444537780948498, 0.042370075864284269, 0.041960360296385307, 0.041033758492069999, 0.039355824397888171, 0.039339975143661055, 0.037614340958251613, 0.037241872681711329, 0.036433902834751071, 0.036024234798450974, 0.035230114811389421, 0.035084597050199985, 0.035022893474141018, 0.034936119028084432, 0.03389320301293814, 0.033508783216108201, 0.032926087563099803, 0.032077088384365404, 0.031921328456955569, 0.030894552880552362, 0.030691718770149267, 0.029797279699852845, 0.029277698153545376, 0.028676874316371399, 0.028218046808675304, 0.027832121651710597, 0.026951294921671348, 0.026051005679674985, 0.025425069687622329, 0.024932192527758079, 0.024865618115877602, 0.024604301936114292, 0.023581690315975276, 0.023011542944812908, 0.022013068839167114, 0.0194155004583456, 0.019356127788426145, 0.018898123829187589, 0.017430621115827032, 0.017322068484622647, 0.010186593599429746, 0.0071103246606383117]
# ordered_names = ['Rtail_GLS', 'CoM1RstdTail_MS_MS', 'MoI1RstdTail_MS_MS', 'CoM0RstdTail_MS_MS', 'devHead_MS', 'CoM0GstdTail_MS_MS', 'MoI0RstdTail_GLS_GLS', 'MoI1GstdTail_GLS_GLS', 'CoM1Ytd_GLS_GLS', 'CoM0YstdTail_GLS_GLS', 'aR_GLS', 'MoI0GstdTail_MS_MS', 'fracR_GLS', 'MoI0YavgRes_GLS_GLS', 'MoI0GavgRes_MS_MS', 'Ytail_GLS', 'CoM0RavgRes_MS_MS', 'MoI0GavgRes_GLS_GLS', 'aG_MS', 'CoM0Ytd_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'maxR_MS', 'aG_GLS', 'MoI0Ytd_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'MoI0GstdTail_GLS_GLS', 'CoM0Rtd_MS_MS', 'MoI1Rtd_MS_MS', 'mRmG_MS', 'MoI0RstdTail_MS_MS', 'CoM1GstdTail_MS_MS', 'MoI1YstdTail_GLS_GLS', 'MoI1GavgRes_MS_MS', 'CoM1Gscale_MS_MS', 'MoI1Gscale_MS_MS', 'Gtail_GLS', 'CoM1Yscale_GLS_GLS', 'tailHead_MS', 'MoI1RavgRes_GLS_GLS', 'MoI1RstdTail_GLS_GLS', 'CoM0Gscale_MS_MS', 'sG_GLS', 'CoM1RavgRes_MS_MS', 'CoM0Rscale_GLS_GLS', 'maxG_GLS', 'maxHead_MS', 'CoM1YstdTail_GLS_GLS', 'MoI1Gscale_GLS_GLS', 'CoM1Rtd_MS_MS', 'maxR_GLS', 'maxG_MS', 'MoI1GavgRes_GLS_GLS', 'tMove_GLS', 'CoM1GavgRes_MS_MS', 'MoI0Gtd_MS_MS', 'MoI1Ytd_GLS_GLS', 'MoI1Gtd_GLS_GLS', 'scaleLength_MS', 'CoM0Gtd_MS_MS', 'sR_GLS', 'MoI0Rscale_GLS_GLS', 'CoM0Rscale_MS_MS', 'MoI1Rtd_GLS_GLS', 'aR_MS', 'MoI1GstdTail_MS_MS', 'CoM1RstdTail_GLS_GLS', 'CoM1Rscale_MS_MS', 'CoM0RavgRes_GLS_GLS', 'devLength_MS', 'sY_GLS', 'MoI0Yscale_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'CoM0Rtd_GLS_GLS', 'CoM1Gtd_MS_MS', 'CoM1Rscale_GLS_GLS', 'CoM0GavgRes_MS_MS', 'MoI0Rtd_MS_MS', 'CoM0Yscale_GLS_GLS', 'MoI0Rtd_GLS_GLS', 'MoI1YavgRes_GLS_GLS', 'CoM0RstdTail_GLS_GLS', 'MoI0Gscale_MS_MS', 'fracG_GLS', 'scaleHead_MS', 'CoM1RavgRes_GLS_GLS', 'MoI0RavgRes_MS_MS', 'tScale_GLS', 'tScale_MS', 'fracY_GLS', 'MoI1RavgRes_MS_MS', 'MoI0YstdTail_GLS_GLS', 'mR_MS', 'sG_MS', 'MoI0Gscale_GLS_GLS', 'CoM1YavgRes_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'tailLength_MS', 'aY_GLS', 'MoI0Rscale_MS_MS', 'mG_MS', 'mRmG_GLS', 'MoI1Rscale_MS_MS', 'MoI0RavgRes_GLS_GLS', 'sR_MS', 'MoI1Gtd_MS_MS', 'MoI1Yscale_GLS_GLS', 'tMove_MS', 'mYmG_GLS']
#
# gls_ordered_names = []
# gls_ordered_vals = []
# ms_ordered_names = []
# ms_ordered_vals = []
#
# for i in range(len(ordered_names)):
#     name = ordered_names[i]
#     if name[-2:] == 'LS':
#         gls_ordered_names.append(name[:-4])
#         gls_ordered_vals.append(ordered_vals[i])
#     else:
#         ms_ordered_names.append(name[:-3])
#         ms_ordered_vals.append(ordered_vals[i])
#
# PARAM_NAMES_USE_GLS = list(np.array(gls_ordered_names)[np.array(gls_ordered_vals) > th])
# PARAM_NAMES_USE_MS = list(np.array(ms_ordered_names)[np.array(ms_ordered_vals) > th])

'''removed params after PAD rank optimization [bug + no rank normalization] 10/12/17'''
# PARAM_NAMES_USE_GLS=['maxG', 'aG', 'Gtail', 'aR', 'Rtail', 'aY', 'sY', 'tMove']+\
#                     ['fracR']+\
#                     ['CoM1R{0}_GLS'.format(par) for par in ['scale']]+\
#                     ['CoM0Y{0}_GLS'.format(par) for par in ['stdTail']]+\
#                     ['CoM1Y{0}_GLS'.format(par) for par in ['td', 'avgRes']]+\
#                     ['MoI0G{0}_GLS'.format(par) for par in ['avgRes']]+\
#                     ['MoI1G{0}_GLS'.format(par) for par in ['scale']]+\
#                     ['MoI0R{0}_GLS'.format(par) for par in ['td', 'avgRes']]+\
#                     ['MoI0Y{0}_GLS'.format(par) for par in ['scale']]
#
# PARAM_NAMES_USE_MS=['sG',\
#                     'maxR', 'aR', 'mR','sR',\
#                     'mRmG']+\
#                    ['maxHead', 'tailHead', 'scaleHead', 'devHead']+\
#                    ['devLength', 'tailLength']+\
#                    ['CoM0G{0}_MS'.format(par) for par in ['td', 'scale']]+\
#                    ['CoM1G{0}_MS'.format(par) for par in ['td', 'avgRes', 'stdTail']]+\
#                    ['CoM0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM1R{0}_MS'.format(par) for par in ['stdTail']]+\
#                    ['MoI0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI0R{0}_MS'.format(par) for par in ['scale', 'avgRes', 'stdTail']]+\
#                    ['MoI1R{0}_MS'.format(par) for par in ['td', 'avgRes', 'stdTail']]

'''removed params after PAD rank optimization 10/16/17'''
# PARAM_NAMES_USE_GLS=['maxR', 'maxG', 'aG', 'Gtail', 'aR', 'Rtail', 'aY', 'Ytail', 'mYmG', 'tMove']+\
#                     ['fracR','fracY', 'tScale']+\
#                     ['CoM0R{0}_GLS'.format(par) for par in ['avgRes']]+\
#                     ['CoM1R{0}_GLS'.format(par) for par in ['scale', 'avgRes']]+\
#                     ['CoM1Y{0}_GLS'.format(par) for par in ['scale', 'avgRes']]+\
#                     ['MoI0G{0}_GLS'.format(par) for par in ['scale', 'avgRes']]+\
#                     ['MoI1G{0}_GLS'.format(par) for par in ['stdTail']]+\
#                     ['MoI0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1R{0}_GLS'.format(par) for par in ['td']]+\
#                     ['MoI0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1Y{0}_GLS'.format(par) for par in ['scale']]
#
# PARAM_NAMES_USE_MS=['maxG', 'sG',\
#                     'maxR', 'aR','sR']+\
#                    ['maxHead', 'tailHead', 'devHead']+\
#                    ['tailLength']+\
#                    ['CoM0G{0}_MS'.format(par) for par in ['scale']]+\
#                    ['CoM1G{0}_MS'.format(par) for par in ['scale', 'avgRes']]+\
#                    ['CoM0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM1R{0}_MS'.format(par) for par in ['avgRes', 'stdTail']]+\
#                    ['MoI0G{0}_MS'.format(par) for par in ['td', 'scale', 'stdTail']]+\
#                    ['MoI1G{0}_MS'.format(par) for par in ['scale', 'stdTail']]+\
#                    ['MoI0R{0}_MS'.format(par) for par in ['scale', 'avgRes', 'stdTail']]+\
#                    ['MoI1R{0}_MS'.format(par) for par in ['td', 'avgRes', 'stdTail']]

'''group rank optimization with th cutoff 10/24/17'''
# th = 2 ###NOTE DIRECTION OF > < in PARAMS USE LINES
#
# ordered_vals = [12,12,12,12,12,12,12,12,12,12,12,12,12,12,11,11,11,11,11,11,10,10,10,9,9,9,8,8,8,8,8,8,8,7,7,7,7,7,7,7,7,7,7,6,6,6,6,6,6,5,5,4,4,4,4,4,4,4,3,3,3,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1]
# ordered_names = ['scaleLength_MS', 'mG_MS', 'tScale_MS', 'mRmG_MS', 'CoM0Ytd_GLS_GLS', 'MoI1Rscale_GLS_GLS', 'CoM0GavgRes_MS_MS', 'CoM1YstdTail_GLS_GLS', 'MoI0GavgRes_MS_MS', 'CoM0Gtd_MS_MS', 'MoI1Gtd_GLS_GLS', 'CoM0YavgRes_GLS_GLS', 'tMove_MS', 'MoI1GavgRes_GLS_GLS', 'MoI1Ytd_GLS_GLS', 'MoI1Rscale_MS_MS', 'MoI1RstdTail_GLS_GLS', 'MoI0Gtd_GLS_GLS', 'MoI1Gscale_GLS_GLS', 'mRmG_GLS', 'devLength_MS', 'MoI0GstdTail_GLS_GLS', 'CoM1RstdTail_GLS_GLS', 'CoM1Rtd_GLS_GLS', 'CoM0GstdTail_MS_MS', 'CoM0Rtd_GLS_GLS', 'fracG_GLS', 'CoM1Rtd_MS_MS', 'MoI1YstdTail_GLS_GLS', 'MoI1Gtd_MS_MS', 'fracY_GLS', 'MoI1RavgRes_GLS_GLS', 'MoI0YstdTail_GLS_GLS', 'aG_MS', 'maxG_MS', 'sG_GLS', 'CoM1Ytd_GLS_GLS', 'MoI0Rscale_GLS_GLS', 'CoM1Rscale_MS_MS', 'CoM0YstdTail_GLS_GLS', 'CoM1GstdTail_MS_MS', 'CoM0Rscale_GLS_GLS', 'MoI1GavgRes_MS_MS', 'sR_GLS', 'CoM0RstdTail_GLS_GLS', 'MoI0GstdTail_MS_MS', 'aG_GLS', 'sG_MS', 'CoM0Yscale_GLS_GLS', 'MoI0Rtd_MS_MS', 'MoI1Rtd_GLS_GLS', 'mR_MS', 'CoM1Gtd_MS_MS', 'sY_GLS', 'maxHead_MS', 'MoI0YavgRes_GLS_GLS', 'MoI0Ytd_GLS_GLS', 'Gtail_GLS', 'scaleHead_MS', 'tScale_GLS', 'tailHead_MS', 'MoI0Rscale_MS_MS', 'CoM1Rscale_GLS_GLS', 'maxR_GLS', 'MoI0RavgRes_MS_MS', 'aY_GLS', 'MoI1Gscale_MS_MS', 'maxG_GLS', 'sR_MS', 'MoI1GstdTail_MS_MS', 'CoM1GavgRes_MS_MS', 'MoI0RstdTail_GLS_GLS', 'tailLength_MS', 'CoM0RavgRes_GLS_GLS', 'CoM1RavgRes_MS_MS', 'CoM1Gscale_MS_MS', 'CoM1Yscale_GLS_GLS', 'MoI1YavgRes_GLS_GLS', 'CoM1RavgRes_GLS_GLS', 'MoI0Yscale_GLS_GLS', 'MoI1Yscale_GLS_GLS', ' Ytail_GLS']
# gls_ordered_names = []
# gls_ordered_vals = []
# ms_ordered_names = []
# ms_ordered_vals = []
#
# for i in range(len(ordered_names)):
#     name = ordered_names[i]
#     if name[-2:] == 'LS':
#         gls_ordered_names.append(name[:-4])
#         gls_ordered_vals.append(ordered_vals[i])
#     else:
#         ms_ordered_names.append(name[:-3])
#         ms_ordered_vals.append(ordered_vals[i])
#
# PARAM_REMOVE_GLS = list(np.array(gls_ordered_names)[np.array(gls_ordered_vals) > th])
# PARAM_REMOVE_MS = list(np.array(ms_ordered_names)[np.array(ms_ordered_vals) > th])
# PARAM_NAMES_CONSIDER_GLS=['maxR', 'maxG', 'aG', 'sG', 'Gtail', 'aR', 'sR', 'Rtail', 'aY', 'sY', 'Ytail', 'mRmG', 'mYmG', 'tMove']+\
#                     ['fracR','fracG','fracY', 'tScale']+\
#                     ['CoM0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['CoM1R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['CoM0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['CoM1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI0G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]#+\
#
# PARAM_NAMES_CONSIDER_MS=['maxG', 'aG', 'mG','sG',\
#                     'maxR', 'aR', 'mR','sR',\
#                     'tMove','mRmG', 'tScale']+\
#                    ['maxHead', 'tailHead', 'scaleHead', 'devHead']+\
#                    ['scaleLength', 'devLength', 'tailLength']+\
#                    ['CoM0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]#+\
#
# PARAM_NAMES_USE_GLS = []
# PARAM_NAMES_USE_MS = []
#
# for i in PARAM_NAMES_CONSIDER_GLS:
#     if i not in PARAM_REMOVE_GLS:
#         PARAM_NAMES_USE_GLS.append(i)
#
# for i in PARAM_NAMES_CONSIDER_MS:
#     if i not in PARAM_REMOVE_MS:
#         PARAM_NAMES_USE_MS.append(i)
# # print len(PARAM_NAMES_USE_GLS), len(PARAM_NAMES_USE_MS)
# # print PARAM_NAMES_USE_MS
# # print PARAM_NAMES_USE_GLS

'''group rank optimization with th cutoff 11/3/17- includes dev delay and WT groups'''
# th = 3 ###NOTE DIRECTION OF > < in PARAMS USE LINES
#
# ordered_vals = [14,14,14,14,14,14,14,14,13,13,13,13,13,12,12,12,11,11,11,10,10,10,9,9,9,9,9,9,8,8,7,7,7,6,6,6,6,5,5,5,5,5,5,5,4,4,4,4,4,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1]
# ordered_names = ['MoI1Gtd_GLS_GLS','mG_MS','tScale_MS','mRmG_MS','CoM0Ytd_GLS_GLS','CoM0Gtd_MS_MS','CoM0GavgRes_MS_MS','CoM1YstdTail_GLS_GLS','mRmG_GLS','MoI0GavgRes_MS_MS','tMove_MS','CoM0YavgRes_GLS_GLS','MoI1Rscale_GLS_GLS','MoI1YstdTail_GLS_GLS','MoI1Gscale_GLS_GLS','fracG_GLS','CoM1Rtd_GLS_GLS','MoI1GavgRes_GLS_GLS','MoI1Gtd_MS_MS','CoM0RstdTail_GLS_GLS','MoI0GstdTail_GLS_GLS','CoM0GstdTail_MS_MS','MoI0GstdTail_MS_MS','MoI1RstdTail_GLS_GLS','devLength_MS','scaleLength_MS','CoM0Rtd_GLS_GLS','aG_MS','MoI1Ytd_GLS_GLS','MoI0Rscale_MS_MS','CoM1Gtd_MS_MS','CoM0YstdTail_GLS_GLS','CoM1RstdTail_GLS_GLS','MoI1Rscale_MS_MS','CoM1Ytd_GLS_GLS','CoM1Rscale_MS_MS','CoM0Rscale_GLS_GLS','MoI1RavgRes_GLS_GLS','sR_GLS','tailHead_MS','MoI0YstdTail_GLS_GLS','MoI1GavgRes_MS_MS','sG_MS','CoM0Yscale_GLS_GLS','CoM1Rtd_MS_MS','fracY_GLS','CoM1GstdTail_MS_MS','MoI0Gtd_GLS_GLS','sY_GLS','MoI1Rtd_GLS_GLS','maxR_GLS','MoI0Ytd_GLS_GLS','CoM1Gscale_MS_MS','aY_GLS','MoI0Rtd_MS_MS','MoI0Gscale_GLS_GLS','aG_GLS','MoI0RavgRes_MS_MS','scaleHead_MS','sG_GLS','maxHead_MS','MoI0YavgRes_GLS_GLS','mR_MS','CoM1Rscale_GLS_GLS','tScale_GLS','MoI1Rtd_MS_MS','CoM0RavgRes_MS_MS','CoM1RavgRes_MS_MS','Gtail_GLS','CoM1RavgRes_GLS_GLS','CoM0Rscale_MS_MS','mYmG_GLS','MoI0Rscale_GLS_GLS','maxG_MS']
# gls_ordered_names = []
# gls_ordered_vals = []
# ms_ordered_names = []
# ms_ordered_vals = []
#
# for i in range(len(ordered_names)):
#     name = ordered_names[i]
#     if name[-2:] == 'LS':
#         gls_ordered_names.append(name[:-4])
#         gls_ordered_vals.append(ordered_vals[i])
#     else:
#         ms_ordered_names.append(name[:-3])
#         ms_ordered_vals.append(ordered_vals[i])
#
# PARAM_REMOVE_GLS = list(np.array(gls_ordered_names)[np.array(gls_ordered_vals) > th])
# PARAM_REMOVE_MS = list(np.array(ms_ordered_names)[np.array(ms_ordered_vals) > th])
# PARAM_NAMES_CONSIDER_GLS=['maxR', 'maxG', 'aG', 'sG', 'Gtail', 'aR', 'sR', 'Rtail', 'aY', 'sY', 'Ytail', 'mRmG', 'mYmG', 'tMove']+\
#                     ['fracR','fracG','fracY', 'tScale']+\
#                     ['CoM0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['CoM1R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['CoM0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['CoM1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI0G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1G{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI0R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1R{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI0Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                     ['MoI1Y{0}_GLS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]#+\
#
# PARAM_NAMES_CONSIDER_MS=['maxG', 'aG', 'mG','sG',\
#                     'maxR', 'aR', 'mR','sR',\
#                     'tMove','mRmG', 'tScale']+\
#                    ['maxHead', 'tailHead', 'scaleHead', 'devHead']+\
#                    ['scaleLength', 'devLength', 'tailLength']+\
#                    ['CoM0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['CoM1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI0G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI1G{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI0R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]+\
#                    ['MoI1R{0}_MS'.format(par) for par in ['td', 'scale', 'avgRes', 'stdTail']]#+\
#
# PARAM_NAMES_USE_GLS = []
# PARAM_NAMES_USE_MS = []
#
# for i in PARAM_NAMES_CONSIDER_GLS:
#     if i not in PARAM_REMOVE_GLS:
#         PARAM_NAMES_USE_GLS.append(i)
#
# for i in PARAM_NAMES_CONSIDER_MS:
#     if i not in PARAM_REMOVE_MS:
#         PARAM_NAMES_USE_MS.append(i)
