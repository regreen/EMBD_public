'''
Created on Jul 25, 2016

@author: renat
'''
from embdFunc import *
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
from matplotlib import gridspec
from varLookup import FOLDER_IN
from emb_handler import loadEmbs

load_file = '''fill in here'''
savefig_path = '''fill in here'''

if __name__ == '__main__':
    #     emb = loadEmbs(load_file)[0]
    emb = loadEmbs(FOLDER_IN + 'cropped/EMBD0028/MS/Emb3/')[0]
    #     emb.thresh= [0,0]
    emb.loadImages()
    dist1 = np.array(emb.getDistFromGreen(1, 1))
    dist2 = np.array(emb.getDistFromGreen(1, 2))
    dist = (dist1 + dist2) / 2.
    d = 0
    #     d = 4
    for t in range(d, 31):
        fig = plt.figure(figsize=(4, 6), facecolor='black')
        rc('axes', linewidth=2)
        rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial'], 'size': 18})
        paramstring = r'\usepackage{bm}'
        matplotlib.rcParams['text.latex.preamble'] = paramstring
        matplotlib.rcParams['svg.fonttype'] = 'none'
        rc('text', usetex=False)
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
        axMovie = fig.add_subplot(gs[0])
        axCM = fig.add_subplot(gs[1])
        imG = np.array(emb.image.getMaxProject(t, 0), dtype=np.float32) / 20000.
        imR = np.array(emb.image.getMaxProject(t, 1), dtype=np.float32) / 30000.
        imG[np.where(imG > 1.)] = 1.
        imR[np.where(imR > 1.)] = 1.
        imB = np.zeros_like(imG)
        xG, yG, zG = emb.image.getIntCenter(t, 0, emb.thresh[0])
        xR, yR, zR = emb.image.getIntCenter(t, 1, emb.thresh[1])
        axMovie.imshow(np.dstack([imR, imG, imB]))
        axMovie.scatter([xG], [yG], color='greenyellow', marker='*', s=100)
        axMovie.scatter([xR], [yR], color='lightcoral', marker='*', s=100)
        axMovie.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
        axMovie.yaxis.set_major_locator(matplotlib.ticker.NullLocator())
        axMovie.axis('off')
        axCM.plot(np.arange(t - d) * 20, dist[d:t], color='r', linewidth=2)
        axCM.set_xlim((0, 500))
        axCM.set_ylim((0, 1.))
        #         axCM.set_yticks([0,0.2,0.4,0.6,0.8,0.1])
        axCM.tick_params(axis='both',  # changes apply to the x-axis and y-axis
                         which='both',  # both major and minor ticks are affected
                         direction='out',
                         left='on',  # ticks along the bottom edge are off
                         right='off',  # ticks along the top edge are off
                         bottom='on',
                         top='off',
                         labelleft='on',
                         labelright='off',
                         labeltop='off',
                         labelbottom='on')
        axCM.set_facecolor('k')
        axCM.spines['bottom'].set_color('white')
        axCM.spines['left'].set_color('white')
        axCM.tick_params(axis='x', colors='white')
        axCM.tick_params(axis='y', colors='white')
        axCM.yaxis.label.set_color('white')
        axCM.xaxis.label.set_color('white')
        #         axMovie.set_title('$pop-1(RNAi$)')
        axMovie.set_title('Control')
        axMovie.title.set_color('white')
        axCM.set_xlabel('Time (min)')
        axMovie.title.set_fontsize(25)
        fig.tight_layout()
        # plt.show()
        #         fig.savefig(savefig_path + '{0:02}.png'.format(t), transparant = True, bbox_inches='tight',pad_inches = 0,facecolor='black')
        fig.savefig(FOLDER_IN + 'embdev_images/CoM/Movie/0028CM/0028CM_{0:02}.png'.format(t), transparant=True,
                    bbox_inches='tight', pad_inches=0, facecolor='black')
        plt.close(fig)
