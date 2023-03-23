'''
Created on Jun 15, 2014

@author: renat

Collection of various commonly used functions
'''

import os, glob, re, cv2, errno, platform
from numpy import pi
import numpy as np
from myFigure import myFigure
import shutil
import scipy.ndimage.filters as filters
import scipy.ndimage as ndimage
# from fitEllipse import create_ellipse
from tifffile import imsave, imread
from PIL import Image
# from lmfit import Parameter, minimize, Parameters
# from time import mktime
from datetime import datetime
import time
from functools import wraps


def clearFolder(folder, subFolder=False):
    if os.path.exists(folder):
        for the_file in os.listdir(folder):
            file_path = os.path.join(folder, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif subFolder and os.path.isdir(file_path): shutil.rmtree(file_path)
            except:
                pass
    else: os.mkdir(folder)


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def loadImFolder(folder, fmt='tif'):
    imNames = glob.glob('{0}*.{1}'.format(folder, fmt))
    if len(imNames) > 0:
        sort_nicely(imNames)
        im = cv2.imread(imNames[0], -1)
        if len(im.shape) == 2:
            images = np.array([im])
        else:
            images = im
        for name in imNames[1:]:
            im = cv2.imread(name, -1)
            if im is not None and len(im.shape) == 2:
                images = np.concatenate((images, np.array([im])))
            else:
                images = np.concatenate((images, im))
        return images
    else:
        print('Msg: Can not find images in {0}'.format(folder))
        return None


def loadImTif(fileName):
    return imread(fileName)


def loadTestSphereIm(shape, zN, R, dz, point=None):
    ims = []
    if point is None:
        zc = zN / 2
        xc = shape[0] / 2
        yc = shape[1] / 2
    else:
        xc, yc, zc = point

    for i in range(zN):
        im = np.zeros(shape)  # creates image of zeros with defined shape
        z = i * dz  # converts z plane to pixels
        r = np.sqrt(R ** 2 - (zc - z) ** 2 * dz ** 2)  # defines R
        x = np.arange(shape[1])  # 1D array defined by shape for x dimension (i.e 1,2,3,4,5,...200 for 200x200 image)
        y = np.arange(shape[0])  # 1D array defined by shape for y dimension (i.e 1,2,3,4,5,...200 for 200x200 image)
        xx, yy = np.meshgrid(x,
                             y)  # creates all possible x,y coordinates (i.e x=[1,2,3,4,5]y[1,1,1,1,1],x=[1,2,3,4,5]y[2,2,2,2,2],etc)
        dist = np.sqrt((xx - xc) ** 2 + (yy - yc) ** 2)  # finds distance for all possible coordinates to origin
        im[np.where(dist <= r)] = 1  # sets all pixels within the radius of the circle = 1
        ims.append(im)  # appends image to ims array
    return ims


def loadTestCylinderIm(shape, zN, R):
    ims = []
    dz = 2. * R / zN
    xc = shape[0] / 2
    yc = shape[1] / 2
    for i in range(zN):
        im = np.zeros(shape)
        z = i * dz
        r = R
        x = np.arange(shape[1])
        y = np.arange(shape[0])
        xx, yy = np.meshgrid(x, y)
        dist = np.sqrt((xx - xc) ** 2 + (yy - yc) ** 2)
        im[np.where(dist <= r)] = 1
        ims.append(im)
    return ims


def loadTestCylinderImEndOn(shape, zN, R):
    ims = []  # empty list, to be populated with individual images from z series
    dz = 2. * R / zN  # change in z step, defined by twice the radius divided by the number of z steps
    zc = zN / 2 * dz
    xc = shape[1] / 2  # x center is the midpoint in the x dimension
    yc = shape[0] / 2  # y center is the midpoint in the y dimension
    for i in range(zN):  # for plane in z range
        im = np.zeros(shape)  # generate a blank image of specified shape
        z = i * dz  # converts z from plane number to pixels
        r = R  # radius
        x = np.arange(shape[1])  # x values
        y = np.arange(shape[0])  # y values
        xx, yy = np.meshgrid(x,
                             y)  # creates all possible x,y coordinates (i.e x=[1,2,3,4,5]y[1,1,1,1,1],x=[1,2,3,4,5]y[2,2,2,2,2],etc)
        xdist = np.sqrt(
            (xx - xc) ** 2)  # measures distances to all possible coordinates in x direction (measure of length)
        # length of the cylinder is 2R
        ydist = np.sqrt((yy - yc) ** 2)  # measures distances to all possible coordinates in y direction (measure of
        # width of the cylinder is 2x which is 2.*np.sqrt(r**2-(z-zc)**2)
        im[(xdist <= r) & (ydist < np.sqrt(r ** 2 - (
        z - zc) ** 2))] = 1  # transforms image array into true/false based on specified requirements (length (x dist)) and width (y dist))
        ims.append(im)  # appends image to ims list
    #     fig =myFigure()#generates an empty figure
    #     fig.imshow(ims[zN/2])#plots middle image in my figure
    #     fig.show()    #shows image
    return ims


def tryint(s):
    try:
        return int(s)
    except:
        return s


def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [tryint(c) for c in re.split('([0-9]+)', s)]


def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)


def maxIntensProject(imList):
    """ Maximum intensity projection.
    INPUT:
    imList: list of images in numpy array form
    OUTPUT:
    single image of the size of input images
    """

    return np.max(imList, axis=0)



def saveImages(imgs, filePrefix, folder):
    """
    saves images
    Input:
    imgs: list of images as numpy arrays
    filePrefix: string prefix to use in front of the file name
    folder: folder to save images into
    """
    for i in range(len(imgs)):
        cv2.imwrite(folder + filePrefix + '{0:0>3}.tif'.format(i), imgs[i])


def saveImagesMulti(imgs, fileName):
    imsave(fileName, np.array(imgs))


def getStrVal(x, xerr):
    import myMath
    if xerr is None or np.isnan(xerr) or np.isinf(xerr):
        return '{0}+/-inf'.format(x)
    elif xerr == 0:
        return '{0}+/-0'.format(x)
    else:
        precis = max(0, max(myMath.getSig(x), myMath.getSig(xerr)))
    if x is None or np.isnan(x):
        return '{1}+/-{2}'.format(precis, x, xerr)
    else:
        return '{1:.{0}f}+/-{2:.{0}f}'.format(precis, np.round(x, precis), np.round(xerr, precis))


def getLargestObj(mask):
    label_im, nb_labels = ndimage.label(mask)
    if nb_labels > 1:
        sizes = ndimage.sum(mask, label_im, range(nb_labels + 1))
        mask_size = sizes < np.max(sizes)
        remove_pixel = mask_size[label_im]
        label_im[remove_pixel] = 0
        labels = np.unique(label_im)
        label_im = np.searchsorted(labels, label_im)
    return label_im


def find_consec_inds(inds, n):
    """ finds consecutive indicies in a list
    INPUT:
    inds: a list of indicies
    n: the number of consecutive indicies to check for
    OUTPUT:
    returns True if consecutive indices are found and returns false if they are not. Returns pos: index position(s) where n consecutive inds are present
    """

    all_rolls = []
    a = np.roll(inds, 1)  # first roll by one
    b = np.abs(inds - a)  # subtract to find which indicies differ by 1
    c = (b == 1)  # converts to  true, false list (0,1,1,0,1)
    all_rolls.append(c)  # appends the first roll
    for i in range(n - 2):
        d = np.roll(c, 1)
        all_rolls.append(d)
        c = d
    sum_all = np.sum(all_rolls, axis=0)
    if n-1 not in sum_all:
        return False, []
    else:
        end_pos = np.where(sum_all >= n-1)[0]
        pos = end_pos - (n-1)
        return True, pos


def readCSV(fileName, delimiter = ' '):
    import csv

    csvFile = csv.reader(open(fileName, 'r'), delimiter=delimiter)
    result=[]
    for row in csvFile:
        result.append([])
        for val in row:
            try:
                if val!='None' and val!='nan' and val!='NA': result[-1].append(float(val))
                else: result[-1].append(np.nan)
            except:
                result[-1].append(val)
    return result


if __name__ == '__main__':
    pass
