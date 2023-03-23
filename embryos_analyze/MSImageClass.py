"""
Created on May 25, 2016

@author: renat

Collection of functions for handling MS embryo images
"""

import os
import numpy as np
from scipy import ndimage
from myFunc import loadImTif, loadImFolder, loadTestCylinderIm, loadTestCylinderImEndOn, loadTestSphereIm, \
    saveImagesMulti, maxIntensProject
from myFigure import myFigure
from scipy.ndimage import gaussian_gradient_magnitude, zoom
from skimage.morphology import skeletonize_3d
from scipy.ndimage.morphology import binary_fill_holes
from scipy.ndimage.measurements import find_objects, label


class MSImage(object):
    """
    class for standard processing and visualization of microscopy images.
    """

    def __init__(self, fileName, z, nT, ch):
        """
        fileName multislice tif file or folder of the file series
        z is a tuple (nZ,dz), where nZ is the number of z planes and dz is the distance between planes in pixels
        t is a tuple (nT,dt), where nT is the number of time points and dt is the interval between time points.
        ch number of channels
        image order has to be z,t,ch
        """

        self.nZ, self.dz = z
        self.nT = nT
        self.nCh = ch
        self.fileName = fileName
        self.images = np.zeros((self.nT, self.nZ, self.nCh, 10, 10)).astype(np.float)
        self.loadIms(fileName)

    def loadIms(self, file_name):
        """
        loads images from file_name into self.image
        :param file_name: Path to a file or folder. Could be a name of a dummy 3d shape generated automatically.
        Note: testS - sphere, testC - cylinder along Z, testCendOn - cylinder along X.
        :return: None
        """
        if os.path.isfile(file_name):
            ims = loadImTif(file_name)
        elif os.path.isdir(file_name):
            ims = loadImFolder(file_name)
        elif file_name == 'testS':
            im = loadTestSphereIm([200, 200], self.nZ, self.nZ * self.dz / 2, self.dz, point=None)
            ims = []
            for i in range(self.nT * self.nCh):
                ims += im
        elif file_name == 'testC':
            im = loadTestCylinderIm([200, 200], self.nZ, self.nZ * self.dz / 2)
            ims = []
            for i in range(self.nT * self.nCh):
                ims += im
        elif file_name == 'testCendOn':
            im = loadTestCylinderImEndOn([200, 200], self.nZ, self.nZ * self.dz / 2)
            ims = []
            for i in range(self.nT * self.nCh):
                ims += im
        else:
            print('Error: wrong filename {0}'.format(file_name))
            return
        if np.array(ims).size == self.nZ * self.nT * self.nCh * ims[0].size:
            self.images = np.reshape(ims, (self.nT, self.nZ, self.nCh, ims[0].shape[0], ims[0].shape[1])).astype(
                np.float)
        else:
            self.images = np.zeros((self.nT, self.nZ, self.nCh, ims[0].shape[0], ims[0].shape[1])).astype(np.float)
            print(
                'Error: number of images (or sizes) does not correspond to z={0}, t={1}, ch={2}'.format(self.nZ,
                                                                                                        self.nT,
                                                                                                        self.nCh))
            raise Exception(
                "Error: number of images (or sizes) does not correspond to z={0}, t={1}, ch={2}".format(self.nZ,
                                                                                                        self.nT,
                                                                                                        self.nCh))

    def saveIms(self, file_name):
        """
        saves images as a multilayered tif file
        :param file_name: file path to save
        :return:
        """
        saveImagesMulti(np.reshape(self.images, (self.nZ * self.nT * self.nCh, self.images[0, 0, 0].shape[0],
                                                 self.images[0, 0, 0].shape[1])), file_name)

    def getEdges(self):
        """
        Uses derivative of gaussian filter to calculate edges in an image
        :return: image stack of edges with the same shape as self.images
        """
        tmp = np.zeros_like(self.images)
        for t in range(self.nT):
            for ch in range(self.nCh):
                for z in range(self.nZ):
                    tmp[t, z, ch] = gaussian_gradient_magnitude(self.images[t, z, ch], sigma=1)
        return tmp

    def getSkeleton(self, ch, th, slide=None):
        """
        finds a 3d skeleton in an image above a threashold
        :param ch: channel
        :param th: threashold for binarization
        :param slide: time slide. If None, all slides are returned
        :return: 3-d image stack with skeletonized image
        """
        tmp = np.zeros((self.nT, self.nZ * self.dz, self.images[0, 0, 0].shape[0], self.images[0, 0, 0].shape[1]))
        if slide is None:
            tRange = range(self.nT)
        else:
            tRange = [slide]
        for t in tRange:
            ims = self.images[t, :, ch].astype(float)
            #             ims = np.array([gaussian_filter(im, 3) for im in ims])
            ims = zoom(ims, (self.dz, 1, 1), order=1)
            #             ims = gaussian_filter(ims, 1)
            ims = (ims > th).astype(np.int)
            ims = binary_fill_holes(ims)
            ims = skeletonize_3d(ims).astype(np.bool)
            #             ims = zoom(ims, (self.dz, 1, 1), order=0)
            tmp[t, :] = ims
        if slide is None:
            return tmp
        else:
            return tmp[slide]

    def getDims(self, ch, th, slide=None):
        """
        Calulates a bounding box in 3d for a stack of images
        :param ch: channel
        :param th: threashold for binarization
        :param slide: time slide. If None, all slides are returned
        :return: (x,y,z) dimensions that encapsulate the fluorescence
        """
        if slide is None:
            tRange = range(self.nT)
        else:
            tRange = [slide]
        dims = []
        for t in tRange:
            ims = self.images[t, :, ch].astype(float)
            #             ims = np.array([gaussian_filter(im, 3) for im in ims])
            if np.sum(ims > th) <= 20:
                dims.append((0, 0, 0))
            else:
                #                 ims = zoom(ims, (self.dz, 1, 1), order=0)
                #             ims = gaussian_filter(ims, 1)
                #                 fig = myFigure()
                #                 fig.imshow(np.max(ims,axis=0))
                #                 fig.show()
                ims = (ims > th).astype(np.int)
                objL, nObj = label(ims)  # label all objects above th
                #                 unique, counts = np.unique(objL, return_counts=True)
                counts = ndimage.sum(ims, objL, range(1, nObj + 1))  # find object areas
                unique = np.arange(1, nObj + 1)
                ims[np.where(np.in1d(objL, unique[np.where(counts < 20)]).reshape(
                    objL.shape))] = 0  # remove objects with small area
                objL[np.where(np.in1d(objL, unique[np.where(counts < 20)]).reshape(
                    objL.shape))] = 0  # remove objects with small area
                bigObj = unique[np.argmax(counts)]  # biggest object
                onEdge = np.unique(objL[:, :, 0], return_counts=False)  # object on the edge
                onEdge = np.concatenate((onEdge, np.unique(objL[:, :, -1], return_counts=False)))
                onEdge = onEdge[onEdge > 0]
                onEdge = onEdge[onEdge != bigObj]  # remove biggest object from consideration
                if onEdge.size > 0:
                    ims[np.where(np.in1d(objL, onEdge).reshape(objL.shape))] = 0  # remove all edge objects from image

                obj = find_objects(ims)
                if len(obj) == 0:
                    dims.append((0, 0, 0))
                else:
                    dims.append(np.array(ims[obj[0]].shape) * np.array([self.dz, 1, 1]))
        return np.array(dims)

    def getMaxProject(self, t, c, th=0):
        """
        Calculates maximum intensity projection
        :param t: time slide
        :param c: channel
        :param th: threashold to truncate intensity
        :return: image
        """
        ims = np.array(self.images[t, :, c])
        ims[np.where(ims < th)] = 0
        return maxIntensProject(ims)

    def calcIntensityMoment(self, t, ch, axis=0, point=None, th=0):
        """
        calculation of an intensity moment - an equivalent of the moment of inertia for pixel intensity from an axis
        passing through a given point in space. Note that axes 2 and 3 are identical, because the moment is calculated
        for a projection of the x axis to bi rotationally independent.

        point: (x, y, z) coordinates. x and y are pixel positions in the image, where x is horizontal axis and y
        vertical. z is the z slice position (may be between slices).
        t: time point
        ch: channel
        """
        if point is None: point = self.getIntCenter(t, ch, th)  # if no center point is assigned, uses center of mass
        x0, y0, z0 = point  # assigns center point to x,y,z coordinate
        z0 *= self.dz  # converts from z planes to pixels
        ims = np.copy(self.images[t, :, ch]).astype(
            np.float32) - th  # makes a copy of ims(to avoid modify in place issues) and subtracts thresh
        ims[np.where(ims < 0)] = 0  # sets region of image that might be negative (after threshold subtraction) to 0
        m = 0.
        for i in range(ims.shape[0]):
            z = i * self.dz  # converts plane to pixels
            x = np.arange(ims[i].shape[1])  # specifies x values for z plane (i)
            y = np.arange(ims[i].shape[0])  # specifies y values for z plane (i)
            xx, yy = np.meshgrid(x, y)  # generates a grid of x,y coordinates with the same shape as image
            if axis == 0:
                rSq = (yy - y0) ** 2 + (
                                           z - z0) ** 2  # defines moments for rotation around x axis, rSq is an array of every vector
            elif axis == 1:
                rSq = (
                          xx - x0) ** 2  # +(z-z0)**2 #gives moments for rotation around y (projection (rotation independent)
            elif axis == 2:
                rSq = (xx - x0) ** 2  # +(yy-y0)**2 #gives moments for rotation around z
            # mArr = np.multiply(rSq, ims[i])
            mArr = np.multiply(rSq, ims[
                i])  # moment of inertia is the sum (all r vectors (rSq) * mass (intensity value))-- sum calculated in next step
            m += np.sum(mArr)  # adds mArr to the value of m (for MoI need to sum Arr
        return m

    def getIntCenter(self, t, ch, th=0):
        '''
        calculation of center of intensity - an equivalent of the center of mass for pixel intensity from a given point in space.
        point: (x, y, z) coordinates. x and y are pixel positions in the image, where x is horizontal axis and y vertical. z is the z slice position (may be between slices).
        t: time point
        ch: channel
        '''
        #         ims = np.copy(self.images[t,:,ch]).astype(np.float32)
        #         ims[np.where(ims<th)] = 0
        ims = np.copy(self.images[t, :, ch]).astype(
            np.float32) - th  # generates a copy of the 3D data array?, converts to 32 but float, subtracts a specified threshold
        ims[np.where(
            ims < 0)] = 0  # finds the index of images where subtraction of threshold results in negative values and sets intensity value equal to 0
        xc, yc, zc = 0, 0, 0  # generates variable for x center, y center, z center
        Itot = np.sum(
            ims)  # sums 3D data (intensity mass) to give the total intensity for all z,t... returns a single value
        if Itot > 0:  # if there is any intensity signal
            x = np.arange(ims[0].shape[
                              1])  # generates a 1D array defined by ims shape for x dimension (i.e 1,2,3,4,5,...200 for 200x200 image)
            y = np.arange(ims[0].shape[0])  # generates a 1D array defined by ims shape for y dimension
            for z in range(ims.shape[0]):  # for each z-step
                xc += np.sum(
                    np.sum(ims[z], axis=0) * x)  # adds up all values in x dimension (for each plane, with for loop)
                yc += np.sum(np.sum(ims[z], axis=1) * y)
                zc += np.sum(ims[z]) * z  # not *dz because z is in plane number dimensions
            xc /= Itot
            yc /= Itot
            zc /= Itot
        return (xc, yc, zc)

    def getTotInt(self, t, ch, th=0):
        im = np.copy(self.images[t, :, ch]).astype(np.float32) - th
        im[np.where(im < 0)] = 0
        return np.sum(im)

    def getAllMoments(self, ch, axis=0, point=None, th=0):
        return np.array([self.calcIntensityMoment(t, ch, axis, point, th) for t in range(self.nT)])

    def getAllIntens(self, ch, th=0):
        return np.array([np.sum(self.images[t, :, ch][np.where(self.images[t, :, ch] >= th)] - th) for t in range(
            self.nT)])  # generates a 1D array of intensity values for each time point (all Z's, for one channel). Only includes timepoints above threshold (np.where gives you the index). Subtracts threshold from each t

    def getIntensProfiles(self, ch):
        IProf = np.zeros([self.nT, self.images[0, 0, 0].shape[1]])
        for slide in range(self.nT):
            IProf[slide] = 1. * np.sum(self.getMaxProject(slide, ch), axis=0) / self.images[0, 0, 0].shape[0]
        return IProf

    def showMaxProj(self, t, ch, th=0):
        fig = myFigure()
        fig.imshow(self.getMaxProject(t, ch, th))
        return fig

    def showImage(self, t, z, ch):
        fig = myFigure()
        fig.imshow(self.images[t, z, ch])
        return fig

    def showMaxProjHist(self, t, ch, th=0, fig=None):
        im = self.getMaxProject(t, ch, th)
        if not fig: fig = myFigure()
        fig.hist(im[np.where(im > 0)].ravel(), bins=30, alpha=0.3)
        return fig


if __name__ == '__main__':
    fileName = 'testC'
    nZ, dz = 18, 8  # 18 planes in 8 pixel steps
    nT, dt = 3, 20  # 3 timepoints, 20 min intervals
    image = MSImage(fileName, (nZ, dz), nT, 1)  # 1 means 1 channel
    R = nZ * dz / 2.  # radius is half of the z-series
    yS, xS = image.images[
        0, 0, 0].shape  # finds the dimensions (shape) of the images. Each image [0 time, 0 z, 0ch] is a 2D array described by y (rows) and x (columns)
    z0, y0, x0 = 18 / 2., yS / 2., xS / 2.  # finds the center point coordinate of the image
    point = (x0, y0, z0)  # redefines the center point with x,y,z convention
    ch = 0  # defines channel to test- channel 0 is green, 1 is red, 2 is brightfield
    th = 0 * 3000
    print('center=', image.getIntCenter(0, 0, 0))
    m0 = image.getAllMoments(ch, axis=0, point=point, th=th)
    m1 = image.getAllMoments(ch, axis=1, point=point, th=th)
    m2 = image.getAllMoments(ch, axis=2, point=point, th=th)
    ints = image.getAllIntens(ch, th)
    fig = myFigure()
    fig.plot(np.arange(nT), m0 / ints / R ** 2, color='g', label='x(0) axis')
    fig.plot(np.arange(nT), m1 / ints / R ** 2, color='r', label='y(1) axis')
    fig.plot(np.arange(nT), m2 / ints / R ** 2, color='b', label='z(2) axis')
    fig.legend(2)
    fig.show()
