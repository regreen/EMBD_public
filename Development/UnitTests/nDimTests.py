'''
Created on Apr 27, 2017

@author: Becky
'''
import unittest
import embdFunc
from embdFunc import getNDAvgVector
import numpy as np
from varLookup import PARAM_NAMES
import random

class nDimTest(unittest.TestCase):
    def testNDdistToCont(self):#tests distance to control
        v2='control'
        d1=[]
        d2=[]
        np.random.seed(0) 
        for v1 in np.random.rand(10,2*len(PARAM_NAMES)): #np.random.rand generates a 10x10 random number array
            ind1= np.array(random.sample(range(len(PARAM_NAMES)),np.random.randint(len(PARAM_NAMES)-1)))#Random
            ind2= np.array(random.sample(range(len(PARAM_NAMES)),np.random.randint(len(PARAM_NAMES)-1)))#Random
            pNamesUse = (np.delete(PARAM_NAMES,ind1),np.delete(PARAM_NAMES,ind2))
            v1[ind1]= np.nan*np.ones(ind1.size)
            v1[(len(PARAM_NAMES)+ind2)]= np.nan*np.ones(ind2.size)
            d1.append(np.linalg.norm(v1[~np.isnan(v1)])) #adds the length of the v1 vector to the d1 list
            d2.append(embdFunc.getNDDistance(v1,v2, pNamesUse))#adds the distance between v1 and control, which should be equivela
        self.assertTrue(d1==d2)
        
    def testNDdistBetweenVectors(self): #tests distance between two genes
        d1 = []
        d2 = []
        np.random.seed(0) 
        for v1,v2 in np.random.rand(10,2,10): #np.random.rand generates 10 tuples of 10 random number arrays (vectors)
            d1.append(np.linalg.norm(v1-v2)) #adds the length of the v1-v2 vector to the d1 list
            d2.append(embdFunc.getNDDistance(v1,v2,(range(10),[])))#adds the distance between v1 and v2, which should be equal to the linalg distance
        self.assertTrue(d1==d2)

    def testNDdimensionalConsistency(self): #when dimensional inconsistency present, checks to see if you flip the dimensional inconsistency if the distance value is the same
        d1 = []
        d2 = []
        np.random.seed(0) 
        for v1,v2 in np.random.rand(10,2,10): #np.random.rand generates 10 tuples of 10 random number arrays (vectors)
            v1+= np.random.choice([0, np.nan],10)#generates an array with values that fall within range defined (can use 2 to make binary or give a choice), that has a length of 10
            d1.append(embdFunc.getNDDistance(v2,v1,(range(10),[]))) #adds the length of the v1-v2 vector to the d1 list
            d2.append(embdFunc.getNDDistance(v1,v2,(range(10),[])))#adds the distance between v1 and v2, which should be equal to the linalg distance
        self.assertTrue(d1==d2)

    def testNDpenalty(self): 
        d1 = []
        d2 = []
        np.random.seed(0) 
        for v1,v2 in np.random.rand(10,2,10):
            vNan= np.random.choice([0, np.nan],10)#generates an array with values that fall within range defined (can use 2 to make binary or give a choice), that has a length of 10
            v1+=vNan
            d1.append(embdFunc.getNDDistance(v1,v2,(range(10),[])))
            v2+=vNan
            d2.append(embdFunc.getNDDistance(v1,v2,(range(10),[])))    
        self.assertTrue((d1>d2).all())
        
    def testNDAvg(self): #tests if average of vectors is properly calculated- this tests equivalent vectors
        d1=[]
        np.random.seed(0) 
        for v1 in np.random.rand(10,10): #makes 10 random vector of length 10
            d1.append(v1) #adds the v1 vector to the d1 list
        vectors= np.array(d1)#makes list of vectors into an array
        v1Avg = np.nanmedian(vectors,axis=0)#calculates the mean of vectors    
        v2Avg = getNDAvgVector(vectors)#uses our function to find the mean of vectors
        self.assertTrue((v1Avg==v2Avg).all())#asserts true if the vectors are the same

    def testNDAvgWnans(self): #tests if average of vectors is properly calculated when one vector contains nans
        d1=[]
        np.random.seed(0) 
        for v1,v2 in np.random.rand(10,2,10): #makes 10 random vector of length 10
            v1+= np.random.choice([0,np.nan],10)#randomly adds nans into the v1 vectors
            d1.append(v1) #adds the v1 vector (which contains nans) to the d1 list
            d1.append(v2) #adds the v2 vector to the d2 list (doesnt contain nans)
        vectors= np.array(d1)#makes list of vectors into an array
        v1Avg = np.nanmedian(vectors,axis=0)#calculates the mean of vectors    
        v2Avg = getNDAvgVector(vectors)#uses our function to find the mean of vectors
        self.assertTrue((v1Avg==v2Avg).all())#asserts true if the vectors are the same

    def testNDAvgAllNans(self): # tests if n vectors are composed entirely of nans, should return nans
        d1= []
        for i in range(10):
            d1.append(np.empty(10) * np.nan)#makes a list of 10 nans
        vectors=np.array(d1)#makes an array of all the nan vectors
        v1Avg = np.isnan(np.empty(10) * np.nan) #np.isnan returns true if nan... np.empty(10) makes a list of 10 all nans    
        v2Avg = np.isnan(getNDAvgVector(vectors))#uses our function to find the mean of vectors, checks to see if nan... if so returns true (this is needed because nan!=nan, but true==true)
        self.assertTrue((v1Avg==v2Avg).all())#asserts true if the vectors are the same
        
    #test PAD


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()