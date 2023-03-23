'''
Created on Mar 30, 2017

@author: Admin
'''
import unittest
from varLookup import *
from Embryos import GSEmbryo

tol = 0.2

class TestEmbAlignment(unittest.TestCase):

    def testEmb1(self):
        tScaleMan = 1.
        embFolder = FOLDER_IN + 'cropped/EMBD0000/GLS/20140515T131046/Emb2/'
        emb = GSEmbryo(embFolder)
        emb.refresh()
        self.assertTrue(tScaleMan-tol<emb.tScale<tScaleMan+tol)

    def testEmb2(self):
        tScaleMan = 1.19
        embFolder = FOLDER_IN + 'cropped/EMBD0000/GLS/20140528T141502/Emb3/'
        emb = GSEmbryo(embFolder)
        emb.refresh()
        self.assertTrue(tScaleMan-tol<emb.tScale<tScaleMan+tol)

    def testEmb3(self):
        tScaleMan = 0.86
        embFolder = FOLDER_IN + 'cropped/EMBD0000/GLS/20140515T131046/Emb1/'
        emb = GSEmbryo(embFolder)
        emb.refresh()
        self.assertTrue(tScaleMan-tol<emb.tScale<tScaleMan+tol)

    def testEmb4(self):
        tScaleMan = 1.12
        embFolder = FOLDER_IN + 'cropped/EMBD0000/GLS/20140516T130559/Emb2/'
        emb = GSEmbryo(embFolder)
        emb.refresh()
        self.assertTrue(tScaleMan-tol<emb.tScale<tScaleMan+tol)

    def testEmb5(self):
        tScaleMan = 0.83
        embFolder = FOLDER_IN + 'cropped/EMBD0000/GLS/20140529T134135/Emb2/'
        emb = GSEmbryo(embFolder)
        emb.refresh()
        self.assertTrue(tScaleMan-tol<emb.tScale<tScaleMan+tol)
    
    def testEmb6(self):
        tScaleMan = 1.
        embFolder = FOLDER_IN + 'cropped/EMBD0000/GLS/20140508T141139/Emb1/'
        emb = GSEmbryo(embFolder)
        emb.refresh()
        self.assertTrue(tScaleMan-tol<emb.tScale<tScaleMan+tol)
    
    def testEmb7(self):
        tScaleMan = 1.
        embFolder = FOLDER_IN + 'cropped/EMBD0000/GLS/20140501T135409/Emb1/'
        emb = GSEmbryo(embFolder)
        emb.refresh()
        self.assertTrue(tScaleMan-tol<emb.tScale<tScaleMan+tol)
        
    def testEmb8(self):
        tScaleMan = 0.8
        embFolder = FOLDER_IN + 'cropped/EMBD0000/GLS/20140701T154139/Emb2/'
        emb = GSEmbryo(embFolder)
        emb.refresh()
        self.assertTrue(tScaleMan-tol<emb.tScale<tScaleMan+tol)
        
    def testEmb9(self):
        tScaleMan = 1.1
        embFolder = FOLDER_IN + 'cropped/EMBD0000/GLS/20140430T140422/Emb1/'
        emb = GSEmbryo(embFolder)
        emb.refresh()
        self.assertTrue(tScaleMan-tol<emb.tScale<tScaleMan+tol)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()