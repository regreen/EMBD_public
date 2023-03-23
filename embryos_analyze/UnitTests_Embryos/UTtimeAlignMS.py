'''
Created on Apr 4, 2017

@author: Admin
'''
import unittest
from varLookup import *
from Embryos import MSEmbryo

tol = 0.2
tol_t0 = 1.5
dt = -2
refresh = False

class TestEmbAlignment(unittest.TestCase):

    def testEmb01(self):
        tScaleMan = 1.05
        t0_man = 6.56
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140515T131046/Emb4/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, tScale_manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
        
    # def testEmb02(self):
    #     tScaleMan = 0.89
    #     embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140509T143953/Emb3/'
    #     emb = MSEmbryo(embFolder)
    #     if refresh: emb.refresh()
    #     else: emb.updateParams()
    #     print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
    #     self.assertTrue(tScaleMan-tol<emb.tScale<tScaleMan+tol)
    #
    # def testEmb03(self):
    #     tScaleMan = 0.83
    #     embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140514T130911/Emb4/'
    #     emb = MSEmbryo(embFolder)
    #     if refresh: emb.refresh()
    #     else: emb.updateParams()
    #     print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
    #     self.assertTrue(tScaleMan-tol<emb.tScale<tScaleMan+tol)
    #
    # def testEmb04(self):
    #     tScaleMan = 1.
    #     embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140425T140657/Emb4/'
    #     emb = MSEmbryo(embFolder)
    #     if refresh: emb.refresh()
    #     else: emb.updateParams()
    #     print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
    #     self.assertTrue(tScaleMan-tol<emb.tScale<tScaleMan+tol)
             
    def testEmb05(self):
        tScaleMan = 0.86
        t0_man = 7.74
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140416T140401/Emb3/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
        
    def testEmb06(self):
        tScaleMan = 0.99
        t0_man = 7.3
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140430T140422/Emb1/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
        
    def testEmb07(self):
        tScaleMan = 0.81
        t0_man = 10.12
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140416T140401/Emb1/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
         
    def testEmb08(self):
        tScaleMan = 1.05
        t0_man = 3.69
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140515T131046/Emb3/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
        
#     def testEmb09(self):
#         tScaleMan = 1.15
    #     t0_man = 8.68
    #     embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140509T143953/Emb4/'
#         emb = MSEmbryo(embFolder)
#         if refresh: emb.refresh()
#         else: emb.updateParams()
#         print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
#         self.assertTrue(tScaleMan-tol<emb.tScale<tScaleMan+tol)
        
    def testEmb10(self):
        tScaleMan = 1.20
        t0_man = 9.09
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140514T130911/Emb5/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
        
    def testEmb11(self):
        tScaleMan = 0.97
        t0_man = 6.45
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140402T140154/Emb2/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
        
    def testEmb12(self):
        tScaleMan = 0.90
        t0_man = 9.21
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140402T140154/Emb1/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
        
    def testEmb13(self):
        tScaleMan = 0.98
        t0_man = 4.11
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140424T140159/Emb1/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
        
    def testEmb14(self):
        tScaleMan = 0.82
        t0_man = 8.36
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140501T135409/Emb1/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
        
    def testEmb15(self):
        tScaleMan = 1.23
        t0_man = 7.38
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140418T143933/Emb1/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
        
    def testEmb16(self):
        tScaleMan = 0.86
        t0_man = 9.36
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140508T141139/Emb1/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
        
    def testEmb17(self):
        tScaleMan = 1.09
        t0_man = 10.41
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140507T121259/Emb1/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
        
    def testEmb18(self):
        tScaleMan = 0.85
        t0_man = 9.28
        embFolder = FOLDER_IN + 'cropped/EMBD0000/MS/20140430T140422/Emb3/'
        emb = MSEmbryo(embFolder)
        if refresh: emb.refresh()
        else: emb.updateParams()
        print('tScale={0}, manual={1}'.format(emb.tScale, tScaleMan))
        print('t0={0}, t0_manual={1}'.format(emb.t0+dt, t0_man))
        self.assertTrue((tScaleMan-tol < emb.tScale < tScaleMan+tol)*(t0_man-tol_t0 < emb.t0 + dt < t0_man+tol_t0))
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()