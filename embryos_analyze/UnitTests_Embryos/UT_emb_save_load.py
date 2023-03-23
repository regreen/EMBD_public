import unittest
from Embryos import MSEmbryo, GSEmbryo
from varLookup import FOLDER_IN
import numpy as np


class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.ms_emb_folder = FOLDER_IN + 'cropped/EMBD0000/MS/20140515T131046/Emb4/'
        self.ms_emb1 = MSEmbryo(self.ms_emb_folder)
        self.ms_emb1.refresh()
        self.ms_emb2 = MSEmbryo(self.ms_emb_folder)

        self.gls_emb_folder = FOLDER_IN + 'cropped/EMBD0000/GLS/20140515T131046/Emb1/'
        self.gls_emb1 = GSEmbryo(self.gls_emb_folder)
        self.gls_emb1.refresh()
        self.gls_emb2 = GSEmbryo(self.gls_emb_folder)

    """tests for ms embs"""
    def test_tInt(self):
        self.assertTrue((self.ms_emb1.tInt == self.ms_emb2.tInt).all(), 'MS tInt are not equal')

    def test_stopInt(self):
        self.assertTrue((self.ms_emb1.stopInt == self.ms_emb2.stopInt), 'MS stopInt is not equal')

    def test_image(self):
        self.assertTrue((self.ms_emb1.image == self.ms_emb2.image), 'MS image is not equal')

    def test_loaded(self):
        self.assertFalse((self.ms_emb1.loaded == self.ms_emb2.loaded), 'MS loaded ARE equal')

    def test_radius(self):
        self.assertTrue((self.ms_emb1.radius == self.ms_emb2.radius), 'MS radius is not equal')

    def test_length(self):
        self.assertTrue((self.ms_emb1.length == self.ms_emb2.length), 'MS length is not equal')

    def test_version(self):
        self.assertTrue((self.ms_emb1.__version__ == self.ms_emb2.__version__), 'MS version is not equal')

    def test_thresh(self):
        self.assertEqual(self.ms_emb1.thresh, self.ms_emb2.thresh, 'MS thresh is not equal')

    def test_t0(self):
        self.assertTrue((self.ms_emb1.t0 == self.ms_emb2.t0), 'MS t0 is not equal')

    def test_tScale(self):
        self.assertTrue((self.ms_emb1.tScale == self.ms_emb2.tScale), 'MS tScale is not equal')

    def test_tMove(self):
        self.assertTrue((self.ms_emb1.tMove == self.ms_emb2.tMove), 'MS tMove is not equal')

    def test_movement(self):
        self.assertTrue((self.ms_emb1.movement == self.ms_emb2.movement), 'MS movement is not equal')

    def test_dead(self):
        self.assertTrue((self.ms_emb1.dead == self.ms_emb2.dead), 'MS dead is not equal')

    def test_startPoint(self):
        self.assertDictEqual(self.ms_emb1.startPoint, self.ms_emb2.startPoint, 'MS startPoint is not equal')

    def test_scale(self):
        self.assertDictEqual(self.ms_emb1.scale, self.ms_emb2.scale, 'MS scale is not equal')

    def test_cutoff(self):
        self.assertDictEqual(self.ms_emb1.cutoff, self.ms_emb2.cutoff, 'MS cutoff is not equal')

    def test_tIntTh(self):
        self.assertTrue((self.ms_emb1.tIntTh == self.ms_emb2.tIntTh).all(), 'MS tIntTh is not equal')

    def test_dims(self):
        self.assertTrue((self.ms_emb1.dims == self.ms_emb2.dims).all(), 'MS dims is not equal')

    def test_headInt(self):
        self.assertTrue((self.ms_emb1.headInt == self.ms_emb2.headInt).all(), 'MS headInt is not equal')

    def test_CMPos(self):
        self.assertTrue((self.ms_emb1.CMPos == self.ms_emb2.CMPos).all(), 'MS CMPos is not equal')

    def test_MoI(self):
        self.assertTrue((self.ms_emb1.MoI == self.ms_emb2.MoI).all(), 'MS MoI is not equal')

    def test_params(self):
        flag = True
        for key in self.ms_emb2.params:
            if self.ms_emb2.params[key] != self.ms_emb1.params[key] and not np.isnan(self.ms_emb2.params[key]) * np.isnan(self.ms_emb1.params[key]):
                print('key =', key, self.ms_emb2.params[key], self.ms_emb1.params[key])
                flag = False
            self.assertTrue((flag == True), 'MS params are not equal')

    def test_embCenter(self):
        self.assertTrue((self.ms_emb1.embCenter == self.ms_emb2.embCenter).all(), 'MS embCenter is not equal')

    def test_nLambda(self):
        self.assertEqual(self.ms_emb1.nLambda, self.ms_emb2.nLambda, 'MS nLambda is not equal')

    def test_tauScale(self):
        self.assertEqual(self.ms_emb1.tauScale, self.ms_emb2.tauScale, 'MS tauScale is not equal')


    """tests for gls embs"""
    def test_gls_tInt(self):
        self.assertTrue((self.gls_emb1.tInt == self.gls_emb2.tInt).all(), 'GLS tInt are not equal')

    def test_gls_stopInt(self):
        self.assertTrue((self.gls_emb1.stopInt == self.gls_emb2.stopInt), 'GLS stopInt is not equal')

    def test_gls_image(self):
        self.assertTrue((self.gls_emb1.image == self.gls_emb2.image), 'GLS image is not equal')

    def test_gls_loaded(self):
        self.assertFalse((self.gls_emb1.loaded == self.gls_emb2.loaded), 'GLS loaded ARE equal')

    def test_gls_radius(self):
        self.assertTrue((self.gls_emb1.radius == self.gls_emb2.radius), 'GLS radius is not equal')

    def test_gls_length(self):
        self.assertTrue((self.gls_emb1.length == self.gls_emb2.length), 'GLS length is not equal')

    def test_gls_version(self):
        self.assertTrue((self.gls_emb1.__version__ == self.gls_emb2.__version__), 'GLS version is not equal')

    def test_gls_thresh(self):
        self.assertEqual(self.gls_emb1.thresh, self.gls_emb2.thresh, 'GLS thresh is not equal')

    def test_gls_t0(self):
        self.assertTrue((self.gls_emb1.t0 == self.gls_emb2.t0), 'GLS t0 is not equal')

    def test_gls_tScale(self):
        self.assertTrue((self.gls_emb1.tScale == self.gls_emb2.tScale), 'GLS tScale is not equal')

    def test_gls_tMove(self):
        self.assertTrue((self.gls_emb1.tMove == self.gls_emb2.tMove), 'GLS tMove is not equal')

    def test_gls_movement(self):
        self.assertTrue((self.gls_emb1.movement == self.gls_emb2.movement), 'GLS movement is not equal')

    def test_gls_dead(self):
        self.assertTrue((self.gls_emb1.dead == self.gls_emb2.dead), 'GLS dead is not equal')

    def test_gls_startPoint(self):
        self.assertDictEqual(self.gls_emb1.startPoint, self.gls_emb2.startPoint, 'GLS startPoint is not equal')

    def test_gls_scale(self):
        self.assertDictEqual(self.gls_emb1.scale, self.gls_emb2.scale, 'GLS scale is not equal')

    def test_gls_cutoff(self):
        self.assertDictEqual(self.gls_emb1.cutoff, self.gls_emb2.cutoff, 'GLS cutoff is not equal')

    def test_gls_tIntTh(self):
        self.assertTrue((self.gls_emb1.tIntTh == self.gls_emb2.tIntTh).all(), 'GLS tIntTh is not equal')

    def test_gls_dims(self):
        self.assertTrue((self.gls_emb1.dims == self.gls_emb2.dims).all(), 'GLS dims is not equal')

    def test_gls_headInt(self):
        self.assertTrue((self.gls_emb1.headInt == self.gls_emb2.headInt).all(), 'GLS headInt is not equal')

    def test_gls_CMPos(self):
        self.assertTrue((self.gls_emb1.CMPos == self.gls_emb2.CMPos).all(), 'GLS CMPos is not equal')

    def test_gls_MoI(self):
        self.assertTrue((self.gls_emb1.MoI[np.invert(np.isnan(self.gls_emb1.MoI))] ==
                         self.gls_emb2.MoI[np.invert(np.isnan(self.gls_emb2.MoI))]).all(), 'GLS MoI is not equal')

    def test_gls_params(self):
        flag = True
        for key in self.gls_emb2.params:
            if self.gls_emb2.params[key] != self.gls_emb1.params[key] and not np.isnan(self.gls_emb2.params[key]) * \
                    np.isnan(self.gls_emb1.params[key]):
                print('key =', key, self.gls_emb2.params[key], self.gls_emb1.params[key])
                flag = False
            self.assertTrue((flag == True), 'GLS params are not equal')

    def test_gls_nLambda(self):
        self.assertEqual(self.gls_emb1.nLambda, self.gls_emb2.nLambda, 'GLS nLambda is not equal')

    def test_gls_tauScale(self):
        self.assertEqual(self.gls_emb1.tauScale, self.gls_emb2.tauScale, 'GLS tauScale is not equal')

    def test_gls_green(self):
        self.assertTrue((self.gls_emb1.green == self.gls_emb2.green).all(), 'GLS green is not equal')

    def test_gls_red(self):
        self.assertTrue((self.gls_emb1.red == self.gls_emb2.red).all(), 'GLS red is not equal')

    def test_gls_yellow(self):
        self.assertTrue((self.gls_emb1.yellow == self.gls_emb2.yellow).all(), 'GLS yellow is not equal')

    def test_gls_greenPos(self):
        self.assertTrue((self.gls_emb1.greenPos == self.gls_emb2.greenPos).all(), 'GLS greenPos is not equal')

    def test_gls_redPos(self):
        self.assertTrue((self.gls_emb1.redPos == self.gls_emb2.redPos).all(), 'GLS redPos is not equal')

    def test_gls_yellowPos(self):
        self.assertTrue((self.gls_emb1.yellowPos == self.gls_emb2.yellowPos).all(), 'GLS yellowPos is not equal')


if __name__ == '__main__':
    unittest.main()
