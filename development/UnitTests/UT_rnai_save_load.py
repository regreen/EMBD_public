import unittest
from RNAiClass import *


class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        i = 63
        self.rnai_1 = RNAiClass(i)
        self.rnai_1.refresh_params_data()  # refreshes and saves params data to mySQL
        self.rnai_2 = RNAiClass(i)  # loads params from mySQL

    def test_version(self):
        self.assertTrue(self.rnai_1.__version__ == self.rnai_2.__version__, 'versions are not equal')

    def test_refreshed(self):
        self.assertFalse(self.rnai_1.refreshed == self.rnai_2.refreshed, 'refreshed ARE equal')

    def test_paramsGLS(self):
        flag = True
        for key in self.rnai_2.paramsGLS:
            if self.rnai_2.paramsGLS[key] != self.rnai_1.paramsGLS[key] and not np.isnan(
                    self.rnai_2.paramsGLS[key]) * np.isnan(self.rnai_1.paramsGLS[key]):
                print('key =', key, self.rnai_2.paramsGLS[key], self.rnai_1.paramsGLS[key])
                flag = False
            self.assertTrue(flag, 'GLS params are not equal')

    def test_prev_MS(self):
        self.assertTrue(self.rnai_1.prev_MS == self.rnai_2.prev_MS, 'prev_MS are not equal')

    def test_prev_GLS(self):
        self.assertTrue(self.rnai_1.prev_GLS == self.rnai_2.prev_GLS, 'prev_GLS are not equal')

    def test_paramsMS(self):
        flag = True
        for key in self.rnai_2.paramsMS:
            if self.rnai_2.paramsMS[key] != self.rnai_1.paramsMS[key] and not np.isnan(
                    self.rnai_2.paramsMS[key]) * np.isnan(self.rnai_1.paramsMS[key]):
                print('key =', key, self.rnai_2.paramsMS[key], self.rnai_1.paramsMS[key])
                flag = False
            self.assertTrue(flag, 'paramsMS are not equal')





if __name__ == '__main__':
    unittest.main()
