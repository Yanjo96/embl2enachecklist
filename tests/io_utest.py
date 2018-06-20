import unittest
import sys, os

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'embl2enachecklists'))
import Embl2enachecklistsMain as EMBL2ENAclMain

class IOTest(unittest.TestCase):

    def test_embl2enachecklists(self):
        self.assertTrue(EMBL2ENAclMain.embl2enachecklists('examples/input/matK.embl', 'examples/output/matK_SubmissionChecklist.tsv', 'trnK_matK'))


if __name__ == '__main__':
    unittest.main()
