#!/usr/bin/env python3

import unittest

from utils import isdnaorproteins, get_pssm_fam


class TestIsDnaOrProteins(unittest.TestCase):
    def test_dna_ok(self):
        """
        Test that the function runs smoothly
        """
        data = "ATCGATCG"
        self.assertEqual(isdnaorproteins(data), "DNA")

    def test_proteins_ok(self):
        """
        Test that the function runs smoothly
        """
        data = "MKRCCSL"
        self.assertEqual(isdnaorproteins(data), "protein")

    def test_not_ok(self):
        """
        Test that the function returns correctly an error
        when needed
        """
        data = "XYZGRTE"
        self.assertEqual(isdnaorproteins(data), "unknown")


class TestGetPSSMFam(unittest.TestCase):
    def test_get_correct_fam(self):
        """
        Test we correctly get the families
        """
        mydict = {"ID1": ["A", "A", "B", "M"], "ID2": ["M", "P", "O1", "O1"]}
        resdict = {"ID1": "A", "ID2": "O1"}

        self.assertEqual(get_pssm_fam(mydict), resdict)

    def test_get_incorrect_fam(self):
        """
        Test we correctly get the families
        """
        mydict = {"ID1": ["A", "A", "B", "M"], "ID2": ["M", "P", "O1", "O1"]}
        resdict = {"ID1": "M", "ID2": "P"}

        self.assertNotEqual(get_pssm_fam(mydict), resdict)


if __name__ == "__main__":
    unittest.main()
