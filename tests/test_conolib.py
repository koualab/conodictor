import unittest

from conodictor.conolib import isdnaorproteins


class TestIsDnaOrProteins(unittest.TestCase):
    """Unit test for isdnaorproteins."""

    def test_dna_sequence(self):
        assert isdnaorproteins("ATCG") == "DNA"
        assert isdnaorproteins("GATTACA") == "DNA"
        assert isdnaorproteins("CGTAGCTAGCTG") == "DNA"

    def test_protein_sequence(self):
        assert isdnaorproteins("MEEPQSDPSV") == "protein"
        assert isdnaorproteins("LIPKAD") == "protein"
        assert isdnaorproteins("ABCDEFGHIKLMNPQRSTVWYZ*X") == "protein"

    def test_mixed_sequence(self):
        assert isdnaorproteins("ATCBXO") == "unknown"
        assert isdnaorproteins("MEEGATTACA") == "unknown"

    def test_invalid_characters(self):
        assert isdnaorproteins("12345") == "unknown"
        assert isdnaorproteins("!@#$%") == "unknown"
        assert isdnaorproteins("ATCG!") == "unknown"


if __name__ == "__main__":
    unittest.main()
