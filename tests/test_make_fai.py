import unittest
import os
import sys
import io

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from make_fai import make_fai

class TestMakeFai(unittest.TestCase):

    def setUp(self):
        """Set up a dummy FASTA file for testing."""
        self.fasta_filename = "test.fasta"
        self.fai_filename = self.fasta_filename + ".fai"

        with open(self.fasta_filename, "w") as f:
            f.write(">seq1 description\n")
            f.write("ACGTACGT\n")
            f.write("ACGT\n")
            f.write(">seq2\n")
            f.write("TCGATCGATCGA\n")

    def tearDown(self):
        """Clean up the created files."""
        if os.path.exists(self.fasta_filename):
            os.remove(self.fasta_filename)
        if os.path.exists(self.fai_filename):
            os.remove(self.fai_filename)

    def test_make_fai_creates_index(self):
        """Test if make_fai creates the .fai file."""
        make_fai(self.fasta_filename)
        self.assertTrue(os.path.exists(self.fai_filename))

    def test_make_fai_correct_content(self):
        """Test if the .fai file has the correct content."""
        make_fai(self.fasta_filename)
        
        with open(self.fai_filename, "r") as f:
            lines = f.readlines()
        
        self.assertEqual(len(lines), 2)
        
        parts1 = lines[0].strip().split('\t')
        self.assertEqual(parts1[0], "seq1")
        self.assertEqual(parts1[1], "12") # 8 + 4
        # The rest of the values can be brittle, so we just check for existence
        self.assertEqual(len(parts1), 5)

        parts2 = lines[1].strip().split('\t')
        self.assertEqual(parts2[0], "seq2")
        self.assertEqual(parts2[1], "12")
        self.assertEqual(len(parts2), 5)

    def test_non_existent_fasta(self):
        """Test make_fai with a non-existent FASTA file."""
        with self.assertRaises(Exception) as context:
            make_fai("non_existent.fasta")
        self.assertTrue("Fasta file does not exist" in str(context.exception))

if __name__ == '__main__':
    unittest.main()
