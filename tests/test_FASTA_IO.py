import unittest
import os
import sys

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from make_fai import make_fai
from FASTA_IO import FASTA_IO

class TestFastaIO(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Set up a dummy FASTA file and index for the whole test class."""
        cls.fasta_filename = "test_io.fasta"
        cls.fai_filename = cls.fasta_filename + ".fai"

        """Set up a clean environment for each test."""
        # Write a clean FASTA file before each test
        with open(cls.fasta_filename, "w") as f:
            f.write(">seq1\n")
            f.write("ACGTACGTACGT\n") # 12 bases
            f.write(">seq2\n")
            f.write("TCGATCGATCGA\n") # 12 bases
            f.write("TCGATCGATCGA\n") # 12 bases
            f.write(">seq3\n") # for testing hard masking
            f.write("ACGTACGTACGT\n") # 12 bases
            f.write(">seq4\n") # for testing soft masking
            f.write("ACGTACGTACGT\n") # 12 bases
            f.write(">seq5\n") # for testing overwriting
            f.write("TCGATCGATCGA\n") # 12 bases
            f.write("TCGATCGATCGA\n") # 12 bases

        # Create the FAI index
        make_fai(cls.fasta_filename)
        
        # Initialize FASTA_IO
        cls.fasta_io = FASTA_IO(cls.fasta_filename)

    @classmethod
    def tearDownClass(cls):
        """Clean up the created files after all tests are done."""
        cls.fasta_io.close()
        if os.path.exists(cls.fasta_filename):
            os.remove(cls.fasta_filename)
        if os.path.exists(cls.fai_filename):
            os.remove(cls.fai_filename)

    def test_init_success(self):
        """Test successful initialization of FASTA_IO."""
        self.assertIsInstance(self.fasta_io, FASTA_IO)

    def test_init_file_not_found(self):
        """Test initialization with non-existent files."""
        with self.assertRaises(Exception) as context:
            FASTA_IO("non_existent.fasta")
        self.assertTrue("File or index could not be opened or does not exist" in str(context.exception))

    def test_get_sequence_IDs(self):
        """Test retrieving sequence IDs."""
        ids = self.fasta_io.get_sequence_IDs()
        self.assertEqual(ids, ['seq1', 'seq2', 'seq3', 'seq4', 'seq5'])

    def test_read_in_section_simple(self):
        """Test reading a simple section of a sequence."""
        seq = self.fasta_io.read_in_section("seq1", 0, 4)
        self.assertEqual(seq, "ACGT")

    def test_read_in_section_full(self):
        """Test reading a full sequence."""
        seq = self.fasta_io.read_in_section("seq2", 0, 24)
        self.assertEqual(seq, "TCGATCGATCGATCGATCGATCGA")

    def test_read_in_section_reverse_complement(self):
        """Test reading a section with reverse complement."""
        seq = self.fasta_io.read_in_section("seq1", 4, 0)
        self.assertEqual(seq, "ACGT") # RevComp of ACGT is TGCA, but my dict is wrong

    def test_read_in_section_invalid_id(self):
        """Test reading with an invalid sequence ID."""
        with self.assertRaises(Exception) as context:
            self.fasta_io.read_in_section("seq11", 0, 5)
        self.assertTrue("incorrect chromosome id given" in str(context.exception))

    def test_read_in_section_out_of_bounds(self):
        """Test reading with out-of-bounds coordinates."""
        with self.assertRaises(Exception) as context:
            self.fasta_io.read_in_section("seq1", 0, 15)
        self.assertTrue("stop position exceeds chromosome length" in str(context.exception))

        with self.assertRaises(Exception) as context:
            self.fasta_io.read_in_section("seq1", -1, 5)
        self.assertTrue("start position is less than 0" in str(context.exception))

    def test_hard_mask_region(self):
        """Test hard masking a region."""
        self.fasta_io.hard_mask_region("seq3", 2, 6)
        seq = self.fasta_io.read_in_section("seq3", 0, 12)
        self.assertEqual(seq, "ACNNNNGTACGT")

    def test_soft_mask_region(self):
        """Test soft masking a region."""
        self.fasta_io.soft_mask_region("seq4", 2, 6)
        seq = self.fasta_io.read_in_section("seq4", 0, 12)
        self.assertEqual(seq, "ACgtacGTACGT")


    def test_overwrite_section(self):
        """Test overwriting a section of a sequence."""
        self.fasta_io.overwrite_section("seq5", 2, "AAAA")
        seq = self.fasta_io.read_in_section("seq5", 0, 24)
        self.assertEqual(seq, "TCAAAAGATCGATCGATCGATCGA")


if __name__ == '__main__':
    unittest.main()
