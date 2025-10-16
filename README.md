# FASTA_IO

This repository provides tools for working with FASTA files, including indexing and efficient sequence access.

## Files

- `FASTA_IO.py`: Contains the `FASTA_IO` class for reading, writing, and masking regions in indexed FASTA files.
- `make_fai.py`: Provides the `make_fai` function to generate a `.fai` index file for a given FASTA file.
- `tests/`: Contains unit tests for the project.
  - `test_make_fai.py`: Unit tests for `make_fai.py`.
  - `test_FASTA_IO.py`: Unit tests for `FASTA_IO.py`.

## Usage

### Creating a FASTA Index

Before using `FASTA_IO`, you need to create a FASTA index file (`.fai`).

```python
from make_fai import make_fai

# This will generate `your_reference.fasta.fai`
make_fai("your_reference.fasta")
```

### Reading and Manipulating FASTA data

The `FASTA_IO` class provides an interface to read and write to an indexed FASTA file.

```python
from FASTA_IO import FASTA_IO

# Initialize with the path to your FASTA file
fasta = FASTA_IO("your_reference.fasta")

# Get a list of all sequence IDs
ids = fasta.get_sequence_IDs()
print(f"Sequence IDs: {ids}")

# Read a section of a sequence (e.g., first 100 bases of the first sequence)
seq = fasta.read_in_section(ids[0], 0, 100)
print(f"First 100 bases: {seq}")

# Hard mask a region with 'N'
fasta.hard_mask_region(ids[0], 10, 20)
masked_seq = fasta.read_in_section(ids[0], 0, 30)
print(f"After hard masking: {masked_seq}")

# Soft mask a region (convert to lowercase)
fasta.soft_mask_region(ids[0], 30, 40)
soft_masked_seq = fasta.read_in_section(ids[0], 20, 50)
print(f"After soft masking: {soft_masked_seq}")

# Overwrite a section with a new sequence
fasta.overwrite_section(ids[0], 50, "ACGT")
overwritten_seq = fasta.read_in_section(ids[0], 40, 60)
print(f"After overwriting: {overwritten_seq}")
```

## Testing

The project includes unit tests to ensure correctness. To run the tests, execute the test files from the root directory of the project:

```bash
python tests/test_make_fai.py
python tests/test_FASTA_IO.py
```

## License

MIT License. See source files for details.