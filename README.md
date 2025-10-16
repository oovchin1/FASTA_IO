# FASTA_IO

This repository provides tools for working with FASTA files, including indexing and efficient sequence access.

## Files

- [`FASTA_IO.py`](FASTA_IO.py): Contains the [`FASTA_IO`](FASTA_IO.py) class for reading, writing, and masking regions in indexed FASTA files.
- [`make_fai.py`](make_fai.py): Provides the [`make_fai`](make_fai.py) function to generate a `.fai` index file for a given FASTA file.

## Usage

### Creating a FASTA Index

To create a `.fai` index for your FASTA file:

```python
from make_fai import make_fai

make_fai("your_reference.fasta")
```

This will generate `your_reference.fasta.fai`.

### Accessing Sequences

Use the [`FASTA_IO`](FASTA_IO.py) class to read and modify sequences:

```python
from FASTA_IO import FASTA_IO

fasta = FASTA_IO("your_reference.fasta")
ids = fasta.get_sequence_IDs()
seq = fasta.read_in_section(ids[0], 0, 100)
fasta.hard_mask_region(ids[0], 0, 100)
```
## License

MIT License. See source files for details.