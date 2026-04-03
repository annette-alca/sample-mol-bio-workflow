## 🧬 This is a molecular biology and bioinformatics workflow that is run from one Jupyter notebook.
- Automating the analysis of multiple Sanger sequencing files (limit of 50 set by NCBI).
- DNA Sanger sequencing data to a CSV file of the top Blast matches, XML files of 50 matches are also saved.

### Who it is for
Bench scientists, students, molecular biology enthusiasts who want to automate batch identification of unknown DNA sequences from Sanger sequencing data. Requires internet access for NCBI BLAST queries.

### What it does
1. **Extract** - Reads multiple .ab1 sequencing files and combines them into a single FASTA file
2. **BLAST** - Searches each sequence against NCBI's nucleotide database (with rate limiting)
3. **Summarize** - Produces a CSV table of top BLAST hits with accession numbers, identity %, and NCBI links

An example run with 10 sequences takes around 8 minutes. Helper functions are in [src/ab1_helpers.py](src/ab1_helpers.py).

### Project structure
```
seq_data/              <- Place your .ab1 files in a subfolder here
  old_ab1s/            <- Example: 10 sample .ab1 files
python_outputs/        <- Results are saved here (auto-created)
  old_ab1s/            <- Example: FASTA, XML, and CSV output
src/
  ab1_helpers.py       <- Helper functions used by the notebooks
tests/
  test_ab1_helpers.py  <- Unit tests
```

### Notebooks
- **ab1_blast_workflow_template.ipynb** - Blank template. Duplicate this and fill in your subfolder and filename to run with your own data.
- **ab1_blast_workflow_example.ipynb** - Pre-run example with 10 .ab1 files from [seq_data/old_ab1s/](seq_data/old_ab1s/) and results in [python_outputs/old_ab1s/](python_outputs/old_ab1s/).

Tested using Python 3.13 in a Jupyter notebook in VSCode.

## Installation

### Using uv (recommended)
```bash
uv sync
```

### Using pip
```bash
pip install -r requirements.txt
```

## Testing

A series of tests that validate the functions in ab1_helpers.py. Think of them as positive and negative controls that verify the code works correctly.

Run all tests (includes NCBI network calls, takes 1-2 minutes):
```bash
uv run pytest tests/ -v
```

Skip slow NCBI tests for a faster run:
```bash
uv run pytest tests/ -m "not slow" -v
```



-Annette
