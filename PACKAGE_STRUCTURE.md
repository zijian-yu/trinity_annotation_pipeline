# Package Structure

## Directory Layout

```
trinity_annotation_pipeline/
├── setup.py                 # Package setup configuration
├── requirements.txt         # Python dependencies
├── MANIFEST.in             # Files to include in package
├── README.md               # Main documentation
├── INSTALL.md              # Installation instructions
├── example_usage.py        # Example usage script
└── trinity_annotation_pipeline/
    ├── __init__.py         # Package initialization
    ├── download_uniprot.py      # Step 0: Download UniProt database
    ├── cat_and_fastq.py         # Step 1: Quality control and merge FASTQ
    ├── run_trinity.py           # Step 2: Trinity assembly
    ├── run_transdecoder.py      # Step 3: TransDecoder ORF prediction
    ├── run_salmon.py            # Step 4: Salmon quantification
    ├── filter_low_expression.py # Step 5: Filter low expression genes
    ├── run_diamond_b2b.py       # Step 6: DIAMOND bidirectional best hit
    ├── new_cds_pep_gff.py       # Step 7: Generate final filtered files
    └── pipeline.py              # Main pipeline runner
```

## Command Line Tools

After installation, the following commands are available:

1. `trinity-download-uniprot` - Download UniProt database
2. `trinity-cat-fastq` - Quality control and merge FASTQ files
3. `trinity-assembly` - Run Trinity transcriptome assembly
4. `trinity-transdecoder` - Run TransDecoder ORF prediction
5. `trinity-salmon` - Run Salmon transcript quantification
6. `trinity-filter-expression` - Filter low-expressed genes
7. `trinity-diamond` - Run DIAMOND bidirectional best hit
8. `trinity-filter-results` - Filter and generate final files
9. `trinity-pipeline` - Run complete pipeline

## Module Functions

Each module can be imported and used programmatically:

```python
from trinity_annotation_pipeline import (
    download_uniprot,
    cat_and_fastq,
    run_trinity,
    run_transdecoder,
    run_salmon,
    filter_low_expression,
    run_diamond_b2b,
    new_cds_pep_gff,
    pipeline
)
```

## Installation

```bash
pip install .
```

## Usage

### Command Line

```bash
# Run complete pipeline
trinity-pipeline -i input_raw_fastq -s Paca

# Run individual steps
trinity-cat-fastq -i input_raw_fastq -s Paca
trinity-assembly -s Paca
# etc.
```

### Python API

```python
from trinity_annotation_pipeline import pipeline

# Run complete pipeline programmatically
pipeline.main()
```

See `example_usage.py` for detailed examples.

