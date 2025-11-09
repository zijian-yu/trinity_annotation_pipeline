# Installation Guide

## Quick Start

1. **Install the package:**
```bash
cd trinity_annotation_pipeline
pip install .
```

2. **Install dependencies:**
```bash
pip install -r requirements.txt
```

3. **Install external tools:**
   - fastp: https://github.com/OpenGene/fastp
   - Trinity: https://github.com/trinityrnaseq/trinityrnaseq
   - TransDecoder: https://github.com/TransDecoder/TransDecoder
   - Salmon: https://github.com/COMBINE-lab/salmon
   - DIAMOND: https://github.com/bbuchfink/diamond

## Development Installation

For development, install in editable mode:

```bash
pip install -e .
```

## Verify Installation

After installation, verify that all commands are available:

```bash
trinity-download-uniprot --help
trinity-cat-fastq --help
trinity-assembly --help
trinity-transdecoder --help
trinity-salmon --help
trinity-filter-expression --help
trinity-diamond --help
trinity-filter-results --help
trinity-pipeline --help
```

## Requirements

### Python Dependencies

- Python >= 3.7
- tqdm >= 4.65.0
- pandas >= 1.5.0
- biopython >= 1.79

### External Tools

All external tools must be installed and available in your PATH:

- **fastp** - Quality control tool
- **Trinity** - Transcriptome assembler
- **TransDecoder** - ORF predictor
- **Salmon** - Transcript quantifier
- **DIAMOND** - Protein aligner

## Troubleshooting

### Command not found

If commands are not found after installation, make sure:
1. The installation completed successfully
2. Your Python environment is activated
3. The bin/scripts directory is in your PATH

### Import errors

If you encounter import errors:
1. Make sure all dependencies are installed: `pip install -r requirements.txt`
2. Reinstall the package: `pip install --force-reinstall .`

### External tool errors

If external tools are not found:
1. Check that tools are installed: `which fastp`, `which Trinity`, etc.
2. Make sure tools are in your PATH
3. Check tool versions are compatible

