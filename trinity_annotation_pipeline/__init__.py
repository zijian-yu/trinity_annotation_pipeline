"""
Trinity Annotation Pipeline
A comprehensive pipeline for transcriptome assembly, annotation, and quantification.
"""

__version__ = "1.0.0"
__author__ = "yuzijian"

from . import (
    download_uniprot,
    cat_and_fastq,
    run_trinity,
    run_transdecoder,
    run_salmon,
    filter_low_expression,
    run_diamond_b2b,
    new_cds_pep_gff,
    pipeline,
)

__all__ = [
    "download_uniprot",
    "cat_and_fastq",
    "run_trinity",
    "run_transdecoder",
    "run_salmon",
    "filter_low_expression",
    "run_diamond_b2b",
    "new_cds_pep_gff",
    "pipeline",
]

