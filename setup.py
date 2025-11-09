#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Setup script for Trinity Annotation Pipeline
"""
from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

# Read requirements
requirements_file = Path(__file__).parent / "requirements.txt"
if requirements_file.exists():
    with open(requirements_file, "r", encoding="utf-8") as f:
        requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]
else:
    requirements = [
        "tqdm>=4.65.0",
        "pandas>=1.5.0",
        "biopython>=1.79",
    ]

setup(
    name="trinity-annotation-pipeline",
    version="1.0.0",
    author="yuzijian",
    author_email="yuzijian1010@163.com",
    description="A comprehensive pipeline for transcriptome assembly, annotation, and quantification",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yuzijian/trinity-annotation-pipeline",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "trinity-download-uniprot=trinity_annotation_pipeline.download_uniprot:main",
            "trinity-cat-fastq=trinity_annotation_pipeline.cat_and_fastq:main",
            "trinity-assembly=trinity_annotation_pipeline.run_trinity:main",
            "trinity-transdecoder=trinity_annotation_pipeline.run_transdecoder:main",
            "trinity-salmon=trinity_annotation_pipeline.run_salmon:main",
            "trinity-filter-expression=trinity_annotation_pipeline.filter_low_expression:main",
            "trinity-diamond=trinity_annotation_pipeline.run_diamond_b2b:main",
            "trinity-filter-results=trinity_annotation_pipeline.new_cds_pep_gff:main",
            "trinity-pipeline=trinity_annotation_pipeline.pipeline:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)

