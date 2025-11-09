#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Assemble the parameter-free transcriptome using Trinity.
"""
import os
import argparse
import subprocess
from pathlib import Path


# 1. 运行Trinity转录组组装
def run_trinity_assembly(input_dir, output_dir, species_name, threads=4, max_memory="50G"):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Define input file paths
    input_r1 = os.path.join(input_dir, f"{species_name}_combined_R1_clean.fq.gz")
    input_r2 = os.path.join(input_dir, f"{species_name}_combined_R2_clean.fq.gz")
    
    # Check if input files exist
    if not os.path.exists(input_r1) or not os.path.exists(input_r2):
        print(f"Error: Input files not found: {input_r1} or {input_r2}")
        return False
    
    # Define Trinity output directory
    trinity_output = os.path.join(output_dir, "trinity_multi_out")
    
    # Build Trinity command
    trinity_cmd = [
        "Trinity",
        "--seqType", "fq",
        "--left", input_r1,
        "--right", input_r2,
        "--CPU", str(threads),
        "--min_kmer_cov", "2",
        "--no_normalize_reads",
        "--full_cleanup",
        "--trimmomatic",
        "--jaccard_clip",
        "--output", trinity_output,
        "--max_memory", max_memory,
        "--normalize_max_read_cov", "50"
    ]
    
    # Run Trinity
    print(f"Starting Trinity assembly with {threads} threads and {max_memory} memory...")
    print(f"Input files: {input_r1}, {input_r2}")
    print(f"Output directory: {trinity_output}")
    
    try:
        subprocess.run(trinity_cmd, check=True)
        print(f"Trinity assembly completed successfully")
        print(f"Results saved to: {trinity_output}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Trinity assembly failed: {e}")
        return False
    except FileNotFoundError:
        print("Trinity command not found, please ensure Trinity is installed and in PATH")
        return False


def main():
    parser = argparse.ArgumentParser(description="Run Trinity transcriptome assembly")
    parser.add_argument("-i", "--input", default="output_1_cat_and_fastq", help="Input directory containing cleaned FASTQ files (default: output_1_cat_and_fastq)")
    parser.add_argument("-o", "--output", default="output_2_trinity", help="Output directory for Trinity results (default: output_2_trinity)")
    parser.add_argument("-s", "--species", required=True, help="Species name for file naming")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")
    parser.add_argument("-mem", "--memory", default="50G", help="Maximum memory to use (default: 50G)")
    
    args = parser.parse_args()
    
    # 运行 Trinity 无参组装
    success = run_trinity_assembly(args.input, args.output, args.species, args.threads, args.memory)
    
    if success:
        print("Transcriptome assembly completed successfully")
    else:
        print("Transcriptome assembly failed")
        exit(1)


if __name__ == "__main__":
    main()

