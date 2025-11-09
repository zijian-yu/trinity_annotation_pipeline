#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run transcript quantification with Salmon
"""
import os
import glob
import argparse
import subprocess
import concurrent.futures
from pathlib import Path


def build_salmon_index(transcript_file, index_dir, kmer_size=31):
    """Build Salmon index"""
    if not os.path.exists(transcript_file):
        raise FileNotFoundError(f"Transcript file not found: {transcript_file}")
    
    Path(index_dir).mkdir(parents=True, exist_ok=True)
    
    cmd = ["salmon", "index", "-t", transcript_file, "-i", index_dir, "-k", str(kmer_size)]
    
    print("Building Salmon index...")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("Index built successfully!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Index building failed: {e.stderr}")
        return False
    except FileNotFoundError:
        print("Error: salmon command not found, please ensure salmon is installed and in PATH")
        return False


def find_fastq_files(fastq_dir):
    """Find paired fastq files in directory"""
    r1_files = glob.glob(os.path.join(fastq_dir, "*_clean_R1.fq.gz"))
    samples = {}
    for r1_file in r1_files:
        r2_file = r1_file.replace("_clean_R1.fq.gz", "_clean_R2.fq.gz")
        if os.path.exists(r2_file):
            sample_name = os.path.basename(r1_file).replace("_clean_R1.fq.gz", "")
            samples[sample_name] = {'r1': r1_file, 'r2': r2_file}
    print(f"Found {len(samples)} paired samples")
    return samples


def run_salmon_quant(index_dir, r1_file, r2_file, output_dir, threads=20):
    """Run Salmon quantification for a single sample"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    cmd = [
        "salmon", "quant", "-i", index_dir, "-l", "A",
        "-1", r1_file, "-2", r2_file, "-o", output_dir,
        "--gcBias", "--validateMappings", "-p", str(threads)
    ]
    
    sample_name = os.path.basename(output_dir)
    print(f"Processing sample: {sample_name}")
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Sample {sample_name} completed")
        return True, sample_name
    except subprocess.CalledProcessError as e:
        print(f"Sample {sample_name} failed: {e.stderr}")
        return False, sample_name


def main():
    parser = argparse.ArgumentParser(description="Salmon transcript quantification pipeline")
    parser.add_argument("-c", "--transcript", default="output_3_TransDecoder/trinity_multi_out.Trinity.fasta.transdecoder.cds",
                       help="Transcript pep file path")
    parser.add_argument("-x", "--index", default="output_4_salmon/transcript_index",
                       help="Salmon index directory")
    parser.add_argument("-f", "--fastq-dir", default="output_1_cat_and_fastq",
                       help="Fastq files directory")
    parser.add_argument("-o", "--output", default="output_4_salmon",
                       help="Output directory")
    parser.add_argument("-k", "--kmer", type=int, default=31,
                       help="k-mer size")
    parser.add_argument("-t", "--threads", type=int, default=20,
                       help="Threads per salmon job")
    parser.add_argument("-m", "--max-workers", type=int, default=15,
                       help="Maximum parallel jobs")
    parser.add_argument("-s", "--skip-index", action="store_true",
                       help="Skip index building step")
    
    args = parser.parse_args()
    
    # Check input directory exists
    if not os.path.exists(args.fastq_dir):
        print(f"Error: fastq directory not found: {args.fastq_dir}")
        return
    
    # Build salmon index if needed
    if not args.skip_index:
        if not build_salmon_index(args.transcript, args.index, args.kmer):
            return
    else:
        print("Skipping index building step")
    
    # Check index exists
    if not os.path.exists(args.index):
        print(f"Error: index directory not found: {args.index}")
        print("Please run without -s/--skip-index to build index")
        return
    
    # Find fastq files
    print(f"Searching for fastq files in {args.fastq_dir}...")
    samples = find_fastq_files(args.fastq_dir)
    
    if not samples:
        print("No paired fastq files found")
        return
    
    # Prepare salmon quantification tasks
    tasks = []
    for sample_name, reads in samples.items():
        sample_output_dir = os.path.join(args.output, f"quant_{sample_name}")
        tasks.append({
            'index_dir': args.index,
            'r1_file': reads['r1'],
            'r2_file': reads['r2'],
            'output_dir': sample_output_dir,
            'threads': args.threads
        })
    
    # Run salmon quantification in parallel
    print(f"Processing {len(tasks)} samples, max parallel: {args.max_workers}")
    
    successful_samples, failed_samples = [], []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        future_to_task = {executor.submit(run_salmon_quant, **task): task for task in tasks}
        
        for future in concurrent.futures.as_completed(future_to_task):
            sample_name = os.path.basename(future_to_task[future]['output_dir']).replace('quant_', '')
            try:
                success, sample = future.result()
                successful_samples.append(sample) if success else failed_samples.append(sample)
            except Exception as e:
                print(f"Sample {sample_name} exception: {e}")
                failed_samples.append(sample_name)
    
    # Print summary
    print("\n" + "="*50)
    print("Analysis summary:")
    print(f"Successful samples: {len(successful_samples)}")
    print(f"Failed samples: {len(failed_samples)}")
    
    if successful_samples:
        print("Successful samples:")
        for sample in successful_samples:
            print(f"  - {sample}")
    
    if failed_samples:
        print("Failed samples:")
        for sample in failed_samples:
            print(f"  - {sample}")


if __name__ == "__main__":
    main()

