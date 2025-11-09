#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Merge all sequencing files and conduct quality control
"""
import os
import glob
import argparse
import subprocess
from pathlib import Path
import concurrent.futures


# 1. 对输入目录中的所有FASTQ文件进行并行质控
def quality_control_parallel(input_dir, output_dir="output_1_cat_and_fastq", max_workers=4, fastp_threads=1):
    # 确保输出目录存在
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    pattern = os.path.join(input_dir, "*.fastq.gz")
    fastq_files = glob.glob(pattern)
    
    if not fastq_files:
        print(f"No FASTQ files found in {input_dir}")
        return []
    
    # 按样本名分组（假设文件命名格式：samplename_1.fastq.gz, samplename_2.fastq.gz）
    samples = {}
    for fq_file in fastq_files:
        basename = os.path.basename(fq_file)
        if "_1.fastq.gz" in basename:
            sample_name = basename.replace("_1.fastq.gz", "")
            samples.setdefault(sample_name, {})["r1"] = fq_file
        elif "_2.fastq.gz" in basename:
            sample_name = basename.replace("_2.fastq.gz", "")
            samples.setdefault(sample_name, {})["r2"] = fq_file

    valid_samples = {name: paths for name, paths in samples.items() 
                    if "r1" in paths and "r2" in paths}
    
    if not valid_samples:
        print("No valid paired-end samples found")
        return []
    
    print(f"Found {len(valid_samples)} samples for quality control")
    
    def process_single_sample(sample_name, r1_file, r2_file):
        """处理单个样本的质控"""
        clean_r1 = os.path.join(output_dir, f"{sample_name}_clean_R1.fq.gz")
        clean_r2 = os.path.join(output_dir, f"{sample_name}_clean_R2.fq.gz")
        json_report = os.path.join(output_dir, f"{sample_name}_fastp.json")
        html_report = os.path.join(output_dir, f"{sample_name}_fastp.html")
        
        # fastp命令，包含线程参数
        cmd = [
            "fastp", "-i", r1_file, "-I", r2_file,
            "-o", clean_r1, "-O", clean_r2,
            "--detect_adapter_for_pe", "--correction",
            "--compression", "6", "-w", str(fastp_threads),
            "-j", json_report, "-h", html_report
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            print(f"Quality control completed for {sample_name}")
            return clean_r1, clean_r2
        except subprocess.CalledProcessError as e:
            print(f"fastp failed for {sample_name}: {e}")
            return None, None
    
    # 并行处理所有样本
    print(f"Running quality control with {max_workers} workers (each using {fastp_threads} threads)...")
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_sample = {
            executor.submit(process_single_sample, name, paths["r1"], paths["r2"]): name
            for name, paths in valid_samples.items()
        }

        for future in concurrent.futures.as_completed(future_to_sample):
            sample_name = future_to_sample[future]
            try:
                result = future.result()
                if result[0] and result[1]:
                    results.append(result)
            except Exception as e:
                print(f"Sample {sample_name} generated an exception: {e}")
    
    print(f"Quality control completed for {len(results)} samples")
    return results


# 2. 合并所有质控后的FASTQ文件
def merge_clean_fastq(clean_files, output_dir, species_name):
    if not clean_files:
        print("No clean files to merge")
        return None, None
    
    r1_files = [f[0] for f in clean_files]
    r2_files = [f[1] for f in clean_files]
    merged_r1 = os.path.join(output_dir, f"{species_name}_combined_R1_clean.fq.gz")
    merged_r2 = os.path.join(output_dir, f"{species_name}_combined_R2_clean.fq.gz")

    try:
        cmd = f"cat {' '.join(r1_files)} > {merged_r1}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"R1 files merged: {merged_r1}")
    except subprocess.CalledProcessError as e:
        print(f"R1 merge failed: {e}")
        return None, None

    try:
        cmd = f"cat {' '.join(r2_files)} > {merged_r2}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"R2 files merged: {merged_r2}")
    except subprocess.CalledProcessError as e:
        print(f"R2 merge failed: {e}")
        return None, None
    
    return merged_r1, merged_r2


def main():
    parser = argparse.ArgumentParser(description="Quality control and merge FASTQ files")
    parser.add_argument("-i", "--input", required=True, help="Input directory with raw FASTQ files")
    parser.add_argument("-s", "--species", required=True, help="Species name for output files")
    parser.add_argument("-o", "--output", default="output_1_cat_and_fastq", help="Output directory for results (default: output_1_cat_and_fastq)")
    parser.add_argument("-m", "--max_workers", type=int, default=2, help="Max workers for parallel QC (default: 2)")
    parser.add_argument("-t", "--fastp_threads", type=int, default=4, help="Threads for each fastp process (default: 4)")
    
    args = parser.parse_args()
    Path(args.output).mkdir(parents=True, exist_ok=True)
    
    # 第一步：并行质控
    print("Step 1: Running parallel quality control...")
    clean_files = quality_control_parallel(args.input, args.output, args.max_workers, args.fastp_threads)
    
    if not clean_files:
        print("No clean files generated, exiting")
        return
    
    # 第二步：合并质控后的文件
    print("Step 2: Merging clean FASTQ files...")
    final_r1, final_r2 = merge_clean_fastq(clean_files, args.output, args.species)
    
    if final_r1 and final_r2:
        print("Processing completed successfully!")
        print(f"Final files: {final_r1}, {final_r2}")
    else:
        print("Processing failed!")
        exit(1)


if __name__ == "__main__":
    main()

