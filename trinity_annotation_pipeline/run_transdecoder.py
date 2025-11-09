#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Call TransDecoder to identify ORF regions and predict coding sequences
"""
import os
import subprocess
import argparse
from pathlib import Path


# 1. 运行TransDecoder以识别开放阅读框（ORFs）并预测编码序列
def run_transdecoder(input_fasta, output_dir, threads=4, min_protein_length=50):
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if not os.path.exists(input_fasta):
        print(f"Error: Input file does not exist: {input_fasta}")
        return False

    input_fasta = os.path.abspath(input_fasta)
    output_dir = os.path.abspath(output_dir)
    
    try:
        # 第一步: 识别ORF区域
        print("Step 1: Use TransDecoder.LongOrfs to identify ORF regions...")

        longorfs_cmd = f"TransDecoder.LongOrfs -t {input_fasta} -m {min_protein_length} -O {output_dir}"
        result = subprocess.run(
            longorfs_cmd, 
            shell=True,
            cwd=output_dir,
            capture_output=True,
            text=True,
            check=True
        )
        
        print("ORF identification successfully completed")
        if result.stdout:
            print(f"Standard output: {result.stdout}")
        if result.stderr:
            print(f"Standard error: {result.stderr}")
        
        # 第二步: 预测编码序列
        print("Step 2: Use TransDecoder.Predict to predict the encoding sequence...")
        predict_cmd = f"TransDecoder.Predict -t {input_fasta} --single_best_only --cpu {threads} -O {output_dir}"

        result = subprocess.run(
            predict_cmd,
            shell=True,
            cwd=output_dir,
            capture_output=True,
            text=True,
            check=True
        )
        
        print("Sequence prediction successfully completed.")
        if result.stdout:
            print(f"Standard output: {result.stdout}")
        if result.stderr:
            print(f"Standard error: {result.stderr}")
        
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"TransDecoder execution failed, exit code: {e.returncode}")
        print(f"Error Output: {e.stderr}")
        print(f"Standard Output: {e.stdout}")
        print("Try using different parameter formats...")
    
    except FileNotFoundError:
        print("The TransDecoder command was not found. Please ensure TransDecoder is installed and added to your PATH.")
        return False
    except Exception as e:
        print(f"Unexpected error: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Run TransDecoder to identify ORFs and predict coding sequences.")
    parser.add_argument("-i", "--input", default="output_2_trinity/trinity_multi_out.Trinity.fasta", help="Input Trinity assembly fasta file (default: output_2_trinity/trinity_multi_out.Trinity.fasta)")
    parser.add_argument("-o", "--output", default="output_3_TransDecoder", help="TransDecoder Output Directory (Default: output_3_TransDecoder)")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (Default: 4)")
    parser.add_argument("-min", "--min_protein_length", type=int, default=50, help="Minimum protein length for TransDecoder.LongOrfs (default: 50)")
    
    args = parser.parse_args()
    
    success = run_transdecoder(args.input, args.output, args.threads, args.min_protein_length)
    
    if success:
        print("TransDecoder analysis completed successfully")
    else:
        print("TransDecoder analysis failure")
        exit(1)


if __name__ == "__main__":
    main()

