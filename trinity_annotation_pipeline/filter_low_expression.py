#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Filter low-expressed genes based on salmon quantification results
"""
import os
import glob
import argparse
import pandas as pd
from Bio import SeqIO
from pathlib import Path


def filter_low_expression_genes(salmon_output_dir, threshold):
    """Filter low expression genes from salmon quant results"""
    salmon_folders = glob.glob(os.path.join(salmon_output_dir, "quant_*"))
    if not salmon_folders:
        print(f"Error: No quant_* folders found in {salmon_output_dir}")
        return None
    
    print(f"Found {len(salmon_folders)} salmon quantification folders")
    
    low_expression_genes = set()
    
    for folder in salmon_folders:
        try:
            quant_file = os.path.join(folder, "quant.sf")
            if not os.path.exists(quant_file):
                print(f"Warning: quant.sf not found in {folder}, skipping")
                continue
                
            df = pd.read_csv(quant_file, sep='\t')
            low_genes = set(df[df['NumReads'] < threshold]['Name'].tolist())
            low_expression_genes.update(low_genes)
            print(f"{os.path.basename(folder)}: {len(low_genes)} low-expressed genes (NumReads < {threshold})")
            
        except Exception as e:
            print(f"Error processing {folder}: {e}")
            continue
    
    return sorted(low_expression_genes)


def filter_fasta_by_expression(input_fasta, low_expression_genes, output_dir):
    """Filter FASTA file by expression level"""
    low_set = set(low_expression_genes)
    high_exp_sequences, low_exp_sequences = [], []
    
    high_output = os.path.join(output_dir, "trinity_high_expression.fasta")
    low_output = os.path.join(output_dir, "trinity_low_expression.fasta")
    
    try:
        for record in SeqIO.parse(input_fasta, "fasta"):
            (low_exp_sequences if record.id in low_set else high_exp_sequences).append(record)
        
        SeqIO.write(high_exp_sequences, high_output, "fasta")
        SeqIO.write(low_exp_sequences, low_output, "fasta")
            
    except Exception as e:
        print(f"Error processing FASTA file: {e}")
        return False
    
    print(f"\nProcessing summary:")
    print(f"Total sequences: {len(high_exp_sequences) + len(low_exp_sequences)}")
    print(f"High-expressed: {len(high_exp_sequences)}")
    print(f"Low-expressed: {len(low_exp_sequences)}")
    print(f"High-expressed output: {high_output}")
    print(f"Low-expressed output: {low_output}")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Filter low-expressed genes based on salmon quantification results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  filter_low_expression -s output_4_salmon -i input.fasta -o output_5_select -e 1
  
  filter_low_expression --salmon-dir output_4_salmon --input input.fasta --output output_5_select --expression-threshold 5
        """
    )
    parser.add_argument("-s", "--salmon-dir", default="output_4_salmon", help="Salmon output directory (default: output_4_salmon)")
    parser.add_argument("-i", "--input",  default="output_3_TransDecoder/trinity_multi_out.Trinity.fasta.transdecoder.pep", help="Input FASTA file (default: output_3_TransDecoder/trinity_multi_out.Trinity.fasta.transdecoder.pep)")
    parser.add_argument("-o", "--output", default="output_5_select", help="Output directory (default: output_5_select)")
    parser.add_argument("-e", "--expression-threshold", type=int, default=1, help="Expression threshold, genes with NumReads < threshold will be filtered (default: 1)")
    
    args = parser.parse_args()
    
    Path(args.output).mkdir(parents=True, exist_ok=True)
    
    print("Filtering low-expressed genes...")
    print(f"Salmon directory: {args.salmon_dir}")
    print(f"Input FASTA: {args.input}")
    print(f"Output directory: {args.output}")
    print(f"Expression threshold: NumReads < {args.expression_threshold}")
    print("-" * 50)
    
    low_genes = filter_low_expression_genes(args.salmon_dir, args.expression_threshold)
    
    if low_genes is None:
        print("Filtering process failed")
        return
    
    print(f"\nTotal low-expressed genes: {len(low_genes)} (NumReads < {args.expression_threshold})")
    
    if not os.path.exists(args.input):
        print(f"Error: Input FASTA file not found: {args.input}")
        return
    
    print("\nFiltering FASTA file...")
    if filter_fasta_by_expression(args.input, low_genes, args.output):
        print(f"\nProcessing completed! Results saved to {args.output}")


if __name__ == "__main__":
    main()

