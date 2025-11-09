#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Filter TransDecoder results to keep high-expression and reciprocal best hit genes
"""
import os
import argparse
from Bio import SeqIO
from pathlib import Path


def read_reciprocal_hits(reciprocal_hits_file):
    """Read reciprocal best hits file and extract gene IDs"""
    gene_list = []
    
    try:
        with open(reciprocal_hits_file, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    gene_list.append(parts[0])
        
        print(f"Read {len(gene_list)} genes from {reciprocal_hits_file}")
        return gene_list
    except Exception as e:
        print(f"Failed to read reciprocal_best_hits.tsv: {e}")
        return []


def read_high_expression_genes(high_expression_fasta):
    """Read high expression FASTA file and extract gene IDs"""
    gene_list = []
    
    try:
        for record in SeqIO.parse(high_expression_fasta, "fasta"):
            gene_list.append(record.id)
        
        print(f"Read {len(gene_list)} high-expression genes from {high_expression_fasta}")
        return gene_list
    except Exception as e:
        print(f"Failed to read high expression FASTA: {e}")
        return []


def combine_gene_lists(reciprocal_genes, high_expression_genes):
    """Combine and deduplicate gene lists"""
    combined_genes = list(set(reciprocal_genes + high_expression_genes))
    print(f"Combined genes: {len(combined_genes)} (high-exp: {len(high_expression_genes)}, reciprocal: {len(reciprocal_genes)})")
    return combined_genes


def filter_fasta_file(input_file, output_file, gene_list, file_type):
    """Filter FASTA file (CDS or PEP) to keep only sequences in gene list"""
    gene_set = set(gene_list)
    kept_sequences = 0
    
    try:
        with open(output_file, 'w') as out_f:
            for record in SeqIO.parse(input_file, "fasta"):
                if record.id in gene_set:
                    SeqIO.write(record, out_f, "fasta")
                    kept_sequences += 1
        
        print(f"{file_type} file filtered: {kept_sequences} sequences kept")
        return True
    except Exception as e:
        print(f"Failed to filter {file_type} file: {e}")
        return False


def filter_bed_file(input_bed, output_bed, gene_list):
    """Filter BED file to keep only entries in gene list"""
    gene_set = set(gene_list)
    kept_entries = 0
    
    try:
        with open(input_bed, 'r') as in_f, open(output_bed, 'w') as out_f:
            for line in in_f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    attributes = parts[3].split(';')
                    if attributes:
                        gene_id_part = attributes[0]
                        if '=' in gene_id_part:
                            gene_id = gene_id_part.split('=')[1]
                            if gene_id in gene_set:
                                out_f.write(line)
                                kept_entries += 1
        
        print(f"BED file filtered: {kept_entries} entries kept")
        return True
    except Exception as e:
        print(f"Failed to filter BED file: {e}")
        return False


def filter_gff3_file(input_gff3, output_gff3, gene_list):
    """Filter GFF3 file to keep only genes in list and their features"""
    gene_set = set(gene_list)
    kept_entries = 0
    current_gene_in_list = False
    
    try:
        with open(input_gff3, 'r') as in_f, open(output_gff3, 'w') as out_f:
            for line in in_f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 9:
                    if parts[2] == "gene":
                        last_column = parts[8]
                        first_part = last_column.split(';')[0]
                        gene_id = first_part.split('~')[-1]
                        
                        if gene_id in gene_set:
                            current_gene_in_list = True
                            out_f.write(line)
                            kept_entries += 1
                        else:
                            current_gene_in_list = False
                    else:
                        if current_gene_in_list:
                            out_f.write(line)
                            kept_entries += 1
        
        print(f"GFF3 file filtered: {kept_entries} entries kept")
        return True
    except Exception as e:
        print(f"Failed to filter GFF3 file: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Filter TransDecoder results to keep high-expression and reciprocal best hit genes")
    
    parser.add_argument("-d", "--diamond-dir", default="output_6_diamond_b2b",
                       help="Directory with reciprocal best hits (default: output_6_diamond_b2b)")
    parser.add_argument("-s", "--select-dir", default="output_5_select",
                       help="Directory with high expression FASTA (default: output_5_select)")
    parser.add_argument("-t", "--transdecoder-dir", default="output_3_TransDecoder",
                       help="TransDecoder output directory (default: output_3_TransDecoder)")
    parser.add_argument("-o", "--output", default="output_7_height_add_uniprot",
                       help="Output directory (default: output_7_height_add_uniprot)")
    
    args = parser.parse_args()
    
    # Create output directory
    Path(args.output).mkdir(parents=True, exist_ok=True)
    
    print("Filtering TransDecoder result files...")
    print(f"Diamond directory: {args.diamond_dir}")
    print(f"Select directory: {args.select_dir}")
    print(f"TransDecoder directory: {args.transdecoder_dir}")
    print(f"Output directory: {args.output}")
    print("-" * 50)
    
    # Read reciprocal hits file
    reciprocal_hits_file = os.path.join(args.diamond_dir, "reciprocal_best_hits.tsv")
    if not os.path.exists(reciprocal_hits_file):
        print(f"Error: File not found: {reciprocal_hits_file}")
        return
    
    reciprocal_genes = read_reciprocal_hits(reciprocal_hits_file)
    if not reciprocal_genes:
        print("Error: Could not read reciprocal_best_hits.tsv")
        return
    
    # Read high expression genes
    high_expression_fasta = os.path.join(args.select_dir, "trinity_high_expression.fasta")
    if not os.path.exists(high_expression_fasta):
        print(f"Error: File not found: {high_expression_fasta}")
        return
    
    high_expression_genes = read_high_expression_genes(high_expression_fasta)
    if not high_expression_genes:
        print("Error: Could not read high expression FASTA file")
        return
    
    # Combine gene lists
    combined_genes = combine_gene_lists(reciprocal_genes, high_expression_genes)
    
    # Define input and output files
    base_name = "trinity_multi_out.Trinity.fasta.transdecoder"
    file_types = ["cds", "pep", "bed", "gff3"]
    
    input_files = {ft: os.path.join(args.transdecoder_dir, f"{base_name}.{ft}") for ft in file_types}
    output_files = {ft: os.path.join(args.output, f"{base_name}.filtered.{ft}") for ft in file_types}
    
    # Check input files exist
    for file_type, file_path in input_files.items():
        if not os.path.exists(file_path):
            print(f"Error: Input file not found: {file_path}")
            return
    
    # Filter files
    print("\nFiltering files...")
    
    # Filter CDS and PEP files
    if not filter_fasta_file(input_files["cds"], output_files["cds"], combined_genes, "CDS"):
        return
    
    if not filter_fasta_file(input_files["pep"], output_files["pep"], combined_genes, "PEP"):
        return
    
    # Filter BED and GFF3 files
    if not filter_bed_file(input_files["bed"], output_files["bed"], combined_genes):
        return
    
    if not filter_gff3_file(input_files["gff3"], output_files["gff3"], combined_genes):
        return
    
    print(f"\nAll files filtered successfully! Results saved to: {args.output}")


if __name__ == "__main__":
    main()

