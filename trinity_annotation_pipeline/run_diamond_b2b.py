#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run DIAMOND bidirectional best hit analysis
"""
import os
import sys
import shutil
import argparse
import subprocess
from pathlib import Path


def check_and_download_uniprot():
    """Check if uniprots.pep exists, download if not"""
    if not os.path.exists("uniprots.pep"):
        print("UniProt protein file not found. Downloading...")
        try:
            # Try different import methods
            try:
                # When used as module
                from .download_uniprot import download_uniprot_pep
            except (ImportError, ValueError):
                try:
                    # When used as installed package
                    from trinity_annotation_pipeline.download_uniprot import download_uniprot_pep
                except ImportError:
                    # When used as script in the same directory
                    import importlib.util
                    spec = importlib.util.spec_from_file_location(
                        "download_uniprot", 
                        os.path.join(os.path.dirname(__file__), "download_uniprot.py")
                    )
                    download_uniprot_module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(download_uniprot_module)
                    download_uniprot_pep = download_uniprot_module.download_uniprot_pep
            
            if download_uniprot_pep(force_download=False):
                print("UniProt database downloaded successfully")
                return True
            else:
                print("Failed to download UniProt database")
                return False
        except Exception as e:
            print(f"Error downloading UniProt database: {e}")
            import traceback
            traceback.print_exc()
            return False
    else:
        print("UniProt protein file already exists")
        return True


def setup_working_directories(input_dir, output_dir):
    """Setup working directories and copy necessary files"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Define file paths
    low_expression_fasta = os.path.join(input_dir, "trinity_low_expression.fasta")
    uniprot_pep = "uniprots.pep"
    low_pep_dest = os.path.join(output_dir, "low.pep")
    uniprot_pep_dest = os.path.join(output_dir, "uniprots.pep")
    
    # Check if input files exist
    if not os.path.exists(low_expression_fasta):
        raise FileNotFoundError(f"Low expression FASTA file not found: {low_expression_fasta}")
    
    if not os.path.exists(uniprot_pep):
        raise FileNotFoundError(f"UniProt protein file not found: {uniprot_pep}")
    
    # Copy files
    shutil.copy2(low_expression_fasta, low_pep_dest)
    shutil.copy2(uniprot_pep, uniprot_pep_dest)
    
    print(f"Files copied:")
    print(f"  {low_expression_fasta} -> {low_pep_dest}")
    print(f"  {uniprot_pep} -> {uniprot_pep_dest}")
    
    return low_pep_dest, uniprot_pep_dest


def run_diamond_makedb(pep_file, db_name, threads=10):
    """Create DIAMOND database for protein file"""
    cmd = ["diamond", "makedb", "--in", pep_file, "-d", db_name, "--threads", str(threads)]
    
    print(f"Creating DIAMOND database for {pep_file}...")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Database created: {db_name}.dmnd")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Database creation failed: {e}")
        return False
    except FileNotFoundError:
        print("Error: diamond command not found, please ensure diamond is installed and in PATH")
        return False


def run_diamond_blastp(query, db, output, e_value="1e-5", threads=10):
    """Run DIAMOND blastp alignment"""
    cmd = [
        "diamond", "blastp", 
        "-q", query, "-d", db, "-o", output, 
        "--evalue", e_value, "--max-target-seqs", "1", 
        "--threads", str(threads), "--outfmt", "6"
    ]
    
    print(f"Running DIAMOND blastp...")
    print(f"Query: {query}")
    print(f"Database: {db}")
    print(f"Output: {output}")
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Alignment completed: {output}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Alignment failed: {e}")
        return False


def parse_diamond_result(result_file):
    """Parse DIAMOND result file, return gene pairs dictionary"""
    gene_pairs = {}
    
    try:
        with open(result_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    query_gene, target_gene = parts[0], parts[1]
                    gene_pairs[query_gene] = target_gene
                    
        print(f"Parsed {result_file}: found {len(gene_pairs)} gene pairs")
        return gene_pairs
    except Exception as e:
        print(f"Failed to parse DIAMOND result file: {e}")
        return {}


def find_reciprocal_best_hits(forward_pairs, reverse_pairs):
    """Find reciprocal best hit gene pairs"""
    reciprocal_pairs = []
    
    for query_gene, target_gene in forward_pairs.items():
        if target_gene in reverse_pairs and reverse_pairs[target_gene] == query_gene:
            reciprocal_pairs.append((query_gene, target_gene))
    
    print(f"Found {len(reciprocal_pairs)} reciprocal best hit gene pairs")
    return reciprocal_pairs


def save_reciprocal_pairs(reciprocal_pairs, output_file):
    """Save reciprocal best hit gene pairs to file"""
    try:
        with open(output_file, 'w') as f:
            f.write("Query_gene\tTarget_gene\n")
            for pair in reciprocal_pairs:
                f.write(f"{pair[0]}\t{pair[1]}\n")
        
        print(f"Reciprocal best hits saved to: {output_file}")
        return True
    except Exception as e:
        print(f"Failed to save reciprocal best hits: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Run DIAMOND bidirectional best hit analysis")
    parser.add_argument("-i", "--input", default="output_5_select",
                       help="Input directory (default: output_5_select)")
    parser.add_argument("-o", "--output", default="output_6_diamond_b2b",
                       help="Output directory (default: output_6_diamond_b2b)")
    parser.add_argument("-t", "--threads", type=int, default=10,
                       help="Number of threads for DIAMOND (default: 10)")
    parser.add_argument("-e", "--evalue", default="1e-5",
                       help="E-value threshold (default: 1e-5)")
    
    args = parser.parse_args()
    
    print("Starting DIAMOND bidirectional best hit analysis...")
    print(f"Input directory: {args.input}")
    print(f"Output directory: {args.output}")
    print(f"Threads: {args.threads}")
    print(f"E-value: {args.evalue}")
    print("-" * 50)
    
    # Check and download UniProt database if needed
    if not check_and_download_uniprot():
        print("Failed to obtain UniProt database, exiting")
        return
    
    # Setup working directories
    try:
        low_pep, uniprot_pep = setup_working_directories(args.input, args.output)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return
    
    # Create DIAMOND databases
    low_db = os.path.join(args.output, "low_db")
    uniprot_db = os.path.join(args.output, "uniprot_db")
    
    if not run_diamond_makedb(low_pep, low_db, args.threads):
        return
    
    if not run_diamond_makedb(uniprot_pep, uniprot_db, args.threads):
        return
    
    # Run bidirectional DIAMOND alignment
    forward_result = os.path.join(args.output, "low_vs_uniprot.tsv")
    if not run_diamond_blastp(low_pep, f"{uniprot_db}.dmnd", forward_result, args.evalue, args.threads):
        return
    
    reverse_result = os.path.join(args.output, "uniprot_vs_low.tsv")
    if not run_diamond_blastp(uniprot_pep, f"{low_db}.dmnd", reverse_result, args.evalue, args.threads):
        return
    
    # Parse results and find reciprocal best hits
    forward_pairs = parse_diamond_result(forward_result)
    reverse_pairs = parse_diamond_result(reverse_result)
    
    if not forward_pairs or not reverse_pairs:
        print("Error: unable to parse DIAMOND results")
        return
    
    reciprocal_pairs = find_reciprocal_best_hits(forward_pairs, reverse_pairs)
    
    if not reciprocal_pairs:
        print("Warning: no reciprocal best hits found")
        return
    
    # Save results
    reciprocal_output = os.path.join(args.output, "reciprocal_best_hits.tsv")
    if save_reciprocal_pairs(reciprocal_pairs, reciprocal_output):
        print(f"\nAnalysis completed! Results saved to: {args.output}")


if __name__ == "__main__":
    main()

