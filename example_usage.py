#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example usage of Trinity Annotation Pipeline
This script demonstrates how to use the pipeline programmatically
"""

import os
from trinity_annotation_pipeline import (
    download_uniprot,
    cat_and_fastq,
    run_trinity,
    run_transdecoder,
    run_salmon,
    filter_low_expression,
    run_diamond_b2b,
    new_cds_pep_gff
)

def example_full_pipeline():
    """Example of running the full pipeline programmatically"""
    
    # Configuration
    input_dir = "input_raw_fastq"
    species_name = "Paca"
    work_dir = "."
    
    # Step 0: Download UniProt database
    print("Step 0: Downloading UniProt database...")
    download_uniprot.download_uniprot_pep(force_download=False)
    
    # Step 1: Quality control and merge FASTQ files
    print("\nStep 1: Quality control and merging FASTQ files...")
    output_1 = "output_1_cat_and_fastq"
    clean_files = cat_and_fastq.quality_control_parallel(
        input_dir, output_1, max_workers=2, fastp_threads=4
    )
    if clean_files:
        cat_and_fastq.merge_clean_fastq(clean_files, output_1, species_name)
    
    # Step 2: Trinity assembly
    print("\nStep 2: Running Trinity assembly...")
    output_2 = "output_2_trinity"
    run_trinity.run_trinity_assembly(
        output_1, output_2, species_name, threads=4, max_memory="50G"
    )
    
    # Step 3: TransDecoder
    print("\nStep 3: Running TransDecoder...")
    output_3 = "output_3_TransDecoder"
    trinity_fasta = os.path.join(output_2, "trinity_multi_out.Trinity.fasta")
    run_transdecoder.run_transdecoder(
        trinity_fasta, output_3, threads=4, min_protein_length=50
    )
    
    # Step 4: Salmon quantification
    print("\nStep 4: Running Salmon quantification...")
    output_4 = "output_4_salmon"
    transcript_file = os.path.join(output_3, "trinity_multi_out.Trinity.fasta.transdecoder.cds")
    index_dir = os.path.join(output_4, "transcript_index")
    
    # Build index
    run_salmon.build_salmon_index(transcript_file, index_dir, kmer_size=31)
    
    # Find samples and run quantification
    samples = run_salmon.find_fastq_files(output_1)
    for sample_name, reads in samples.items():
        sample_output_dir = os.path.join(output_4, f"quant_{sample_name}")
        run_salmon.run_salmon_quant(
            index_dir, reads['r1'], reads['r2'], sample_output_dir, threads=20
        )
    
    # Step 5: Filter low expression
    print("\nStep 5: Filtering low expression genes...")
    output_5 = "output_5_select"
    input_pep = os.path.join(output_3, "trinity_multi_out.Trinity.fasta.transdecoder.pep")
    
    low_genes = filter_low_expression.filter_low_expression_genes(
        output_4, threshold=1
    )
    if low_genes:
        filter_low_expression.filter_fasta_by_expression(
            input_pep, low_genes, output_5
        )
    
    # Step 6: DIAMOND bidirectional best hit
    print("\nStep 6: Running DIAMOND bidirectional best hit...")
    output_6 = "output_6_diamond_b2b"
    
    # Setup directories
    low_pep, uniprot_pep = run_diamond_b2b.setup_working_directories(
        output_5, output_6
    )
    
    # Create databases
    low_db = os.path.join(output_6, "low_db")
    uniprot_db = os.path.join(output_6, "uniprot_db")
    run_diamond_b2b.run_diamond_makedb(low_pep, low_db, threads=10)
    run_diamond_b2b.run_diamond_makedb(uniprot_pep, uniprot_db, threads=10)
    
    # Run alignments
    forward_result = os.path.join(output_6, "low_vs_uniprot.tsv")
    reverse_result = os.path.join(output_6, "uniprot_vs_low.tsv")
    run_diamond_b2b.run_diamond_blastp(
        low_pep, f"{uniprot_db}.dmnd", forward_result, e_value="1e-5", threads=10
    )
    run_diamond_b2b.run_diamond_blastp(
        uniprot_pep, f"{low_db}.dmnd", reverse_result, e_value="1e-5", threads=10
    )
    
    # Find reciprocal best hits
    forward_pairs = run_diamond_b2b.parse_diamond_result(forward_result)
    reverse_pairs = run_diamond_b2b.parse_diamond_result(reverse_result)
    reciprocal_pairs = run_diamond_b2b.find_reciprocal_best_hits(
        forward_pairs, reverse_pairs
    )
    
    # Save results
    reciprocal_output = os.path.join(output_6, "reciprocal_best_hits.tsv")
    run_diamond_b2b.save_reciprocal_pairs(reciprocal_pairs, reciprocal_output)
    
    # Step 7: Filter and generate final files
    print("\nStep 7: Filtering and generating final files...")
    output_7 = "output_7_height_add_uniprot"
    
    reciprocal_hits_file = os.path.join(output_6, "reciprocal_best_hits.tsv")
    high_expression_fasta = os.path.join(output_5, "trinity_high_expression.fasta")
    
    reciprocal_genes = new_cds_pep_gff.read_reciprocal_hits(reciprocal_hits_file)
    high_expression_genes = new_cds_pep_gff.read_high_expression_genes(
        high_expression_fasta
    )
    combined_genes = new_cds_pep_gff.combine_gene_lists(
        reciprocal_genes, high_expression_genes
    )
    
    # Filter files
    base_name = "trinity_multi_out.Trinity.fasta.transdecoder"
    file_types = ["cds", "pep", "bed", "gff3"]
    
    input_files = {
        ft: os.path.join(output_3, f"{base_name}.{ft}") for ft in file_types
    }
    output_files = {
        ft: os.path.join(output_7, f"{base_name}.filtered.{ft}") for ft in file_types
    }
    
    os.makedirs(output_7, exist_ok=True)
    
    # Filter CDS and PEP
    new_cds_pep_gff.filter_fasta_file(
        input_files["cds"], output_files["cds"], combined_genes, "CDS"
    )
    new_cds_pep_gff.filter_fasta_file(
        input_files["pep"], output_files["pep"], combined_genes, "PEP"
    )
    
    # Filter BED and GFF3
    new_cds_pep_gff.filter_bed_file(
        input_files["bed"], output_files["bed"], combined_genes
    )
    new_cds_pep_gff.filter_gff3_file(
        input_files["gff3"], output_files["gff3"], combined_genes
    )
    
    print("\nPipeline completed successfully!")
    print(f"Final results are in: {output_7}")


if __name__ == "__main__":
    # Note: This is an example script. 
    # For actual usage, use the command-line interface or the pipeline module.
    print("This is an example script showing programmatic usage.")
    print("For command-line usage, use: trinity-pipeline -i input -s species")
    print("\nTo run this example, uncomment the following line:")
    print("# example_full_pipeline()")

