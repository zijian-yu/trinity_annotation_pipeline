#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Main pipeline script to run all steps of the Trinity annotation pipeline
"""
import os
import sys
import argparse
from pathlib import Path

# Import all pipeline modules
from . import (
    download_uniprot,
    cat_and_fastq,
    run_trinity,
    run_transdecoder,
    run_salmon,
    filter_low_expression,
    run_diamond_b2b,
    new_cds_pep_gff
)


def run_step(step_name, step_func, *args, **kwargs):
    """Run a pipeline step with error handling"""
    print("\n" + "="*70)
    print(f"Running Step: {step_name}")
    print("="*70)
    try:
        result = step_func(*args, **kwargs)
        if result is False:
            print(f"ERROR: {step_name} failed!")
            return False
        print(f"SUCCESS: {step_name} completed!")
        return True
    except Exception as e:
        print(f"ERROR: {step_name} failed with exception: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Run complete Trinity annotation pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline with default parameters
  trinity-pipeline -i input_raw_fastq -s Paca
  
  # Run with custom output directories
  trinity-pipeline -i input_raw_fastq -s Paca --work-dir ./work
  
  # Run with custom threads and memory
  trinity-pipeline -i input_raw_fastq -s Paca -t 20 --max-workers 2 --memory 100G
        """
    )
    
    # Required arguments
    parser.add_argument("-i", "--input", required=True,
                       help="Input directory with raw FASTQ files")
    parser.add_argument("-s", "--species", required=True,
                       help="Species name for output files")
    
    # Optional arguments for output directories
    parser.add_argument("--work-dir", default=".",
                       help="Working directory for all outputs (default: current directory)")
    parser.add_argument("--output-1", default="output_1_cat_and_fastq",
                       help="Output directory for step 1 (default: output_1_cat_and_fastq)")
    parser.add_argument("--output-2", default="output_2_trinity",
                       help="Output directory for step 2 (default: output_2_trinity)")
    parser.add_argument("--output-3", default="output_3_TransDecoder",
                       help="Output directory for step 3 (default: output_3_TransDecoder)")
    parser.add_argument("--output-4", default="output_4_salmon",
                       help="Output directory for step 4 (default: output_4_salmon)")
    parser.add_argument("--output-5", default="output_5_select",
                       help="Output directory for step 5 (default: output_5_select)")
    parser.add_argument("--output-6", default="output_6_diamond_b2b",
                       help="Output directory for step 6 (default: output_6_diamond_b2b)")
    parser.add_argument("--output-7", default="output_7_height_add_uniprot",
                       help="Output directory for step 7 (default: output_7_height_add_uniprot)")
    
    # Pipeline control
    parser.add_argument("--skip-step", action="append", type=int, choices=range(0, 8),
                       help="Skip specific steps (0-7). Can be used multiple times.")
    parser.add_argument("--start-from", type=int, choices=range(0, 8), default=0,
                       help="Start from specific step (0-7, default: 0)")
    parser.add_argument("--stop-at", type=int, choices=range(1, 8),
                       help="Stop at specific step (1-7)")
    
    # Step 0: Download UniProt
    parser.add_argument("--force-download", action="store_true",
                       help="Force re-download of UniProt database")
    parser.add_argument("--uniprot-file", default="uniprots.pep",
                       help="UniProt database filename (default: uniprots.pep)")
    
    # Unified parameters
    parser.add_argument("-t", "--threads", type=int, default=20,
                       help="Number of threads for all steps (default: 20)")
    parser.add_argument("--max-workers", type=int, default=2,
                       help="Maximum parallel workers for parallel tasks (default: 15)")
    
    # Step 2: Trinity assembly
    parser.add_argument("--memory", default="50G",
                       help="Maximum memory for Trinity (default: 50G)")
    
    # Step 3: TransDecoder
    parser.add_argument("--min-protein-length", type=int, default=50,
                       help="Minimum protein length for TransDecoder (default: 50)")
    
    # Step 4: Salmon
    parser.add_argument("--salmon-kmer", type=int, default=31,
                       help="k-mer size for Salmon (default: 31)")
    parser.add_argument("--skip-salmon-index", action="store_true",
                       help="Skip Salmon index building")
    
    # Step 5: Filter expression
    parser.add_argument("--expression-threshold", type=int, default=1,
                       help="Expression threshold for filtering (default: 1)")
    
    # Step 6: DIAMOND
    parser.add_argument("--diamond-evalue", default="1e-5",
                       help="E-value threshold for DIAMOND (default: 1e-5)")
    
    args = parser.parse_args()
    
    # Set working directory
    work_dir = Path(args.work_dir).absolute()
    work_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(work_dir)
    
    # Normalize output paths
    output_dirs = {
        1: args.output_1,
        2: args.output_2,
        3: args.output_3,
        4: args.output_4,
        5: args.output_5,
        6: args.output_6,
        7: args.output_7,
    }
    
    # Determine which steps to skip
    skip_steps = set(args.skip_step or [])
    
    # Determine which steps to run
    start_step = args.start_from
    stop_step = args.stop_at if args.stop_at else 7
    
    print("\n" + "="*70)
    print("Trinity Annotation Pipeline")
    print("="*70)
    print(f"Input directory: {args.input}")
    print(f"Species name: {args.species}")
    print(f"Working directory: {work_dir}")
    print(f"Steps to run: {start_step} to {stop_step}")
    if skip_steps:
        print(f"Steps to skip: {sorted(skip_steps)}")
    print("="*70)
    
    # Step 0: Download UniProt database
    if 0 >= start_step and 0 <= stop_step and 0 not in skip_steps:
        if not run_step(
            "Step 0: Download UniProt database",
            download_uniprot.download_uniprot_pep,
            force_download=args.force_download,
            output_file=args.uniprot_file
        ):
            print("Pipeline failed at step 0")
            sys.exit(1)
    else:
        print(f"\nSkipping step 0 (download UniProt)")
    
    # Step 1: Quality control and merge FASTQ files
    if 1 >= start_step and 1 <= stop_step and 1 not in skip_steps:
        print("\n" + "="*70)
        print("Step 1: Quality control and merge FASTQ files")
        print("="*70)
        
        clean_files = cat_and_fastq.quality_control_parallel(
            args.input, output_dirs[1], args.max_workers, args.threads
        )
        
        if not clean_files:
            print("No clean files generated, exiting")
            sys.exit(1)
        
        final_r1, final_r2 = cat_and_fastq.merge_clean_fastq(
            clean_files, output_dirs[1], args.species
        )
        
        if not final_r1 or not final_r2:
            print("Pipeline failed at step 1 (merge)")
            sys.exit(1)
        
        print("SUCCESS: Step 1 completed!")
    else:
        print(f"\nSkipping step 1 (quality control)")
    
    # Step 2: Trinity assembly
    if 2 >= start_step and 2 <= stop_step and 2 not in skip_steps:
        if not run_step(
            "Step 2: Trinity transcriptome assembly",
            run_trinity.run_trinity_assembly,
            output_dirs[1],
            output_dirs[2],
            args.species,
            args.threads,
            args.memory
        ):
            print("Pipeline failed at step 2")
            sys.exit(1)
    else:
        print(f"\nSkipping step 2 (Trinity assembly)")
    
    # Step 3: TransDecoder
    if 3 >= start_step and 3 <= stop_step and 3 not in skip_steps:
        trinity_fasta = os.path.join(output_dirs[2], "trinity_multi_out.Trinity.fasta")
        if not run_step(
            "Step 3: TransDecoder ORF prediction",
            run_transdecoder.run_transdecoder,
            trinity_fasta,
            output_dirs[3],
            args.threads,
            args.min_protein_length
        ):
            print("Pipeline failed at step 3")
            sys.exit(1)
    else:
        print(f"\nSkipping step 3 (TransDecoder)")
    
    # Step 4: Salmon quantification
    if 4 >= start_step and 4 <= stop_step and 4 not in skip_steps:
        transcript_file = os.path.join(output_dirs[3], "trinity_multi_out.Trinity.fasta.transdecoder.cds")
        index_dir = os.path.join(output_dirs[4], "transcript_index")
        
        if not args.skip_salmon_index:
            if not run_step(
                "Step 4: Build Salmon index",
                run_salmon.build_salmon_index,
                transcript_file,
                index_dir,
                args.salmon_kmer
            ):
                print("Pipeline failed at step 4 (index building)")
                sys.exit(1)
        
        # Find fastq files and run quantification
        samples = run_salmon.find_fastq_files(output_dirs[1])
        if not samples:
            print("No paired fastq files found for Salmon quantification")
            sys.exit(1)
        
        tasks = []
        for sample_name, reads in samples.items():
            sample_output_dir = os.path.join(output_dirs[4], f"quant_{sample_name}")
            tasks.append({
                'index_dir': index_dir,
                'r1_file': reads['r1'],
                'r2_file': reads['r2'],
                'output_dir': sample_output_dir,
                'threads': args.threads
            })
        
        print(f"\nRunning Salmon quantification for {len(tasks)} samples...")
        successful_samples, failed_samples = [], []
        import concurrent.futures
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.max_workers) as executor:
            future_to_task = {executor.submit(run_salmon.run_salmon_quant, **task): task for task in tasks}
            
            for future in concurrent.futures.as_completed(future_to_task):
                sample_name = os.path.basename(future_to_task[future]['output_dir']).replace('quant_', '')
                try:
                    success, sample = future.result()
                    successful_samples.append(sample) if success else failed_samples.append(sample)
                except Exception as e:
                    print(f"Sample {sample_name} exception: {e}")
                    failed_samples.append(sample_name)
        
        if failed_samples:
            print(f"Warning: {len(failed_samples)} samples failed in Salmon quantification")
    else:
        print(f"\nSkipping step 4 (Salmon quantification)")
    
    # Step 5: Filter low expression
    if 5 >= start_step and 5 <= stop_step and 5 not in skip_steps:
        print("\n" + "="*70)
        print("Step 5: Filter low expression genes")
        print("="*70)
        
        input_pep = os.path.join(output_dirs[3], "trinity_multi_out.Trinity.fasta.transdecoder.pep")
        
        low_genes = filter_low_expression.filter_low_expression_genes(
            output_dirs[4], args.expression_threshold
        )
        
        if low_genes is None:
            print("Pipeline failed at step 5 (filter genes)")
            sys.exit(1)
        
        Path(output_dirs[5]).mkdir(parents=True, exist_ok=True)
        
        if not filter_low_expression.filter_fasta_by_expression(
            input_pep, low_genes, output_dirs[5]
        ):
            print("Pipeline failed at step 5")
            sys.exit(1)
        
        print("SUCCESS: Step 5 completed!")
    else:
        print(f"\nSkipping step 5 (filter expression)")
    
    # Step 6: DIAMOND bidirectional best hit
    if 6 >= start_step and 6 <= stop_step and 6 not in skip_steps:
        # Check and download UniProt if needed
        if not os.path.exists(args.uniprot_file):
            download_uniprot.download_uniprot_pep(force_download=False, output_file=args.uniprot_file)
        
        try:
            low_pep, uniprot_pep = run_diamond_b2b.setup_working_directories(
                output_dirs[5], output_dirs[6]
            )
        except FileNotFoundError as e:
            print(f"Error: {e}")
            sys.exit(1)
        
        low_db = os.path.join(output_dirs[6], "low_db")
        uniprot_db = os.path.join(output_dirs[6], "uniprot_db")
        
        if not run_step(
            "Step 6: Create DIAMOND databases",
            run_diamond_b2b.run_diamond_makedb,
            low_pep,
            low_db,
            args.threads
        ):
            print("Pipeline failed at step 6 (database creation)")
            sys.exit(1)
        
        if not run_step(
            "Step 6: Create UniProt DIAMOND database",
            run_diamond_b2b.run_diamond_makedb,
            uniprot_pep,
            uniprot_db,
            args.threads
        ):
            print("Pipeline failed at step 6 (UniProt database creation)")
            sys.exit(1)
        
        forward_result = os.path.join(output_dirs[6], "low_vs_uniprot.tsv")
        if not run_step(
            "Step 6: DIAMOND forward alignment",
            run_diamond_b2b.run_diamond_blastp,
            low_pep,
            f"{uniprot_db}.dmnd",
            forward_result,
            args.diamond_evalue,
            args.threads
        ):
            print("Pipeline failed at step 6 (forward alignment)")
            sys.exit(1)
        
        reverse_result = os.path.join(output_dirs[6], "uniprot_vs_low.tsv")
        if not run_step(
            "Step 6: DIAMOND reverse alignment",
            run_diamond_b2b.run_diamond_blastp,
            uniprot_pep,
            f"{low_db}.dmnd",
            reverse_result,
            args.diamond_evalue,
            args.threads
        ):
            print("Pipeline failed at step 6 (reverse alignment)")
            sys.exit(1)
        
        forward_pairs = run_diamond_b2b.parse_diamond_result(forward_result)
        reverse_pairs = run_diamond_b2b.parse_diamond_result(reverse_result)
        
        if not forward_pairs or not reverse_pairs:
            print("Error: unable to parse DIAMOND results")
            sys.exit(1)
        
        reciprocal_pairs = run_diamond_b2b.find_reciprocal_best_hits(forward_pairs, reverse_pairs)
        
        reciprocal_output = os.path.join(output_dirs[6], "reciprocal_best_hits.tsv")
        if not run_step(
            "Step 6: Save reciprocal best hits",
            run_diamond_b2b.save_reciprocal_pairs,
            reciprocal_pairs,
            reciprocal_output
        ):
            print("Pipeline failed at step 6 (save results)")
            sys.exit(1)
    else:
        print(f"\nSkipping step 6 (DIAMOND)")
    
    # Step 7: Filter and generate final files
    if 7 >= start_step and 7 <= stop_step and 7 not in skip_steps:
        reciprocal_hits_file = os.path.join(output_dirs[6], "reciprocal_best_hits.tsv")
        high_expression_fasta = os.path.join(output_dirs[5], "trinity_high_expression.fasta")
        
        if not os.path.exists(reciprocal_hits_file):
            print(f"Error: File not found: {reciprocal_hits_file}")
            sys.exit(1)
        
        if not os.path.exists(high_expression_fasta):
            print(f"Error: File not found: {high_expression_fasta}")
            sys.exit(1)
        
        reciprocal_genes = new_cds_pep_gff.read_reciprocal_hits(reciprocal_hits_file)
        high_expression_genes = new_cds_pep_gff.read_high_expression_genes(high_expression_fasta)
        combined_genes = new_cds_pep_gff.combine_gene_lists(reciprocal_genes, high_expression_genes)
        
        base_name = "trinity_multi_out.Trinity.fasta.transdecoder"
        file_types = ["cds", "pep", "bed", "gff3"]
        
        input_files = {ft: os.path.join(output_dirs[3], f"{base_name}.{ft}") for ft in file_types}
        output_files = {ft: os.path.join(output_dirs[7], f"{base_name}.filtered.{ft}") for ft in file_types}
        
        Path(output_dirs[7]).mkdir(parents=True, exist_ok=True)
        
        for file_type in ["cds", "pep"]:
            if not run_step(
                f"Step 7: Filter {file_type.upper()} file",
                new_cds_pep_gff.filter_fasta_file,
                input_files[file_type],
                output_files[file_type],
                combined_genes,
                file_type.upper()
            ):
                print(f"Pipeline failed at step 7 (filter {file_type})")
                sys.exit(1)
        
        if not run_step(
            "Step 7: Filter BED file",
            new_cds_pep_gff.filter_bed_file,
            input_files["bed"],
            output_files["bed"],
            combined_genes
        ):
            print("Pipeline failed at step 7 (filter bed)")
            sys.exit(1)
        
        if not run_step(
            "Step 7: Filter GFF3 file",
            new_cds_pep_gff.filter_gff3_file,
            input_files["gff3"],
            output_files["gff3"],
            combined_genes
        ):
            print("Pipeline failed at step 7 (filter gff3)")
            sys.exit(1)
    else:
        print(f"\nSkipping step 7 (filter results)")
    
    print("\n" + "="*70)
    print("Pipeline completed successfully!")
    print("="*70)
    print(f"Final results are in: {output_dirs[7]}")
    print("="*70)


if __name__ == "__main__":
    main()

