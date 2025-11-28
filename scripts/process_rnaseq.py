#!/usr/bin/env python3
"""
PACE RNA-seq Data Processing Module

Process RNA-seq data for integration with ABC model:
1. Calculate gene expression (TPM/FPKM)
2. Identify enhancer RNAs (eRNAs)
3. Filter genes by expression level

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd


def read_gene_annotation(gtf_file):
    """Read gene annotation from GTF file"""
    print(f"Reading gene annotation: {gtf_file}")
    
    genes = {}
    
    if gtf_file.endswith('.gz'):
        f = gzip.open(gtf_file, 'rt')
    else:
        f = open(gtf_file, 'r')
    
    for line in f:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        
        chrom, source, feature, start, end, score, strand, frame, attributes = fields
        
        if feature != 'gene':
            continue
        
        # Parse attributes
        attrs = {}
        for attr in attributes.strip().split(';'):
            attr = attr.strip()
            if not attr:
                continue
            if ' "' in attr:
                key, value = attr.split(' "', 1)
                value = value.rstrip('"')
                attrs[key.strip()] = value
        
        gene_id = attrs.get('gene_id', '')
        gene_name = attrs.get('gene_name', gene_id)
        gene_type = attrs.get('gene_biotype', attrs.get('gene_type', 'unknown'))
        
        genes[gene_id] = {
            'gene_id': gene_id,
            'gene_name': gene_name,
            'chrom': chrom,
            'start': int(start),
            'end': int(end),
            'strand': strand,
            'gene_type': gene_type,
            'length': int(end) - int(start) + 1
        }
    
    f.close()
    print(f"  Loaded {len(genes)} genes")
    return genes


def read_count_matrix(count_file, genes):
    """Read gene count matrix from featureCounts or htseq-count output"""
    print(f"Reading count matrix: {count_file}")
    
    # Try to detect format
    with open(count_file, 'r') as f:
        first_line = f.readline()
        if first_line.startswith('#'):
            # featureCounts format
            header = f.readline().strip().split('\t')
            sample_cols = header[6:]  # Columns after Geneid, Chr, Start, End, Strand, Length
            skip_cols = 6
        else:
            # Simple format: gene_id \t count
            header = first_line.strip().split('\t')
            if len(header) == 2:
                sample_cols = ['sample']
                skip_cols = 0
            else:
                sample_cols = header[1:]
                skip_cols = 0
    
    # Read counts
    counts = {}
    gene_lengths = {}
    
    with open(count_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('Geneid'):
                continue
            
            fields = line.strip().split('\t')
            gene_id = fields[0]
            
            if skip_cols == 6:  # featureCounts format
                gene_lengths[gene_id] = int(fields[5])
                count_values = [float(x) for x in fields[6:]]
            else:
                count_values = [float(x) for x in fields[1:]]
                if gene_id in genes:
                    gene_lengths[gene_id] = genes[gene_id]['length']
                else:
                    gene_lengths[gene_id] = 1000  # Default length
            
            counts[gene_id] = count_values
    
    # Convert to DataFrame
    df = pd.DataFrame.from_dict(counts, orient='index', columns=sample_cols)
    
    print(f"  Loaded counts for {len(df)} genes, {len(sample_cols)} samples")
    return df, gene_lengths


def calculate_tpm(counts_df, gene_lengths):
    """Calculate TPM (Transcripts Per Million)"""
    print("Calculating TPM values...")
    
    # Get gene lengths in same order as counts
    lengths = np.array([gene_lengths.get(gene_id, 1000) for gene_id in counts_df.index])
    
    # RPK (Reads Per Kilobase)
    rpk = counts_df.values / (lengths[:, np.newaxis] / 1000)
    
    # TPM
    scaling_factor = rpk.sum(axis=0) / 1e6
    tpm = rpk / scaling_factor
    
    tpm_df = pd.DataFrame(tpm, index=counts_df.index, columns=counts_df.columns)
    
    print(f"  TPM calculation complete")
    return tpm_df


def calculate_fpkm(counts_df, gene_lengths, library_sizes=None):
    """Calculate FPKM (Fragments Per Kilobase per Million)"""
    print("Calculating FPKM values...")
    
    lengths = np.array([gene_lengths.get(gene_id, 1000) for gene_id in counts_df.index])
    
    if library_sizes is None:
        library_sizes = counts_df.sum(axis=0).values
    
    # FPKM = (count * 10^9) / (length * total_reads)
    fpkm = (counts_df.values * 1e9) / (lengths[:, np.newaxis] * library_sizes)
    
    fpkm_df = pd.DataFrame(fpkm, index=counts_df.index, columns=counts_df.columns)
    
    print(f"  FPKM calculation complete")
    return fpkm_df


def identify_expressed_genes(expression_df, min_tpm=1.0, min_samples=1):
    """Identify expressed genes based on TPM threshold"""
    print(f"Identifying expressed genes (TPM >= {min_tpm} in >= {min_samples} samples)...")
    
    # Count samples where gene is expressed
    expressed_mask = expression_df >= min_tpm
    sample_counts = expressed_mask.sum(axis=1)
    
    expressed_genes = expression_df.index[sample_counts >= min_samples].tolist()
    
    print(f"  Found {len(expressed_genes)} expressed genes")
    return expressed_genes


def identify_enhancer_rna(bam_file, enhancer_regions, output_file, min_reads=5):
    """
    Identify enhancer RNAs from RNA-seq BAM file
    
    This requires samtools and bedtools to be installed
    """
    import subprocess
    
    print(f"Identifying enhancer RNAs from: {bam_file}")
    
    # Count reads in enhancer regions using bedtools
    cmd = f"bedtools coverage -a {enhancer_regions} -b {bam_file} -counts"
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        return None
    
    # Parse results
    erna_counts = []
    for line in result.stdout.strip().split('\n'):
        fields = line.split('\t')
        chrom, start, end, name = fields[:4]
        count = int(fields[-1])
        
        erna_counts.append({
            'chrom': chrom,
            'start': int(start),
            'end': int(end),
            'name': name,
            'eRNA_count': count,
            'has_eRNA': count >= min_reads
        })
    
    erna_df = pd.DataFrame(erna_counts)
    erna_df.to_csv(output_file, sep='\t', index=False)
    
    n_with_erna = erna_df['has_eRNA'].sum()
    print(f"  Found {n_with_erna} enhancers with detectable eRNA")
    
    return erna_df


def create_gene_expression_file(expression_df, genes, output_file, expressed_genes=None):
    """Create gene expression file for PACE"""
    print(f"Creating gene expression file: {output_file}")
    
    records = []
    for gene_id in expression_df.index:
        if gene_id not in genes:
            continue
        
        gene = genes[gene_id]
        
        # Calculate mean expression across samples
        mean_tpm = expression_df.loc[gene_id].mean()
        
        # Check if expressed
        is_expressed = gene_id in expressed_genes if expressed_genes else mean_tpm >= 1.0
        
        records.append({
            'gene_id': gene_id,
            'gene_name': gene['gene_name'],
            'chrom': gene['chrom'],
            'start': gene['start'],
            'end': gene['end'],
            'strand': gene['strand'],
            'gene_type': gene['gene_type'],
            'mean_TPM': mean_tpm,
            'is_expressed': is_expressed
        })
    
    df = pd.DataFrame(records)
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"  Wrote expression data for {len(df)} genes")
    return df


def main():
    parser = argparse.ArgumentParser(
        description='PACE RNA-seq Data Processing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Process gene counts and calculate TPM
    python process_rnaseq.py \\
        --counts gene_counts.txt \\
        --gtf annotation.gtf.gz \\
        --output_dir rnaseq_processed/
    
    # With enhancer RNA identification
    python process_rnaseq.py \\
        --counts gene_counts.txt \\
        --gtf annotation.gtf.gz \\
        --bam aligned.bam \\
        --enhancers candidate_regions.bed \\
        --output_dir rnaseq_processed/
        """
    )
    
    parser.add_argument('--counts', required=True, 
                        help='Gene count matrix (featureCounts or htseq-count format)')
    parser.add_argument('--gtf', required=True, 
                        help='GTF annotation file')
    parser.add_argument('--output_dir', required=True, 
                        help='Output directory')
    parser.add_argument('--bam', 
                        help='RNA-seq BAM file for eRNA detection (optional)')
    parser.add_argument('--enhancers', 
                        help='Candidate enhancer regions BED file (for eRNA detection)')
    parser.add_argument('--min_tpm', type=float, default=1.0,
                        help='Minimum TPM for expressed genes (default: 1.0)')
    parser.add_argument('--min_samples', type=int, default=1,
                        help='Minimum samples with expression (default: 1)')
    parser.add_argument('--min_erna_reads', type=int, default=5,
                        help='Minimum reads for eRNA detection (default: 5)')
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("=" * 60)
    print("PACE RNA-seq Data Processing")
    print("=" * 60)
    
    # Read gene annotation
    genes = read_gene_annotation(args.gtf)
    
    # Read count matrix
    counts_df, gene_lengths = read_count_matrix(args.counts, genes)
    
    # Calculate TPM
    tpm_df = calculate_tpm(counts_df, gene_lengths)
    tpm_file = os.path.join(args.output_dir, 'gene_expression_TPM.tsv')
    tpm_df.to_csv(tpm_file, sep='\t')
    print(f"  Saved TPM matrix: {tpm_file}")
    
    # Calculate FPKM
    fpkm_df = calculate_fpkm(counts_df, gene_lengths)
    fpkm_file = os.path.join(args.output_dir, 'gene_expression_FPKM.tsv')
    fpkm_df.to_csv(fpkm_file, sep='\t')
    print(f"  Saved FPKM matrix: {fpkm_file}")
    
    # Identify expressed genes
    expressed_genes = identify_expressed_genes(tpm_df, args.min_tpm, args.min_samples)
    expressed_file = os.path.join(args.output_dir, 'expressed_genes.txt')
    with open(expressed_file, 'w') as f:
        for gene_id in expressed_genes:
            f.write(f"{gene_id}\n")
    print(f"  Saved expressed gene list: {expressed_file}")
    
    # Create gene expression file for PACE
    gene_expr_file = os.path.join(args.output_dir, 'gene_expression_summary.tsv')
    create_gene_expression_file(tpm_df, genes, gene_expr_file, expressed_genes)
    
    # Identify enhancer RNAs (if BAM provided)
    if args.bam and args.enhancers:
        erna_file = os.path.join(args.output_dir, 'enhancer_RNA.tsv')
        identify_enhancer_rna(args.bam, args.enhancers, erna_file, args.min_erna_reads)
    
    print("\n" + "=" * 60)
    print("Processing complete!")
    print("=" * 60)
    print(f"\nOutput files:")
    print(f"  TPM matrix:        {tpm_file}")
    print(f"  FPKM matrix:       {fpkm_file}")
    print(f"  Expressed genes:   {expressed_file}")
    print(f"  Gene summary:      {gene_expr_file}")
    if args.bam and args.enhancers:
        print(f"  Enhancer RNA:      {erna_file}")


if __name__ == '__main__':
    main()
