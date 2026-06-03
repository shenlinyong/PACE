#!/usr/bin/env python3
"""
Generate TSS region BED file from gene BED file

Usage:
    python prepare_tss_bed.py --genes genes.bed --chrom_sizes genome.chrom.sizes --output genes.TSS500bp.bed

Author: Linyong Shen @ Northwest A&F University
"""

import argparse


def main():
    parser = argparse.ArgumentParser(
        description='Generate TSS region BED file from gene BED file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python prepare_tss_bed.py --genes genes.bed --chrom_sizes genome.chrom.sizes --output genes.TSS500bp.bed
    python prepare_tss_bed.py --genes genes.bed --chrom_sizes genome.chrom.sizes --slop 1000 --output genes.TSS1kb.bed
        """
    )
    
    parser.add_argument('--genes', required=True, help='Input gene BED file')
    parser.add_argument('--chrom_sizes', required=True, help='Chromosome sizes file')
    parser.add_argument('--output', required=True, help='Output TSS BED file')
    parser.add_argument('--slop', type=int, default=500, help='TSS extension distance (default: 500)')
    
    args = parser.parse_args()
    
    # Read chromosome sizes
    chrom_sizes = {}
    print(f"Reading chromosome sizes: {args.chrom_sizes}")
    with open(args.chrom_sizes, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chrom_sizes[parts[0]] = int(parts[1])
    
    print(f"  Found {len(chrom_sizes)} chromosomes")
    
    # Process gene file
    print(f"Processing gene file: {args.genes}")
    gene_count = 0
    
    with open(args.genes, 'r') as infile, open(args.output, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
            
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3]
            score = fields[4]
            strand = fields[5]
            
            # Determine TSS position
            if strand == '+':
                tss = start
            else:
                tss = end
            
            # Calculate TSS region
            tss_start = max(0, tss - args.slop)
            tss_end = tss + args.slop
            
            # Ensure within chromosome boundaries
            if chrom in chrom_sizes:
                tss_end = min(tss_end, chrom_sizes[chrom])
            
            outfile.write(f"{chrom}\t{tss_start}\t{tss_end}\t{name}\t{score}\t{strand}\n")
            gene_count += 1
    
    print(f"Done! Processed {gene_count} genes")
    print(f"Output: {args.output}")
    print(f"TSS region: ±{args.slop}bp")


if __name__ == '__main__':
    main()
