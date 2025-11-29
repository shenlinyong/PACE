#!/usr/bin/env python3
"""
Generate example/test data for PACE
Creates a small synthetic dataset for testing the pipeline

Author: Linyong Shen @ Northwest A&F University
"""

import random
import gzip
import os

# Set random seed for reproducibility
random.seed(42)

# Configuration
OUTPUT_DIR = "example"
CHROM_NAME = "chr1"
CHROM_LENGTH = 5000000  # 5 Mb
NUM_GENES = 50
NUM_ATAC_READS = 100000
NUM_H3K27AC_READS = 80000

def generate_random_sequence(length):
    """Generate random DNA sequence"""
    return ''.join(random.choices('ACGT', k=length))

def generate_fasta(output_file):
    """Generate a small reference genome FASTA"""
    print(f"Generating reference genome: {output_file}")
    
    with gzip.open(output_file, 'wt') as f:
        f.write(f">{CHROM_NAME}\n")
        seq = generate_random_sequence(CHROM_LENGTH)
        # Write in 80-character lines
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")
    
    print(f"  Created {CHROM_NAME}: {CHROM_LENGTH:,} bp")

def generate_genes():
    """Generate random gene positions"""
    genes = []
    gene_id = 1
    pos = 10000  # Start position
    
    while pos < CHROM_LENGTH - 50000 and len(genes) < NUM_GENES:
        # Gene length: 1kb - 50kb
        gene_length = random.randint(1000, 50000)
        strand = random.choice(['+', '-'])
        
        # Number of exons: 1-10
        num_exons = random.randint(1, min(10, gene_length // 500))
        
        gene_start = pos
        gene_end = pos + gene_length
        
        # Generate exons
        exons = []
        exon_pos = gene_start
        for i in range(num_exons):
            remaining = gene_end - exon_pos
            if remaining < 200:
                break
            max_exon_len = min(500, remaining // 2)
            if max_exon_len < 100:
                break
            exon_length = random.randint(100, max_exon_len)
            exons.append((exon_pos, exon_pos + exon_length))
            exon_pos += exon_length + random.randint(100, min(2000, remaining - exon_length))  # Intron
        
        if exons:
            genes.append({
                'gene_id': f'ENSPIG{gene_id:011d}',
                'gene_name': f'GENE{gene_id}',
                'chrom': CHROM_NAME,
                'start': gene_start,
                'end': gene_end,
                'strand': strand,
                'exons': exons,
                'gene_type': 'protein_coding'
            })
            gene_id += 1
        
        # Gap between genes
        pos = gene_end + random.randint(5000, 50000)
    
    return genes

def generate_gtf(genes, output_file):
    """Generate GTF annotation file"""
    print(f"Generating GTF annotation: {output_file}")
    
    with gzip.open(output_file, 'wt') as f:
        f.write('##gtf-version 3\n')
        f.write(f'##genome-build Example Genome\n')
        
        for gene in genes:
            # Gene entry
            attrs = f'gene_id "{gene["gene_id"]}"; gene_name "{gene["gene_name"]}"; gene_biotype "{gene["gene_type"]}";'
            f.write(f'{gene["chrom"]}\tPACE\tgene\t{gene["start"]+1}\t{gene["end"]}\t.\t{gene["strand"]}\t.\t{attrs}\n')
            
            # Transcript entry
            transcript_id = gene["gene_id"].replace("ENSPIG", "ENSPIT")
            attrs = f'gene_id "{gene["gene_id"]}"; transcript_id "{transcript_id}"; gene_name "{gene["gene_name"]}"; gene_biotype "{gene["gene_type"]}";'
            f.write(f'{gene["chrom"]}\tPACE\ttranscript\t{gene["start"]+1}\t{gene["end"]}\t.\t{gene["strand"]}\t.\t{attrs}\n')
            
            # Exon entries
            for i, (exon_start, exon_end) in enumerate(gene["exons"], 1):
                attrs = f'gene_id "{gene["gene_id"]}"; transcript_id "{transcript_id}"; exon_number "{i}"; gene_name "{gene["gene_name"]}"; gene_biotype "{gene["gene_type"]}";'
                f.write(f'{gene["chrom"]}\tPACE\texon\t{exon_start+1}\t{exon_end}\t.\t{gene["strand"]}\t.\t{attrs}\n')
    
    print(f"  Created {len(genes)} genes")

def generate_peaks(genes):
    """Generate peak positions (enhancers and promoters)"""
    peaks = []
    
    # Add promoter peaks (near TSS)
    for gene in genes:
        if gene["strand"] == '+':
            tss = gene["start"]
        else:
            tss = gene["end"]
        
        # Promoter peak
        peak_center = tss + random.randint(-200, 200)
        peak_center = max(100, min(peak_center, CHROM_LENGTH - 100))
        peaks.append({
            'center': peak_center,
            'strength': random.uniform(0.7, 1.0),
            'type': 'promoter'
        })
        
        # Add enhancer peaks (5-500kb from TSS)
        num_enhancers = random.randint(1, 5)
        for _ in range(num_enhancers):
            distance = random.randint(5000, 500000) * random.choice([-1, 1])
            enh_center = tss + distance
            if 100 < enh_center < CHROM_LENGTH - 100:
                peaks.append({
                    'center': enh_center,
                    'strength': random.uniform(0.3, 0.8),
                    'type': 'enhancer'
                })
    
    # Add some random background peaks
    for _ in range(50):
        peaks.append({
            'center': random.randint(1000, CHROM_LENGTH - 1000),
            'strength': random.uniform(0.1, 0.3),
            'type': 'background'
        })
    
    return peaks

def generate_tagalign(peaks, output_file, num_reads, read_length=50):
    """Generate tagAlign file (simulated ATAC-seq reads)"""
    print(f"Generating tagAlign: {output_file}")
    
    reads = []
    
    for peak in peaks:
        # Number of reads proportional to peak strength
        n_reads = int(num_reads * peak['strength'] / len(peaks))
        
        for _ in range(n_reads):
            # Position follows normal distribution around peak center
            pos = int(random.gauss(peak['center'], 75))
            if pos < 0 or pos + read_length > CHROM_LENGTH:
                continue
            
            strand = random.choice(['+', '-'])
            if strand == '+':
                start = pos
                end = pos + read_length
            else:
                start = pos - read_length
                end = pos
                if start < 0:
                    continue
            
            reads.append((CHROM_NAME, start, end, 'N', 1000, strand))
    
    # Add some background reads
    for _ in range(num_reads // 10):
        start = random.randint(0, CHROM_LENGTH - read_length)
        end = start + read_length
        strand = random.choice(['+', '-'])
        reads.append((CHROM_NAME, start, end, 'N', 1000, strand))
    
    # Sort by position
    reads.sort(key=lambda x: x[1])
    
    # Write to file
    with gzip.open(output_file, 'wt') as f:
        for read in reads:
            f.write(f"{read[0]}\t{read[1]}\t{read[2]}\t{read[3]}\t{read[4]}\t{read[5]}\n")
    
    print(f"  Created {len(reads):,} reads")

def generate_chrom_sizes(output_file):
    """Generate chromosome sizes file"""
    print(f"Generating chrom.sizes: {output_file}")
    with open(output_file, 'w') as f:
        f.write(f"{CHROM_NAME}\t{CHROM_LENGTH}\n")

def main():
    os.makedirs(f"{OUTPUT_DIR}/data", exist_ok=True)
    os.makedirs(f"{OUTPUT_DIR}/reference", exist_ok=True)
    
    print("=" * 60)
    print("PACE Example Data Generator")
    print("=" * 60)
    print(f"\nGenerating example data for testing...")
    print(f"Chromosome: {CHROM_NAME}, Length: {CHROM_LENGTH:,} bp\n")
    
    # Generate reference genome
    generate_fasta(f"{OUTPUT_DIR}/reference/example_genome.fa.gz")
    
    # Generate genes
    genes = generate_genes()
    
    # Generate GTF
    generate_gtf(genes, f"{OUTPUT_DIR}/reference/example_annotation.gtf.gz")
    
    # Generate peaks
    peaks = generate_peaks(genes)
    
    # Generate ATAC-seq tagAlign
    generate_tagalign(peaks, f"{OUTPUT_DIR}/data/example_ATAC.tagAlign.gz", NUM_ATAC_READS)
    
    # Generate H3K27ac tagAlign (similar but different)
    h3k27ac_peaks = [p for p in peaks if p['type'] != 'background']
    generate_tagalign(h3k27ac_peaks, f"{OUTPUT_DIR}/data/example_H3K27ac.tagAlign.gz", NUM_H3K27AC_READS)
    
    # Generate chromosome sizes
    generate_chrom_sizes(f"{OUTPUT_DIR}/reference/example.chrom.sizes")
    
    print("\n" + "=" * 60)
    print("Example data generation complete!")
    print("=" * 60)
    print(f"\nGenerated files:")
    print(f"  Reference genome: {OUTPUT_DIR}/reference/example_genome.fa.gz")
    print(f"  GTF annotation:   {OUTPUT_DIR}/reference/example_annotation.gtf.gz")
    print(f"  Chrom sizes:      {OUTPUT_DIR}/reference/example.chrom.sizes")
    print(f"  ATAC-seq data:    {OUTPUT_DIR}/data/example_ATAC.tagAlign.gz")
    print(f"  H3K27ac data:     {OUTPUT_DIR}/data/example_H3K27ac.tagAlign.gz")

if __name__ == '__main__':
    main()
