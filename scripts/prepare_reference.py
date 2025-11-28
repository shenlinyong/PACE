#!/usr/bin/env python3
"""
PACE Reference File Preparation Tool

Generate all required reference files for PACE from GTF and FASTA files.

Usage:
    python prepare_reference.py --gtf annotation.gtf.gz --fasta genome.fa.gz --output_dir reference/my_species

This script generates:
    1. genome.chrom.sizes - Chromosome sizes file
    2. genes.bed - Gene annotation file
    3. genes.TSS500bp.bed - TSS region file
    4. blacklist.bed - Empty blacklist file (update manually if needed)

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import gzip
import os
import subprocess
import sys


def run_command(cmd, description=""):
    """Run shell command"""
    print(f"  Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {description}")
        print(f"  stderr: {result.stderr}")
        sys.exit(1)
    return result.stdout


def generate_chrom_sizes(fasta_file, output_file):
    """Generate chromosome sizes file from FASTA"""
    print("\n[1/4] Generating chromosome sizes file...")
    
    fai_file = fasta_file + '.fai'
    
    if not os.path.exists(fai_file):
        print(f"  Creating FASTA index: {fai_file}")
        if fasta_file.endswith('.gz'):
            print("  Note: FASTA file is compressed, using alternative method...")
            generate_chrom_sizes_from_gz(fasta_file, output_file)
            return
        else:
            run_command(f"samtools faidx {fasta_file}", "Failed to create FASTA index")
    
    run_command(f"cut -f1,2 {fai_file} > {output_file}", "Failed to generate chromosome sizes")
    
    with open(output_file) as f:
        line_count = sum(1 for _ in f)
    print(f"  Done: Found {line_count} chromosomes/scaffolds")


def generate_chrom_sizes_from_gz(fasta_gz, output_file):
    """Generate chromosome sizes from gzipped FASTA"""
    chrom_sizes = {}
    current_chrom = None
    current_size = 0
    
    print(f"  Parsing compressed FASTA file...")
    
    with gzip.open(fasta_gz, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chrom:
                    chrom_sizes[current_chrom] = current_size
                current_chrom = line[1:].split()[0]
                current_size = 0
            else:
                current_size += len(line)
        
        if current_chrom:
            chrom_sizes[current_chrom] = current_size
    
    with open(output_file, 'w') as out:
        for chrom, size in sorted(chrom_sizes.items()):
            out.write(f"{chrom}\t{size}\n")
    
    print(f"  Done: Found {len(chrom_sizes)} chromosomes/scaffolds")


def parse_gtf_for_genes(gtf_file, output_file, allowed_types, add_chr_prefix=False):
    """Extract gene information from GTF file"""
    print("\n[2/4] Generating gene BED file from GTF...")
    
    def parse_attributes(attr_string):
        attributes = {}
        for attr in attr_string.strip().split(';'):
            attr = attr.strip()
            if not attr:
                continue
            if ' "' in attr:
                key, value = attr.split(' "', 1)
                value = value.rstrip('"')
            elif '=' in attr:
                key, value = attr.split('=', 1)
            else:
                continue
            attributes[key.strip()] = value
        return attributes
    
    genes = {}
    
    if gtf_file.endswith('.gz'):
        f = gzip.open(gtf_file, 'rt')
    else:
        f = open(gtf_file, 'r')
    
    line_count = 0
    for line in f:
        if line.startswith('#'):
            continue
        
        line_count += 1
        if line_count % 100000 == 0:
            print(f"  Processed {line_count} lines...")
        
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        
        chrom, source, feature, start, end, score, strand, frame, attributes = fields
        
        if feature != 'gene':
            continue
        
        attrs = parse_attributes(attributes)
        
        gene_id = attrs.get('gene_id', '')
        gene_name = attrs.get('gene_name', attrs.get('gene_symbol', gene_id))
        gene_type = attrs.get('gene_biotype', attrs.get('gene_type', 'unknown'))
        
        if allowed_types and gene_type not in allowed_types:
            continue
        
        if add_chr_prefix and not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        
        if '_' in chrom or 'Un' in chrom or 'random' in chrom.lower():
            continue
        
        start = int(start) - 1  # BED is 0-based
        end = int(end)
        
        gene_key = f"{chrom}:{gene_id}"
        
        if gene_key not in genes:
            genes[gene_key] = {
                'chrom': chrom, 'start': start, 'end': end,
                'gene_name': gene_name, 'gene_id': gene_id,
                'strand': strand, 'gene_type': gene_type
            }
    
    f.close()
    
    with open(output_file, 'w') as out:
        for gene_key in sorted(genes.keys()):
            gene = genes[gene_key]
            name = f"{gene['gene_name']};{gene['gene_id']}"
            out.write(f"{gene['chrom']}\t{gene['start']}\t{gene['end']}\t{name}\t0\t{gene['strand']}\t{gene['gene_id']}\t{gene['gene_type']}\n")
    
    print(f"  Done: Found {len(genes)} genes")


def generate_tss_bed(genes_file, chrom_sizes_file, output_file, slop=500):
    """Generate TSS region BED file"""
    print("\n[3/4] Generating TSS region BED file...")
    
    chrom_sizes = {}
    with open(chrom_sizes_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chrom_sizes[parts[0]] = int(parts[1])
    
    gene_count = 0
    
    with open(genes_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
            
            chrom, start, end, name, score, strand = fields[:6]
            start, end = int(start), int(end)
            
            tss = start if strand == '+' else end
            
            tss_start = max(0, tss - slop)
            tss_end = tss + slop
            
            if chrom in chrom_sizes:
                tss_end = min(tss_end, chrom_sizes[chrom])
            
            outfile.write(f"{chrom}\t{tss_start}\t{tss_end}\t{name}\t{score}\t{strand}\n")
            gene_count += 1
    
    print(f"  Done: Processed {gene_count} genes (TSS ± {slop}bp)")


def create_empty_blacklist(output_file):
    """Create empty blacklist file"""
    print("\n[4/4] Creating blacklist file...")
    with open(output_file, 'w') as f:
        f.write("# Blacklist regions file\n")
        f.write("# Add regions to exclude here\n")
        f.write("# Format: chr\tstart\tend\n")
    print("  Done: Created empty blacklist file (update manually if needed)")


def main():
    parser = argparse.ArgumentParser(
        description='PACE Reference File Preparation Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python prepare_reference.py \\
        --gtf Sus_scrofa.Sscrofa11.1.111.gtf.gz \\
        --fasta Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz \\
        --output_dir reference/pig
    
    # With species name and gene types
    python prepare_reference.py \\
        --gtf annotation.gtf.gz \\
        --fasta genome.fa \\
        --output_dir reference/cattle \\
        --species_name bosTau9 \\
        --gene_types protein_coding,lncRNA
        """
    )
    
    parser.add_argument('--gtf', required=True, help='GTF annotation file (supports gzip)')
    parser.add_argument('--fasta', required=True, help='Reference genome FASTA file (supports gzip)')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--species_name', default='genome', help='Species name for file naming')
    parser.add_argument('--gene_types', default='protein_coding', 
                        help='Gene types to include, comma-separated')
    parser.add_argument('--tss_slop', type=int, default=500, help='TSS extension distance')
    parser.add_argument('--add_chr_prefix', action='store_true', help='Add chr prefix to chromosomes')
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("=" * 60)
    print("PACE Reference File Preparation Tool")
    print("=" * 60)
    print(f"\nInput files:")
    print(f"  GTF: {args.gtf}")
    print(f"  FASTA: {args.fasta}")
    print(f"\nOutput directory: {args.output_dir}")
    
    chrom_sizes_file = os.path.join(args.output_dir, f"{args.species_name}.chrom.sizes")
    genes_file = os.path.join(args.output_dir, "genes.bed")
    tss_file = os.path.join(args.output_dir, f"genes.TSS{args.tss_slop}bp.bed")
    blacklist_file = os.path.join(args.output_dir, "blacklist.bed")
    
    generate_chrom_sizes(args.fasta, chrom_sizes_file)
    
    allowed_types = set(args.gene_types.split(',')) if args.gene_types else None
    parse_gtf_for_genes(args.gtf, genes_file, allowed_types, args.add_chr_prefix)
    
    generate_tss_bed(genes_file, chrom_sizes_file, tss_file, args.tss_slop)
    
    create_empty_blacklist(blacklist_file)
    
    print("\n" + "=" * 60)
    print("Generated files:")
    print("=" * 60)
    print(f"  Chromosome sizes: {chrom_sizes_file}")
    print(f"  Gene annotation:  {genes_file}")
    print(f"  TSS regions:      {tss_file}")
    print(f"  Blacklist:        {blacklist_file}")
    
    print("\n" + "=" * 60)
    print("Next steps:")
    print("=" * 60)
    print(f"""
1. Verify generated files

2. Update blacklist if needed: {blacklist_file}

3. Update config/config.yaml:
   
   ref:
     chrom_sizes: "{chrom_sizes_file}"
     regions_blocklist: "{blacklist_file}"
     genes: "{genes_file}"
     genome_tss: "{tss_file}"

4. Prepare sample configuration: config/config_biosamples.tsv

5. Run PACE:
   snakemake --cores 8
""")


if __name__ == '__main__':
    main()
