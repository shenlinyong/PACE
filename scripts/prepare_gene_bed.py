#!/usr/bin/env python3
"""
Generate gene BED file from GTF annotation

Usage:
    python prepare_gene_bed.py --gtf annotation.gtf --output genes.bed

Output format:
    chr  start  end  name;Ensembl_ID  score  strand  Ensembl_ID  gene_type

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import gzip
import sys


def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string"""
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


def main():
    parser = argparse.ArgumentParser(
        description='Generate gene BED file from GTF annotation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python prepare_gene_bed.py --gtf Sus_scrofa.gtf.gz --output genes.bed
    python prepare_gene_bed.py --gtf annotation.gtf --output genes.bed --gene_types protein_coding,lncRNA
        """
    )
    
    parser.add_argument('--gtf', required=True, help='Input GTF file (supports gzip)')
    parser.add_argument('--output', required=True, help='Output BED file')
    parser.add_argument('--gene_types', default='protein_coding', 
                        help='Gene types to include, comma-separated (default: protein_coding)')
    parser.add_argument('--feature_type', default='gene', 
                        help='GTF feature type (default: gene)')
    parser.add_argument('--add_chr_prefix', action='store_true',
                        help='Add chr prefix to chromosome names')
    parser.add_argument('--remove_chr_prefix', action='store_true',
                        help='Remove chr prefix from chromosome names')
    
    args = parser.parse_args()
    
    allowed_gene_types = set(args.gene_types.split(','))
    genes = {}
    
    if args.gtf.endswith('.gz'):
        f = gzip.open(args.gtf, 'rt')
    else:
        f = open(args.gtf, 'r')
    
    print(f"Reading GTF file: {args.gtf}")
    
    line_count = 0
    gene_count = 0
    
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
        
        if feature != args.feature_type:
            continue
        
        attrs = parse_gtf_attributes(attributes)
        
        gene_id = attrs.get('gene_id', '')
        gene_name = attrs.get('gene_name', attrs.get('gene_symbol', gene_id))
        gene_type = attrs.get('gene_biotype', attrs.get('gene_type', 'unknown'))
        
        if gene_type not in allowed_gene_types:
            continue
        
        if args.add_chr_prefix and not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        if args.remove_chr_prefix and chrom.startswith('chr'):
            chrom = chrom[3:]
        
        if '_' in chrom or 'Un' in chrom or 'random' in chrom.lower():
            continue
        
        start = int(start) - 1
        end = int(end)
        
        gene_key = f"{chrom}:{gene_id}"
        
        if gene_key not in genes:
            genes[gene_key] = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'gene_name': gene_name,
                'gene_id': gene_id,
                'strand': strand,
                'gene_type': gene_type
            }
            gene_count += 1
    
    f.close()
    
    print(f"Found {gene_count} genes matching criteria")
    print(f"Writing output file: {args.output}")
    
    with open(args.output, 'w') as out:
        for gene_key in sorted(genes.keys()):
            gene = genes[gene_key]
            name = f"{gene['gene_name']};{gene['gene_id']}"
            line = f"{gene['chrom']}\t{gene['start']}\t{gene['end']}\t{name}\t0\t{gene['strand']}\t{gene['gene_id']}\t{gene['gene_type']}\n"
            out.write(line)
    
    print(f"Done! Output: {args.output}")


if __name__ == '__main__':
    main()
