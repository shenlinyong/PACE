#!/usr/bin/env python3
"""GTF -> deduplicated, strand-aware BED0 TSS catalogue for PACE.

All gene biotypes are kept; gene-boundary fallback is explicitly labelled.
This does not discover missing genes or validate annotated promoters.
"""
import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


def prepare(path):
    genes={}; transcripts=defaultdict(dict)
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            cols=line.rstrip('\n').split('\t')
            if len(cols)!=9 or cols[2] not in ('gene','transcript'):continue
            ch,_,feature,start,end,_,strand,_,attrs=cols
            a=dict(re.findall(r'(\S+)\s+"([^"]*)"',attrs))
            if 'gene_id' not in a or strand not in ('+','-'):continue
            tss=int(start)-1 if strand=='+' else int(end)-1
            record=dict(chr=ch,gene_id=a['gene_id'],gene_name=a.get('gene_name',a['gene_id']),
                        TSS=tss,strand=strand,tss_source='annotated_transcript' if feature=='transcript' else 'gene_boundary_fallback',
                        tss_quality='')
            key=(ch,a['gene_id'])
            if feature=='gene':genes[key]=record
            else:transcripts[key][tss]=record
    records=[]
    for key in sorted(set(genes)|set(transcripts)):
        sites=list(transcripts[key].values()) or [genes[key]]
        for r in sorted(sites,key=lambda x:x['TSS']):
            r['tss_weight']=1/len(sites);records.append(r)
    return records


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--gtf',required=True)
    p.add_argument('--output',required=True);a=p.parse_args();records=prepare(a.gtf)
    if not records:raise ValueError('No annotated genes or transcripts found')
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    with open(a.output,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
