#!/usr/bin/env python3
"""Annotate gene-matched eQTL support. Lack of support is UNLABELLED, not false."""
import argparse
import pandas as pd


def annotate(predictions, eqtls):
    required=['chr','pos0','gene_id']
    if not set(required).issubset(eqtls):
        raise ValueError('eQTL table requires chr, pos0 (BED0), gene_id; prefilter independently')
    d=predictions.copy()
    if 'TargetGeneEnsemblID' not in d:raise ValueError('Predictions require stable gene IDs')
    index={(str(ch),str(g)):sorted(set(x.pos0)) for (ch,g),x in eqtls.groupby(['chr','gene_id'])}
    import bisect
    counts=[]
    for _,r in d.iterrows():
        positions=index.get((str(r['chr']),str(r.TargetGeneEnsemblID)),[])
        counts.append(bisect.bisect_left(positions,r.end)-bisect.bisect_left(positions,r.start))
    d['eqtl_variant_count']=counts
    d['eqtl_support']=d.eqtl_variant_count.map(lambda n:'gene_matched_association' if n else 'unlabelled')
    return d


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--predictions',required=True)
    p.add_argument('--eqtls',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    annotate(pd.read_csv(a.predictions,sep='\t'),pd.read_csv(a.eqtls,sep='\t')).to_csv(a.output,sep='\t',index=False)
