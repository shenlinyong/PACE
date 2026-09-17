#!/usr/bin/env python3
"""ABC-only comparator with a fixed labelled universe and explicit missing edges.

Thresholds must be frozen using independent validation, not selected here.
AP_end_to_end integrates only scored predictions against ALL labelled positives;
unscored positives are never artificially retrieved at a final tied threshold.
"""
import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd

KEY=['chr','start','end','TargetGeneEnsemblID']


def evaluate(labels, predictions, score_column, threshold):
    key=KEY + (['sample_id'] if 'sample_id' in labels else [])
    if labels.duplicated(key).any() or predictions.duplicated(key).any():
        raise ValueError('Duplicate edges; resolve labels or TSS aggregation first')
    if not labels.label.dropna().isin([0,1]).all():
        raise ValueError('Labels must be 0, 1, or missing for unlabelled')
    known=labels.dropna(subset=['label'])
    if known.empty or known.label.nunique()!=2:
        raise ValueError('Both experimentally labelled classes are required')
    if score_column not in predictions:
        raise ValueError('Score column not found')
    d=known[key+['label']].merge(predictions[key+[score_column]],on=key,how='left',validate='one_to_one')
    if np.isinf(d[score_column]).any(): raise ValueError('Infinite scores')
    observed=d[d[score_column].notna()].copy()
    positives=int(d.label.sum())
    if observed.empty: ap=0.
    else:
        groups=observed.groupby(score_column).label.agg(['sum','count']).sort_index(ascending=False)
        precision=groups['sum'].cumsum()/groups['count'].cumsum()
        ap=float((precision*groups['sum']/positives).sum())
    selected=d[score_column].ge(threshold)&d[score_column].notna()
    tp=int(d.loc[selected,'label'].sum());n=int(selected.sum())
    result=dict(n_labelled=len(d),n_positive=positives,n_unlabelled=int(labels.label.isna().sum()),
                n_scored=len(observed),scored_fraction=len(observed)/len(d),
                positive_scored_fraction=float(observed.label.sum()/positives),
                AP_end_to_end=ap,threshold=threshold,n_selected=n,true_positives=tp,
                precision=tp/n if n else None,recall=tp/positives)
    return result,d


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--labels',required=True)
    p.add_argument('--abc',required=True)
    p.add_argument('--pace',required=True)
    p.add_argument('--abc-score',default='ABC.Score')
    p.add_argument('--pace-score',default='PACE.Score')
    p.add_argument('--abc-threshold',type=float,required=True)
    p.add_argument('--pace-threshold',type=float,required=True)
    p.add_argument('--output',required=True)
    a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    labels=pd.read_csv(a.labels,sep='\t');results=[]
    for name,path,col,threshold in [('ABC',a.abc,a.abc_score,a.abc_threshold),('PACE',a.pace,a.pace_score,a.pace_threshold)]:
        r,d=evaluate(labels,pd.read_csv(path,sep='\t'),col,threshold)
        r['method']=name;results.append(r)
        d.to_csv(out/f'{name}_labelled_universe.tsv',sep='\t',index=False)
    pd.DataFrame(results).to_csv(out/'ABC_PACE_metrics.tsv',sep='\t',index=False)
    (out/'metrics.json').write_text(json.dumps(results,indent=2,allow_nan=False))


if __name__=='__main__':main()
