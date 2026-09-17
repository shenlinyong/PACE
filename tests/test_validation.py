import sys
from pathlib import Path
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from benchmark_compare import evaluate
from eqtl_pair_support import annotate


def labels():
    return pd.DataFrame(dict(chr=['1']*3,start=[0,100,200],end=[10,110,210],TargetGeneEnsemblID=['G']*3,label=[1,0,1]))


def test_missing_positive_penalizes_end_to_end_ap():
    d=labels();pred=d.iloc[:2].assign(score=[.9,.1]);m,u=evaluate(d,pred,'score',.5)
    assert m['AP_end_to_end']==pytest.approx(.5)
    assert m['recall']==.5 and m['precision']==1
    assert len(u)==3 and m['n_scored']==2


def test_changing_target_destroys_eqtl_support():
    d=labels();q=pd.DataFrame(dict(chr=['1'],pos0=[5],gene_id=['G']))
    assert annotate(d,q).eqtl_variant_count.sum()==1
    d.TargetGeneEnsemblID='WRONG'
    assert annotate(d,q).eqtl_variant_count.sum()==0
    assert set(annotate(d,q).eqtl_support)=={'unlabelled'}
