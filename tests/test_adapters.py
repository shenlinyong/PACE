import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'workflow'/'scripts'))
sys.path.insert(0,str(ROOT/'scripts'))
from predictor import PACEPredictor, run_predictions
from calculate_pace_score import PACEScoreCalculator
from multiomics_activity import MultiOmicsActivityCalculator, _normalize
from neighborhoods import NeighborhoodAnalyzer
from pace_core import aggregate_activity
from prepare_tss import prepare


def data():
    e=pd.DataFrame(dict(chr=['1','1'],start=[0,1000],end=[100,1100],name=['e1','e2'],activity=[2.,8.]))
    g=pd.DataFrame(dict(chr=['1','1','1'],gene_id=['g1','g1','g2'],gene_name=['X','X','X'],TSS=[10000,10500,20000],strand=['+','+','-']))
    return e,g


def test_both_prediction_entrypoints_identical(monkeypatch):
    e,g=data()
    c=PACEScoreCalculator({})
    monkeypatch.setattr(c,'add_signals_from_sample',lambda x: None)
    monkeypatch.setattr(c,'calculate_activity',lambda x: (e.activity.to_numpy(),{}))
    a=c.calculate_pace_scores(e,g,{})
    b=PACEPredictor().predict(e,g)
    pd.testing.assert_frame_equal(a,b)


def test_both_activity_entrypoints_identical():
    e,g=data(); e['ATAC']=[3,np.nan];e['H3K27ac']=[0,8]
    c=MultiOmicsActivityCalculator()
    a=c._aggregate_signals({'ATAC':e.ATAC.values,'H3K27ac':e.H3K27ac.values},{},{'ATAC':1,'H3K27ac':1},{})
    n=NeighborhoodAnalyzer(e,g,{'1':100000})
    n.calculate_activity(accessibility_column='ATAC',histone_columns=['H3K27ac'],method='missing_geometric')
    np.testing.assert_allclose(a,n.regions.activity)


def test_output_all_not_overwritten_by_filter(tmp_path):
    e,g=data(); ef=tmp_path/'e.tsv';gf=tmp_path/'g.tsv';of=tmp_path/'out.tsv'
    e.to_csv(ef,sep='\t',index=False);g.to_csv(gf,sep='\t',index=False)
    run_predictions(str(ef),str(gf),str(of),score_threshold=.99)
    assert len(pd.read_csv(of,sep='\t'))==4
    assert len(pd.read_csv(str(of)+'.filtered.tsv',sep='\t'))==0


def test_unknown_expression_not_filtered():
    e,g=data(); x=pd.DataFrame({'gene_id':['g1'],'TPM':[0.]})
    p=PACEPredictor();a=p.predict(e,g);b=p.predict(e,g,expression=x,expression_weight=True)
    np.testing.assert_allclose(a['PACE.Score'],b['PACE.Score'])
    assert set(b.loc[b.TargetGeneEnsemblID=='g2','expression_status'])=={'unknown'}


def test_contact_metadata_reaches_shared_kernel():
    e,g=data();p=PACEPredictor();pairs=p.create_pairs(e,g)
    cols=['chr','start','end','TargetGeneEnsemblID','TargetGeneTSS']
    md=pairs[cols].assign(contact_observed=2.,contact_expected=1.,contact_reliability=.5,contact_source='surrogate')
    r=p.predict(e,g,contact_metadata=md)
    assert set(r.contact_state)=={'surrogate_shrunk'}


def test_normalization_preserves_nan():
    x=_normalize(np.array([0.,np.nan,4.]),'rpkm')
    assert np.isnan(x[1]) and x[0]==0 and x[2]==1e6


def test_tss_gtf_coordinates_and_dedup(tmp_path):
    p=tmp_path/'g.gtf'
    p.write_text('1\tx\tgene\t101\t500\t.\t-\t.\tgene_id "G";\n'
                 '1\tx\ttranscript\t101\t500\t.\t-\t.\tgene_id "G"; transcript_id "a";\n'
                 '1\tx\ttranscript\t101\t500\t.\t-\t.\tgene_id "G"; transcript_id "b";\n'
                 '1\tx\ttranscript\t101\t600\t.\t-\t.\tgene_id "G"; transcript_id "c";\n'
                 '1\tx\tgene\t51\t80\t.\t+\t.\tgene_id "H";\n')
    r=prepare(p)
    assert [x['TSS'] for x in r]==[499,599,50]
    assert [x['tss_weight'] for x in r]==[.5,.5,1.]
    assert r[-1]['tss_source']=='gene_boundary_fallback'


def test_conflicting_duplicate_tss_weights_rejected_in_either_order():
    e, g = data()
    g['tss_weight'] = [.5, .5, 1.]
    extra = g.iloc[[0]].copy()
    extra['tss_weight'] = .2
    conflict = pd.concat([g, extra], ignore_index=True)
    for rows in [conflict, conflict.iloc[::-1]]:
        with pytest.raises(ValueError, match='Conflicting tss_weight'):
            PACEPredictor().predict(e, rows)


def test_duplicate_transcripts_with_consistent_metadata_are_harmless():
    e, g = data()
    extra = g.iloc[[0]].copy()
    extra['transcript_id'] = 'another_transcript'
    a = PACEPredictor().predict(e, g)
    b = PACEPredictor().predict(e, pd.concat([g, extra], ignore_index=True))
    pd.testing.assert_frame_equal(a, b)
