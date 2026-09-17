import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'))
from pace_core import aggregate_activity, estimate_contact, score_pairs, ScoreConfig, summarize_resamples


def candidates():
    return pd.DataFrame([
        dict(chr='1', start=e, end=e+100, TargetGeneEnsemblID=g,
             TargetGeneTSS=t, activity=a, contact_prior=c)
        for e,a in [(0,2.), (1000,8.)]
        for g,t,c in [('G1',10000,1.),('G2',20000,2.)]])


def test_unknown_differs_from_measured_zero():
    a = aggregate_activity({'ATAC':[3,3,np.nan,0], 'H3K27ac':[np.nan,0,np.nan,0]})
    assert a.activity.iloc[0] == pytest.approx(3)
    assert a.activity.iloc[1] == pytest.approx(1)
    assert np.isnan(a.activity.iloc[2])
    assert a.activity.iloc[3] == 0
    assert np.isnan(a.activity_quality.iloc[0])


def test_weight_scale_invariance_and_repression_constant():
    s = {'a':[3.,8.], 'b':[0.,3.]}
    a = aggregate_activity(s, {'a':1,'b':2})
    b = aggregate_activity(s, {'a':10,'b':20})
    np.testing.assert_allclose(a.activity,b.activity)
    r = aggregate_activity(s, {'a':1,'b':2}, inhibitory={'meth':[.9,.9]}, inhibition_strength=1)
    np.testing.assert_allclose(r.activity, a.activity*np.exp(-.9))


def test_contact_shrinkage_and_unknown_qc():
    d = estimate_contact([2]*4, [0,8,8,np.nan], [4]*4, [1,0,np.nan,1], ['matched']*4)
    assert list(d.contact) == [0,2,2,2]
    assert d.contact_state.iloc[2] == 'quality_unknown'
    d = estimate_contact([2], [8], [4], [.5], ['surrogate'])
    assert d.contact.iloc[0] == 3
    assert d.contact_state.iloc[0] == 'surrogate_shrunk'


def test_simple_mass_and_gene_id_not_symbol():
    d = candidates().assign(TargetGene='same_symbol')
    r = score_pairs(d)
    assert len(r) == 4
    np.testing.assert_allclose(r.groupby('TargetGeneEnsemblID')['PACE.Score'].sum(),1)
    assert set(r.evidence_status) == {'provisional'}
    assert r.evidence_quality.isna().all()


def test_multiple_tss_duplicates_do_not_inflate():
    d = candidates()
    r = score_pairs(d)
    twice = pd.concat([d,d],ignore_index=True)
    pd.testing.assert_frame_equal(r,score_pairs(twice))
    more = d.copy(); more.TargetGeneTSS += 500
    r2 = score_pairs(pd.concat([d,more],ignore_index=True))
    np.testing.assert_allclose(r.contact_gene,r2.contact_gene)
    np.testing.assert_allclose(r['PACE.Score'],r2['PACE.Score'])


def test_explicit_tss_weights_are_not_renormalized_after_window():
    d = candidates().assign(tss_weight=.25)
    r = score_pairs(d)
    assert sorted(r.contact_gene) == [.25,.25,.5,.5]


def test_missing_activity_is_not_deleted_from_audit():
    d = candidates(); d.loc[d.start==0,'activity']=np.nan
    r = score_pairs(d)
    assert len(r)==4
    assert r.loc[r.start==0,'PACE.Score'].isna().all()
    assert (r.unscored_candidates==1).all()
    assert set(r.loc[r.start==0,'evidence_status']) == {'insufficient'}


def test_residual_mass_prevents_unit_score_when_known():
    d = candidates().query('start==0').assign(unassigned_mass=10.)
    r = score_pairs(d)
    assert (r['PACE.Score']<1).all()
    assert set(r.score_scope)=={'residual_adjusted'}


def test_expression_cannot_change_structural_score():
    d = candidates(); r=score_pairs(d)
    d['Expression']=[1,10000,1,10000]
    np.testing.assert_allclose(r['PACE.Score'],score_pairs(d)['PACE.Score'])


def test_qc_and_score_are_separate():
    d = candidates().assign(activity_quality=1.,tss_quality=1.,catalogue_quality=1.,
                             contact_observed=1.,contact_expected=1.,contact_reliability=1.,contact_source='matched')
    r=score_pairs(d)
    assert set(r.evidence_status)=={'sufficient_input_evidence'}
    d.catalogue_quality=np.nan
    r2=score_pairs(d)
    np.testing.assert_allclose(r['PACE.Score'],r2['PACE.Score'])
    assert set(r2.evidence_status)=={'provisional'}


@pytest.mark.parametrize('column,value',[('activity',-1),('contact_prior',0),('contact_reliability',1.2),('TargetGeneEnsemblID','')])
def test_invalid_input_fails(column,value):
    d=candidates(); d[column]=value
    with pytest.raises(ValueError): score_pairs(d)


def test_conflicting_transcript_records_fail():
    d=candidates(); extra=d.iloc[[0]].copy();extra.contact_prior=4
    with pytest.raises(ValueError):score_pairs(pd.concat([d,extra]))


def test_abc_normalization_limit():
    d=candidates(); r=score_pairs(d,ScoreConfig(competition_power=0))
    expected=d.activity*d.contact_prior
    expected=expected/expected.groupby(d.TargetGeneEnsemblID).transform('sum')
    check=r.merge(d.assign(expected=expected),on=['chr','start','end','TargetGeneEnsemblID'])
    np.testing.assert_allclose(check['PACE.Score'],check.expected)


def test_resampling_denominator_includes_missing_edges():
    r=score_pairs(candidates()); short=r[r.start!=0]
    s=summarize_resamples([r,short])
    assert set(s.loc[s.start==0,'scored_fraction'])=={.5}
    assert set(s.loc[s.start!=0,'scored_fraction'])=={1.}


def test_measured_zero_contact_is_zero_support_not_missing():
    d = candidates().assign(contact_expected=1., contact_reliability=1., contact_source='matched')
    d['contact_observed'] = np.where(d.start == 0, 0., 1.)
    r = score_pairs(d)
    assert r.loc[r.start == 0, 'PACE.Score'].eq(0).all()
    assert r.unscored_candidates.eq(0).all()
    assert r.loc[r.start == 0, 'target_share'].isna().all()
    # With no support anywhere the gene normalization remains undefined.
    d.contact_observed = 0.
    r = score_pairs(d)
    assert r.raw_support.eq(0).all()
    assert r['PACE.Score'].isna().all()
    # Zero contact must not disguise a missing activity measurement.
    d.loc[d.start == 0, 'activity'] = np.nan
    r = score_pairs(d)
    assert r.loc[r.start == 0, 'raw_support'].isna().all()


@pytest.mark.parametrize('column', ['activity_quality', 'tss_quality', 'catalogue_quality'])
def test_quality_conflicts_across_shared_entities_fail(column):
    d = candidates()
    d[column] = 1.
    d.loc[0, column] = np.nan
    with pytest.raises(ValueError, match=column):
        score_pairs(d)
