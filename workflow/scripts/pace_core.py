"""PACE shared, label-free scoring kernel. Scores are not probabilities.

One input row per enhancer/gene/TSS within a single sample and assembly.
All candidate TSS weights must be supplied before distance filtering if a
gene has TSSs outside the candidate window. Missing measurements remain NaN.
See docs/FORMULA.md for the exact definitions and limits.
"""
from dataclasses import dataclass
import numpy as np
import pandas as pd

VERSION = '2026.09.17'


@dataclass(frozen=True)
class ScoreConfig:
    competition_power: float = 1.0
    evidence_threshold: float = 0.5

    def __post_init__(self):
        for name in ('competition_power', 'evidence_threshold'):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f'{name} must be in [0, 1]')


def _nonnegative(x, name, positive=False):
    a = np.asarray(x, dtype=float)
    if a.ndim != 1:
        raise ValueError(f'{name} must be a one-dimensional array')
    if np.isinf(a).any() or (a[np.isfinite(a)] < 0).any():
        raise ValueError(f'{name} must contain finite nonnegative values or NaN')
    if positive and (a[np.isfinite(a)] <= 0).any():
        raise ValueError(f'{name} must be positive when supplied')
    return a


def _quality(x, name):
    a = _nonnegative(x, name)
    if (a[np.isfinite(a)] > 1).any():
        raise ValueError(f'{name} must be in [0, 1] or NaN')
    return a


def aggregate_activity(signals, weights=None, qualities=None,
                       inhibitory=None, inhibitory_weights=None,
                       inhibition_strength=0.0):
    """Missing-aware shifted geometric mean of pre-normalized signals.

    Weights are design priors, not fitted feature importances. Unknown quality
    uses weight 1 for provisional ranking, but remains unknown in QC output.
    Repression is opt-in and requires fraction-scale [0,1] input.
    """
    if not signals:
        raise ValueError('At least one activity modality must be specified')
    names = list(signals)
    s = np.column_stack([_nonnegative(signals[k], k) for k in names])
    w = np.array([(weights or {}).get(k, 1.0) for k in names], dtype=float)
    if not np.isfinite(w).all() or (w < 0).any() or w.sum() <= 0:
        raise ValueError('Activity weights must be finite, nonnegative and sum > 0')
    q = np.column_stack([_quality((qualities or {}).get(k, np.full(len(s), np.nan)),
                                 f'{k} quality') for k in names])
    if len(q) != len(s):
        raise ValueError('Activity quality arrays must have the same length as signals')
    present = np.isfinite(s) & (w > 0)
    effective = present * w * np.where(np.isfinite(q), q, 1.0)
    denom = effective.sum(axis=1)
    logsum = (effective * np.log1p(np.where(present, s, 0))).sum(axis=1)
    activity = np.full(len(s), np.nan)
    np.divide(logsum, denom, out=activity, where=denom > 0)
    activity = np.expm1(activity)
    support = (present * w).sum(axis=1) / w.sum()
    quality_known = ((~present) | np.isfinite(q)).all(axis=1)
    quality = (np.where(present & np.isfinite(q), q, 0) * w).sum(axis=1) / w.sum()
    quality[~quality_known] = np.nan
    if not np.isfinite(inhibition_strength) or inhibition_strength < 0:
        raise ValueError('inhibition_strength must be finite and nonnegative')
    repression = np.zeros(len(s))
    repression_observed = np.zeros(len(s), dtype=int)
    if inhibition_strength and inhibitory:
        rnames = list(inhibitory)
        r = np.column_stack([_quality(inhibitory[k], k) for k in rnames])
        if len(r) != len(s):
            raise ValueError('Inhibitory arrays must have the same length as signals')
        rw = np.array([(inhibitory_weights or {}).get(k, 1.0) for k in rnames])
        if not np.isfinite(rw).all() or (rw < 0).any() or rw.sum() <= 0:
            raise ValueError('Invalid inhibitory weights')
        seen = np.isfinite(r) & (rw > 0)
        rd = (seen * rw).sum(axis=1)
        np.divide((np.where(seen, r, 0) * rw).sum(axis=1), rd,
                  out=repression, where=rd > 0)
        repression_observed = seen.sum(axis=1)
        activity *= np.exp(-inhibition_strength * repression)
    return pd.DataFrame({'activity': activity, 'activity_quality': quality,
                         'activity_observed_fraction': support,
                         'activity_modalities_observed': present.sum(axis=1),
                         'repression_observed': repression_observed,
                         'repression': repression})


def estimate_contact(prior, observed=None, expected=None, reliability=None,
                     source=None):
    """Shrink observed/expected contact ratios toward a distance prior.

    Expected contacts MUST match the observed normalization, sample,
    resolution and genomic distance. Unscaled observations are never mixed
    with the prior. Reliability is local QC, not the observed pair count.
    """
    p = _nonnegative(prior, 'contact_prior', positive=True)
    if not np.isfinite(p).all():
        raise ValueError('Every candidate requires a finite positive contact_prior')
    n = len(p)
    o = _nonnegative(observed if observed is not None else np.full(n, np.nan), 'contact_observed')
    e = _nonnegative(expected if expected is not None else np.full(n, np.nan), 'contact_expected', True)
    q = _quality(reliability if reliability is not None else np.full(n, np.nan), 'contact_reliability')
    src = np.asarray(source if source is not None else ['distance_prior'] * n, dtype=str)
    if not (len(o) == len(e) == len(q) == len(src) == n):
        raise ValueError('Contact arrays must have equal lengths')
    if not np.isin(src, ['matched', 'surrogate', 'distance_prior', 'unknown']).all():
        raise ValueError('Unknown contact source')
    usable = np.isfinite(o) & np.isfinite(e) & np.isfinite(q) & np.isin(src, ['matched', 'surrogate'])
    lam = np.where(usable, q, 0.0)
    ratio = np.ones(n)
    np.divide(o, e, out=ratio, where=usable)
    contact = p * ((1 - lam) + lam * ratio)
    if not np.isfinite(contact).all():
        raise ValueError('Contact overflow; check normalization and expected contacts')
    state = np.full(n, 'prior_only', dtype=object)
    state[np.isfinite(o) & ~np.isfinite(e)] = 'unscaled_observation'
    state[np.isfinite(o) & np.isfinite(e) & ~np.isfinite(q)] = 'quality_unknown'
    state[usable] = np.char.add(src[usable], '_shrunk')
    return pd.DataFrame({'contact': contact, 'contact_weight': lam,
                         'contact_state': state})


def score_pairs(pairs, config=None):
    """Collapse weighted TSSs, allocate enhancer activity, then score genes.

    Canonical columns: chr, start, end, TargetGeneEnsemblID, TargetGeneTSS,
    activity, contact_prior. Optional sample_id scopes normalization. Optional
    tss_weight is already normalized over ALL distinct TSSs of each gene.
    unassigned_mass is a pre-estimated gene-level nonnegative residual in the
    same units as raw_support; NaN means unknown, not estimated zero.
    """
    cfg = config or ScoreConfig()
    d = pairs.copy().reset_index(drop=True)
    required = ['chr', 'start', 'end', 'TargetGeneEnsemblID', 'TargetGeneTSS', 'activity', 'contact_prior']
    missing = set(required) - set(d.columns)
    if missing:
        raise ValueError(f'Missing columns: {sorted(missing)}')
    if d.empty:
        for c in ['PACE.Score', 'raw_support', 'evidence_quality', 'evidence_status']:
            d[c] = pd.Series(dtype=str if c == 'evidence_status' else float)
        return d
    if 'sample_id' not in d:
        d['sample_id'] = 'single_sample'
    for c in ['sample_id', 'chr', 'TargetGeneEnsemblID']:
        if d[c].isna().any() or d[c].astype(str).str.strip().eq('').any():
            raise ValueError(f'{c} cannot be missing; use stable IDs, not symbols')
    for c in ['start', 'end', 'TargetGeneTSS']:
        a = _nonnegative(d[c], c)
        if not np.isfinite(a).all() or (a != np.floor(a)).any():
            raise ValueError(f'{c} must be finite integer coordinates')
    if (d.end <= d.start).any():
        raise ValueError('Enhancer intervals must be BED0 half-open with end > start')
    ek = ['sample_id', 'chr', 'start', 'end']
    gk = ['sample_id', 'chr', 'TargetGeneEnsemblID']
    pk = ek + ['TargetGeneEnsemblID']
    tk = pk + ['TargetGeneTSS']
    # Exact duplicate transcript records are harmless; conflicts are errors.
    d = d.drop_duplicates().reset_index(drop=True)
    if d.duplicated(tk).any():
        raise ValueError('Conflicting or redundant enhancer/gene/TSS records; deduplicate TSSs first')
    d['activity'] = _nonnegative(d.activity, 'activity')
    if (d.groupby(ek).activity.nunique(dropna=False) > 1).any():
        raise ValueError('Activity must be enhancer-specific, not target-specific')
    for c in ['activity_quality', 'tss_quality', 'catalogue_quality']:
        if c not in d:
            d[c] = np.nan
        d[c] = _quality(d[c], c)
    for c, scope in [('activity_quality', ek),
                     ('tss_quality', gk + ['TargetGeneTSS']),
                     ('catalogue_quality', gk)]:
        if (d.groupby(scope)[c].nunique(dropna=False) > 1).any():
            raise ValueError(f'{c} must be consistent within its biological entity')
    if 'unassigned_mass' not in d:
        d['unassigned_mass'] = np.nan
    d['unassigned_mass'] = _nonnegative(d.unassigned_mass, 'unassigned_mass')
    if (d.groupby(gk).unassigned_mass.nunique(dropna=False) > 1).any():
        raise ValueError('unassigned_mass must be constant within a gene')
    if 'tss_weight' not in d:
        ts = d[gk + ['TargetGeneTSS']].drop_duplicates()
        ts['tss_weight'] = 1 / ts.groupby(gk).TargetGeneTSS.transform('size')
        d = d.merge(ts, on=gk + ['TargetGeneTSS'], how='left', validate='many_to_one')
    d['tss_weight'] = _quality(d.tss_weight, 'tss_weight')
    if not np.isfinite(d.tss_weight).all():
        raise ValueError('tss_weight cannot be missing')
    ts = d[gk + ['TargetGeneTSS', 'tss_weight']].drop_duplicates()
    if ts.duplicated(gk + ['TargetGeneTSS']).any():
        raise ValueError('TSS weights must be gene-specific, not enhancer-specific')
    total_weights = ts.groupby(gk).tss_weight.sum()
    if ((total_weights <= 0) | (total_weights > 1 + 1e-9)).any():
        raise ValueError('Unique TSS weights must sum to (0,1] for each gene')
    con = estimate_contact(d.contact_prior, d.get('contact_observed'),
                           d.get('contact_expected'), d.get('contact_reliability'),
                           d.get('contact_source'))
    for c in con:
        d[c] = con[c].to_numpy()
    d['_weighted_contact'] = d.contact * d.tss_weight
    d['_weighted_q'] = d.contact_weight * d.tss_weight
    single_tss = not d.duplicated(pk).any()
    if single_tss:
        # Avoid millions of one-row Python group callbacks in atlas inference.
        # This is algebraically identical to the general aggregation below.
        out = d.drop(columns=['_weighted_contact', '_weighted_q']).copy()
        out['contact_gene'] = d['_weighted_contact'].to_numpy()
        out['contact_quality'] = d['_weighted_q'].to_numpy()
        out['n_tss'] = 1
    else:
        group = d.groupby(pk, sort=False, dropna=False)
        # Preserve explicitly unknown values; pandas first() would skip NaN.
        out = group.nth(0).copy().reset_index(drop=True)
        agg = group.agg(contact_gene=('_weighted_contact', 'sum'),
                        contact_quality=('_weighted_q', 'sum'),
                        n_tss=('TargetGeneTSS', 'size')).reset_index()
        out = out.drop(columns=['_weighted_contact', '_weighted_q']).merge(agg, on=pk, validate='one_to_one')
        if d.contact_state.nunique() == 1:
            out['contact_state'] = d.contact_state.iloc[0]
        else:
            states = group.contact_state.agg(lambda x: '|'.join(sorted(set(x)))).reset_index()
            out = out.drop(columns=['contact_state']).merge(states, on=pk, validate='one_to_one')
        # Built-in reductions preserve unknown QC using an explicit count check.
        for c in ['tss_quality', 'activity_quality', 'catalogue_quality']:
            qc = group[c].agg(['min', 'count', 'size'])
            qc[c] = qc['min'].where(qc['count'].eq(qc['size']))
            out = out.drop(columns=[c]).merge(qc[[c]].reset_index(), on=pk, validate='one_to_one')
    target_total = out.groupby(ek).contact_gene.transform('sum')
    if not np.isfinite(target_total).all():
        raise ValueError('Contact-allocation overflow; reduce input contact scales')
    out['target_share'] = out.contact_gene.div(target_total.where(target_total > 0))
    factor = out.target_share ** cfg.competition_power if cfg.competition_power else 1.0
    out['raw_support'] = out.activity * out.contact_gene * factor
    # Allocation is undefined when all candidate contacts are measured zero,
    # but finite activity times zero contact has known zero support. Preserve
    # missing activity as unknown and leave target_share itself undefined.
    known_zero = np.isfinite(out.activity) & out.contact_gene.eq(0)
    out.loc[known_zero, 'raw_support'] = 0.0
    if np.isinf(out.raw_support).any():
        raise ValueError('Support overflow; check input signal scales')
    gene_support = out.groupby(gk).raw_support
    mass = gene_support.transform('sum', min_count=1)
    if ((gene_support.transform('count') > 0) & ~np.isfinite(mass)).any():
        raise ValueError('Gene-normalization overflow; reduce input signal scales')
    out['unscored_candidates'] = gene_support.transform('size') - gene_support.transform('count')
    out['observed_mass'] = mass
    denom = mass + out.unassigned_mass.fillna(0)
    if np.isinf(denom).any():
        raise ValueError('Residual-normalization overflow; reduce input support scales')
    out['PACE.Score'] = out.raw_support.div(denom.where(denom > 0))
    out['score_scope'] = np.where(out.unassigned_mass.notna(), 'residual_adjusted', 'observed_candidates_only')
    qc_cols = ['activity_quality', 'contact_quality', 'tss_quality', 'catalogue_quality']
    out['evidence_quality'] = out[qc_cols].min(axis=1, skipna=False)
    good = out[qc_cols].ge(cfg.evidence_threshold).all(axis=1) & (out.unscored_candidates == 0)
    out['evidence_status'] = np.where(good, 'sufficient_input_evidence', 'provisional')
    insufficient = out['PACE.Score'].isna() | out.activity_quality.eq(0) | out.tss_quality.eq(0)
    out.loc[insufficient, 'evidence_status'] = 'insufficient'
    reasons = pd.Series('', index=out.index)
    for c in qc_cols:
        reasons += np.where(out[c].isna(), c + '_unknown;',
                            np.where(out[c] < cfg.evidence_threshold, c + '_low;', ''))
    for mask, tag in [(out.unscored_candidates.gt(0), 'unscored_candidates'),
                      (out.unassigned_mass.isna(), 'unassigned_mass_unknown'),
                      (out['PACE.Score'].isna(), 'unscorable'),
                      (out.contact_state.str.contains('surrogate'), 'surrogate_contact')]:
        reasons += np.where(mask, tag + ';', '')
    out['evidence_reasons'] = reasons.str.rstrip(';').replace('', 'none')
    out['model_version'] = VERSION
    # Keep all alternatives explicitly; the representative TSS is not a claim
    # that this particular promoter mediates the inferred gene-level effect.
    if single_tss:
        out['TargetGeneTSSs'] = out.TargetGeneTSS.astype(str)
    else:
        alternatives = group.TargetGeneTSS.agg(lambda x: ','.join(map(str, sorted(set(x))))).reset_index(name='TargetGeneTSSs')
        out = out.merge(alternatives, on=pk, validate='one_to_one')
    return out.sort_values('PACE.Score', ascending=False, na_position='last').reset_index(drop=True)


def summarize_resamples(frames):
    """Empirical rerun stability, NOT a posterior/causal confidence interval.

    Each frame must come from an actual independent resample and full rescore.
    Edges absent/unscorable in a run are tracked as missing, never dropped from
    the detection-frequency denominator.
    """
    if len(frames) < 2:
        raise ValueError('At least two independently rescored replicates required')
    key = ['sample_id', 'chr', 'start', 'end', 'TargetGeneEnsemblID']
    all_rows = []
    for i, f in enumerate(frames):
        if f.duplicated(key).any():
            raise ValueError('Each resample must contain unique gene-level edges')
        all_rows.append(f[key + ['PACE.Score']].assign(resample=i))
    d = pd.concat(all_rows, ignore_index=True)
    g = d.groupby(key)['PACE.Score']
    result = g.agg(score_median='median', n_scored='count').reset_index()
    result['scored_fraction'] = result.n_scored / len(frames)
    for level, name in [(0.05, 'score_p05'), (0.95, 'score_p95')]:
        result = result.merge(g.quantile(level).reset_index(name=name), on=key)
    return result
