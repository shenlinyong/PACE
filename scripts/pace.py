#!/usr/bin/env python3
"""Score quantified candidate/TSS tables with the shared PACE kernel."""
import argparse
import json
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'))
from pace_core import ScoreConfig, aggregate_activity, score_pairs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pairs', required=True, help='TSV of enhancer/gene/TSS candidates')
    p.add_argument('--output', required=True)
    p.add_argument('--activity-config', help='JSON with signals: {column: {weight, scale, quality_column}}')
    p.add_argument('--competition-power', type=float, default=1.0)
    p.add_argument('--evidence-threshold', type=float, default=0.5)
    args = p.parse_args()
    d = pd.read_csv(args.pairs, sep='\t')
    if args.activity_config:
        cfg = json.loads(Path(args.activity_config).read_text())
        key = [c for c in ['sample_id', 'chr', 'start', 'end'] if c in d]
        columns = list(cfg['signals'])
        for name, spec in cfg['signals'].items():
            if name not in d: d[name] = float('nan')
            if spec.get('quality_column'):
                quality_column = spec['quality_column']
                if quality_column not in d:
                    d[quality_column] = float('nan')
                columns.append(quality_column)
        columns = list(dict.fromkeys(columns))
        u = d[key + columns].drop_duplicates()
        if u.duplicated(key).any():
            raise ValueError('Conflicting signal values for a shared enhancer')
        s, w, q = {}, {}, {}
        for name, spec in cfg['signals'].items():
            scale = float(spec.get('scale', 1))
            if not 0 < scale < float('inf'):
                raise ValueError('Signal scale must be finite and positive')
            s[name] = u[name].to_numpy() / scale
            w[name] = spec.get('weight', 1)
            if spec.get('quality_column'): q[name] = u[spec['quality_column']].to_numpy()
        a = aggregate_activity(s, w, q)
        u = pd.concat([u[key].reset_index(drop=True), a], axis=1)
        d = d.drop(columns=[c for c in a if c in d]).merge(u, on=key, validate='many_to_one')
    scored = score_pairs(d, ScoreConfig(args.competition_power, args.evidence_threshold))
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    scored.to_csv(args.output, sep='\t', index=False)
    print(f'Wrote {len(scored)} gene-level edges; all scores are uncalibrated.')


if __name__ == '__main__':
    main()
