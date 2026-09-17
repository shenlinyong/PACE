#!/usr/bin/env python3
"""Run a small, independently checked PACE command-line pipeline.

Author: 申林用 (Linyong Shen), Northwest A&F University.
The input files are software fixtures, not biological measurements.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KEY = ['start', 'end', 'TargetGeneEnsemblID']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    out = Path(args.output_dir).resolve()
    data = out / 'inputs'
    logs = out / 'logs'
    data.mkdir(parents=True, exist_ok=True)
    logs.mkdir(exist_ok=True)
    env = dict(os.environ)
    env.pop('PYTHONPATH', None)
    env['MPLBACKEND'] = 'Agg'
    env['MPLCONFIGDIR'] = str(out / 'matplotlib')
    commands = []

    def run(name, script, *arguments):
        cmd = [sys.executable, str(ROOT / script), *map(str, arguments)]
        result = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=120)
        (logs / f'{name}.log').write_text(result.stdout + result.stderr)
        commands.append({'step': name, 'command': cmd, 'exit_code': result.returncode})
        if result.returncode:
            raise RuntimeError(f'{name} failed; see {logs / (name + ".log")}')
        print(f'{name}: passed', flush=True)

    regions = data / 'regions.bed'
    regions.write_text('1\t0\t100\tsame\t0\t.\n1\t1000\t1101\tsame\t0\t.\n1\t2000\t2100\tzero\t0\t.\n')
    reads = data / 'reads with spaces.bed'
    reads.write_text('1\t0\t20\n1\t10\t30\n1\t1000\t1020\n')
    sizes = data / 'sizes.tsv'
    sizes.write_text('1\t30000\n')
    gtf = data / 'genes.gtf'
    gtf.write_text('1\ttest\ttranscript\t10001\t11000\t.\t+\t.\tgene_id "G1"; gene_name "A"; transcript_id "t1";\n'
                   '1\ttest\ttranscript\t10501\t11000\t.\t+\t.\tgene_id "G1"; gene_name "A"; transcript_id "t2";\n'
                   '1\ttest\ttranscript\t19001\t20001\t.\t-\t.\tgene_id "G2"; gene_name "B"; transcript_id "t3";\n')
    genes = out / 'tss' / 'genes.tsv'
    run('prepare_tss', 'scripts/prepare_tss.py', '--gtf', gtf, '--output', genes)
    assert pd.read_csv(genes, sep='\t').TSS.tolist() == [10000, 10500, 20000]
    samples = data / 'samples.tsv'
    samples.write_text(f'biosample\tATAC\nsmall\t{reads}\n')
    config = data / 'config.json'
    config.write_text(json.dumps({'biosamplesTable': str(samples), 'activity_method': 'missing_geometric'}))
    standalone = out / 'standalone' / 'predictions.tsv'
    run('standalone', 'scripts/calculate_pace_score.py', '--config', config, '--sample', 'small',
        '--candidates', regions, '--genes', genes, '--output', standalone)
    neighborhood = out / 'neighborhoods'
    run('neighborhoods', 'workflow/scripts/pace_neighborhoods.py', '--candidate_regions', regions,
        '--genes', genes, '--chrom_sizes', sizes, '--output_dir', neighborhood,
        '--accessibility_file', reads)
    prediction = out / 'AllPutative_data' / 'EnhancerPredictionsAllPutative.tsv.gz'
    run('prediction', 'workflow/scripts/pace_predict.py', '--enhancers', neighborhood / 'EnhancerList.txt',
        '--genes', neighborhood / 'GeneList.txt', '--output', prediction)
    a = pd.read_csv(standalone, sep='\t').sort_values(KEY).reset_index(drop=True)
    b = pd.read_csv(prediction, sep='\t').sort_values(KEY).reset_index(drop=True)
    assert len(a) == len(b) == 6
    np.testing.assert_allclose(a['PACE.Score'], b['PACE.Score'], rtol=1e-12, atol=1e-14)
    # Reference calculation deliberately does not import the PACE kernel.
    midpoint = np.array([50., 1050.5, 2050.])
    activity = np.array([20., 1000 / 101, 0.])
    prior = 5.9594510043736655 / (np.abs(midpoint[:, None] - np.array([10000, 10500, 20000])) + 5000.) ** 1.024238616787792
    contact = np.column_stack([prior[:, :2].mean(axis=1), prior[:, 2]])
    raw = activity[:, None] * contact * (contact / contact.sum(axis=1, keepdims=True))
    expected = raw / raw.sum(axis=0, keepdims=True)
    error = float(np.max(np.abs(b['PACE.Score'].to_numpy() - expected.ravel())))
    np.testing.assert_allclose(b['PACE.Score'], expected.ravel(), rtol=1e-12, atol=1e-14)
    assert b.loc[b.start.eq(2000), 'PACE.Score'].eq(0).all()
    assert b.activity_quality.isna().all()
    assert b.evidence_status.eq('provisional').all()

    # Measured contacts: zero is observed; omitted BEDPE entries remain missing.
    contacts = data / 'contacts.bedpe'
    contacts.write_text('1\t0\t100\t1\t10000\t10100\t0\n'
                        '1\t0\t100\t1\t10500\t10600\t4\n'
                        '1\t1000\t1100\t1\t10000\t10100\t2\n')
    metadata = pd.DataFrame([dict(chr='1', start=start, end=end, TargetGeneEnsemblID=gene,
                                 TargetGeneTSS=tss, contact_expected=1., contact_reliability=.5,
                                 contact_source='matched')
                             for start, end in [(0, 100), (1000, 1101), (2000, 2100)]
                             for gene, tss in [('G1', 10000), ('G1', 10500), ('G2', 20000)]])
    md = data / 'contacts_qc.tsv'
    metadata.to_csv(md, sep='\t', index=False)
    measured = out / 'measured' / 'predictions.tsv'
    run('measured_contact', 'workflow/scripts/pace_predict.py', '--enhancers', neighborhood / 'EnhancerList.txt',
        '--genes', genes, '--output', measured, '--hic_file', contacts, '--hic_type', 'bedpe',
        '--hic_resolution', 100, '--contact_metadata', md)
    measured_df = pd.read_csv(measured, sep='\t').sort_values(KEY).reset_index(drop=True)
    modified = prior.copy()
    modified[0, 0] *= .5
    modified[0, 1] *= 2.5
    modified[1, 0] *= 1.5
    mc = np.column_stack([modified[:, :2].mean(axis=1), modified[:, 2]])
    mr = activity[:, None] * mc * mc / mc.sum(axis=1, keepdims=True)
    np.testing.assert_allclose(measured_df['PACE.Score'], (mr / mr.sum(axis=0)).ravel(), rtol=1e-12)
    assert measured_df.iloc[0].contact_observed == 0
    assert measured_df.loc[measured_df.TargetGeneEnsemblID.eq('G2'), 'contact_observed'].isna().all()

    quant = out / 'quantified' / 'predictions.tsv'
    run('quantified_example', 'scripts/pace.py', '--pairs', ROOT / 'example_quantified/candidates.tsv',
        '--activity-config', ROOT / 'example_quantified/activity.json', '--output', quant)
    example = pd.read_csv(quant, sep='\t')
    assert len(example) == 6 and example['PACE.Score'].isna().sum() == 2
    filtered = out / 'filtered' / 'predictions.tsv'
    run('filter', 'workflow/scripts/pace_filter.py', '--predictions', prediction,
        '--output', filtered, '--threshold', .02)
    assert len(pd.read_csv(filtered, sep='\t')) == 4
    assert len(pd.read_csv(filtered.with_name('predictions_Full.tsv'), sep='\t')) == 4
    run('metrics', 'workflow/scripts/pace_metrics.py', '--predictions', prediction,
        '--output_dir', out / 'metrics', '--sample_name', 'small')
    qc = pd.read_csv(out / 'metrics/QCSummary_small.tsv', sep='\t').iloc[0]
    assert qc.unique_enhancers == 3 and qc.unique_genes == 2
    report = {'status': 'PASS', 'author': '申林用 (Linyong Shen)',
              'fixture': '3 enhancers, 2 genes, 3 distinct TSSs; software test data only',
              'gene_level_edges': 6, 'filtered_edges': 4, 'max_score_error': error,
              'independent_formula_check': True, 'standalone_and_workflow_agree': True,
              'measured_zero_and_missing_contact_checked': True,
              'commands': commands,
              'limits': ['No FASTQ alignment or peak calling in this small-data run.',
                         'Snakemake scheduling and optional binary Hi-C formats are not exercised.']}
    (out / 'smoke_report.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    files = [p for p in out.rglob('*') if p.is_file() and p.name != 'SHA256SUMS.txt' and 'matplotlib' not in p.parts]
    (out / 'SHA256SUMS.txt').write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}\n' for p in sorted(files)))
    print(json.dumps({k: report[k] for k in ['status', 'gene_level_edges', 'filtered_edges', 'max_score_error']}, indent=2))


if __name__ == '__main__':
    main()
