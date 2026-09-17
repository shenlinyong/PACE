"""Small-file integration and failure cases for the public PACE interfaces."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'workflow/scripts'), str(ROOT / 'scripts')]
from neighborhoods import NeighborhoodAnalyzer, run_neighborhoods
from multiomics_activity import quantify_signal_over_regions
from predictor import run_predictions


def cli(script, *args):
    env = dict(os.environ)
    env.pop('PYTHONPATH', None)
    return subprocess.run([sys.executable, str(ROOT / script), *map(str, args)],
                          capture_output=True, text=True, env=env, timeout=60)


def small_files(folder):
    folder.mkdir(parents=True, exist_ok=True)
    e = pd.DataFrame({'chr': ['1', '1', '1'], 'start': [0, 1000, 2000],
                      'end': [100, 1101, 2100], 'name': ['same', 'same', 'zero'],
                      'activity': [20., 1000 / 101, 0.]})
    g = pd.DataFrame({'chr': ['1', '1', '1'], 'gene_id': ['G1', 'G1', 'G2'],
                      'gene_name': ['A', 'A', 'B'], 'TSS': [10000, 10500, 20000],
                      'strand': ['+', '+', '-'], 'tss_weight': [.5, .5, 1.]})
    e.to_csv(folder / 'enhancers.tsv', sep='\t', index=False)
    g.to_csv(folder / 'genes.tsv', sep='\t', index=False)
    e[['chr', 'start', 'end', 'name']].to_csv(folder / 'regions.bed', sep='\t', index=False, header=False)
    (folder / 'reads with spaces.bed').write_text('1\t0\t20\n1\t10\t30\n1\t1000\t1020\n')
    (folder / 'sizes.tsv').write_text('1\t30000\n')
    return e, g


def test_missing_planned_signal_and_quality_columns_are_unknown(tmp_path):
    d = pd.read_csv(ROOT / 'example_quantified/candidates.tsv', sep='\t')
    d = d.drop(columns=['H3K27ac', 'H3K27ac_quality', 'ATAC_quality'])
    inp, out = tmp_path / 'pairs.tsv', tmp_path / 'result.tsv'
    d.to_csv(inp, sep='\t', index=False)
    p = cli('scripts/pace.py', '--pairs', inp, '--activity-config',
            ROOT / 'example_quantified/activity.json', '--output', out)
    assert p.returncode == 0, p.stderr
    result = pd.read_csv(out, sep='\t')
    assert len(result) == 6
    assert result.loc[result.activity.notna(), 'activity_quality'].isna().all()
    assert result.loc[result.activity.isna(), 'evidence_status'].eq('insufficient').all()
    assert not result.evidence_status.eq('sufficient_input_evidence').any()


def test_neighborhood_read_quantification_matches_standalone(tmp_path):
    if not shutil.which('bedtools'):
        pytest.skip('bedtools unavailable')
    e, g = small_files(tmp_path)
    n = NeighborhoodAnalyzer(e, g, {'1': 30000})
    n.quantify_signal_from_bam(str(tmp_path / 'reads with spaces.bed'), 'ATAC')
    expected = np.array([20., 1000 / 101, 0.])
    np.testing.assert_allclose(n.regions.ATAC, expected)
    np.testing.assert_allclose(quantify_signal_over_regions(e, str(tmp_path / 'reads with spaces.bed')), expected)


def test_duplicate_region_names_do_not_merge_counts(tmp_path):
    if not shutil.which('bedtools'):
        pytest.skip('bedtools unavailable')
    e, g = small_files(tmp_path)
    reads = tmp_path / 'reads.bed'
    reads.write_bytes((tmp_path / 'reads with spaces.bed').read_bytes())
    n = NeighborhoodAnalyzer(e, g, {'1': 30000})
    n.quantify_signal_from_bam(str(reads), 'ATAC', normalize=False)
    np.testing.assert_array_equal(n.regions.ATAC, [2., 1., 0.])


def test_numeric_chromosome_bigwig_matches_exact_summary(tmp_path):
    bwmod = pytest.importorskip('pyBigWig')
    e, g = small_files(tmp_path)
    e['chr'] = 1
    path = str(tmp_path / 'signal.bw')
    with bwmod.open(path, 'w') as bw:
        bw.addHeader([('1', 30000)])
        bw.addEntries(['1', '1'], [0, 1000], ends=[100, 1101], values=[2., 4.])
    n = NeighborhoodAnalyzer(e, g, {'1': 30000})
    n.quantify_signal_from_bigwig(path, 'ATAC')
    np.testing.assert_allclose(n.regions.ATAC, [2., 4., np.nan], equal_nan=True)


def test_bed6_candidates_preserve_coordinates(tmp_path):
    if not shutil.which('bedtools'):
        pytest.skip('bedtools unavailable')
    e, _ = small_files(tmp_path)
    reads = tmp_path / 'reads.bed'
    reads.write_bytes((tmp_path / 'reads with spaces.bed').read_bytes())
    e[['chr', 'start', 'end', 'name']].assign(score=0, strand='.').to_csv(
        tmp_path / 'regions.bed', sep='\t', index=False, header=False)
    ef, _ = run_neighborhoods(str(tmp_path / 'regions.bed'), str(tmp_path / 'genes.tsv'),
        str(tmp_path / 'sizes.tsv'), str(tmp_path / 'neighborhoods'), str(reads))
    result = pd.read_csv(ef, sep='\t')
    np.testing.assert_array_equal(result.start, e.start)
    np.testing.assert_array_equal(result.end, e.end)


def test_explicit_missing_expression_file_fails(tmp_path):
    small_files(tmp_path)
    with pytest.raises(FileNotFoundError):
        run_predictions(str(tmp_path / 'enhancers.tsv'), str(tmp_path / 'genes.tsv'),
                        str(tmp_path / 'result.tsv'), expression_file=str(tmp_path / 'absent.tsv'))


def test_explicit_missing_histone_file_fails(tmp_path):
    small_files(tmp_path)
    p = cli('workflow/scripts/pace_neighborhoods.py', '--candidate_regions', tmp_path / 'regions.bed',
            '--genes', tmp_path / 'genes.tsv', '--chrom_sizes', tmp_path / 'sizes.tsv',
            '--output_dir', tmp_path / 'result', '--accessibility_file', tmp_path / 'reads with spaces.bed',
            '--H3K27ac', tmp_path / 'absent.bw')
    assert p.returncode != 0
    assert 'absent.bw' in p.stderr


def test_filtered_output_keeps_parent_directory(tmp_path):
    folder = tmp_path / 'AllPutative_data'
    small_files(folder)
    out = folder / 'EnhancerPredictionsAllPutative.tsv.gz'
    run_predictions(str(folder / 'enhancers.tsv'), str(folder / 'genes.tsv'), str(out))
    assert out.is_file()
    assert (folder / 'EnhancerPredictionsFiltered.tsv').is_file()


def test_standalone_creates_output_directory(tmp_path):
    if not shutil.which('bedtools'):
        pytest.skip('bedtools unavailable')
    small_files(tmp_path)
    (tmp_path / 'samples.tsv').write_text('biosample\tATAC\nsmall\t' + str(tmp_path / 'reads with spaces.bed') + '\n')
    (tmp_path / 'config.json').write_text(json.dumps({'biosamplesTable': str(tmp_path / 'samples.tsv')}))
    out = tmp_path / 'new' / 'result.tsv'
    p = cli('scripts/calculate_pace_score.py', '--config', tmp_path / 'config.json', '--sample', 'small',
            '--candidates', tmp_path / 'regions.bed', '--genes', tmp_path / 'genes.tsv', '--output', out)
    assert p.returncode == 0, p.stderr
    assert len(pd.read_csv(out, sep='\t')) == 6


def test_tss_preparation_creates_output_directory(tmp_path):
    gtf = tmp_path / 'genes.gtf'
    gtf.write_text('1\tsmall\tgene\t101\t500\t.\t-\t.\tgene_id "G";\n')
    out = tmp_path / 'new' / 'tss.tsv'
    p = cli('scripts/prepare_tss.py', '--gtf', gtf, '--output', out)
    assert p.returncode == 0, p.stderr
    assert pd.read_csv(out, sep='\t').TSS.tolist() == [499]


def test_qc_counts_coordinates_and_stable_gene_ids(tmp_path, monkeypatch):
    monkeypatch.setenv('MPLCONFIGDIR', str(tmp_path / 'mpl'))
    from metrics import PACEMetrics
    # Identical display names must not collapse distinct regions or genes.
    d = pd.DataFrame({'chr': ['1'] * 4, 'start': [0, 0, 100, 100],
                      'end': [50, 50, 150, 150], 'name': ['same'] * 4,
                      'TargetGene': ['same'] * 4, 'TargetGeneEnsemblID': ['A', 'B'] * 2,
                      'PACE.Score': [.6, .2, .4, .8]})
    result = PACEMetrics(d).calculate_all_metrics()
    assert result['unique_enhancers'] == 2
    assert result['unique_genes'] == 2
    assert result['mean_enhancers_per_gene'] == 2
    assert result['genes_with_score_ge_0_02'] == 2


def test_filter_writes_explicit_full_and_slim_paths(tmp_path):
    inp = tmp_path / 'predictions.tsv'
    pd.DataFrame({'PACE.Score': [.8, .01], 'raw_support': [8., .1],
                  'TargetGeneEnsemblID': ['G1', 'G2']}).to_csv(inp, sep='\t', index=False)
    slim, full = tmp_path / 'slim.tsv', tmp_path / 'full' / 'table.tsv'
    p = cli('workflow/scripts/pace_filter.py', '--predictions', inp, '--output', slim,
            '--full_output_file', full, '--threshold', .02)
    assert p.returncode == 0, p.stderr
    assert 'raw_support' in pd.read_csv(full, sep='\t')
    assert 'raw_support' not in pd.read_csv(slim, sep='\t')
    assert len(pd.read_csv(full, sep='\t')) == len(pd.read_csv(slim, sep='\t')) == 1
