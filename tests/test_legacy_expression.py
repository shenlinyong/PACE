"""Check the expression-weighted historical model independently of the current kernel."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


SOURCE = Path(__file__).resolve().parents[1] / 'legacy/scripts/sensitivity_analysis.py'
SPEC = importlib.util.spec_from_file_location('legacy_sensitivity', SOURCE)
LEGACY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(LEGACY)


def candidates():
    return pd.DataFrame({
        'gene': ['g1', 'g2', 'g1', 'g2'],
        'contact': [1., 3., 2., 1.],
        'ExpressionWeight': [.2, .8, .2, .8],
    }), np.array([2., 1., 3., 3.])


def test_expression_weight_rescales_each_gene_after_normalization():
    frame, activity = candidates()
    score = LEGACY.compute_pace_score(
        frame, activity, 'contact', 'gene', 'ExpressionWeight')
    np.testing.assert_allclose(score, [.05, .4, .15, .4])
    totals = pd.Series(score).groupby(frame['gene']).sum()
    np.testing.assert_allclose(totals, [.2, .8])


@pytest.mark.parametrize('expression_column', [None, 'absent'])
def test_optional_expression_preserves_unweighted_normalization(expression_column):
    frame, activity = candidates()
    score = LEGACY.compute_pace_score(
        frame, activity, 'contact', 'gene', expression_column)
    np.testing.assert_allclose(score, [.25, .5, .75, .5])


def test_zero_expression_and_zero_support_remain_finite_zero():
    frame, activity = candidates()
    activity[[0, 2]] = 0.
    frame.loc[frame['gene'].eq('g2'), 'ExpressionWeight'] = 0.
    score = LEGACY.compute_pace_score(
        frame, activity, 'contact', 'gene', 'ExpressionWeight')
    np.testing.assert_array_equal(score, np.zeros(4))
