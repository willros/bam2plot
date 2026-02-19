import numpy as np
import polars as pl
import pytest
from bam2plot.main import (
    coverage_for_value,
    coverage_for_many_values,
    weighted_median_rle,
    compute_gini_rle,
    compute_lorenz_data_rle,
)


# ---- weighted_median_rle tests ----

def test_weighted_median_uniform():
    """Uniform depth across all bases → median equals the depth."""
    depths = np.array([10, 10, 10])
    weights = np.array([100, 100, 100])
    assert weighted_median_rle(depths, weights) == 10.0


def test_weighted_median_skewed():
    """Most bases at depth 0, small region at high depth."""
    depths = np.array([0, 100])
    weights = np.array([900, 100])
    # 90% of bases at depth 0 → median should be 0
    assert weighted_median_rle(depths, weights) == 0.0


def test_weighted_median_single_interval():
    depths = np.array([42])
    weights = np.array([1000])
    assert weighted_median_rle(depths, weights) == 42.0


def test_weighted_median_empty():
    depths = np.array([], dtype=np.int64)
    weights = np.array([], dtype=np.int64)
    assert weighted_median_rle(depths, weights) == 0.0


# ---- compute_gini_rle tests ----

def test_gini_uniform():
    """All bases at the same depth → Gini should be 0."""
    depths = np.array([10, 10, 10])
    weights = np.array([100, 100, 100])
    assert compute_gini_rle(depths, weights) == pytest.approx(0.0, abs=0.01)


def test_gini_maximally_unequal():
    """Almost all bases at 0, one base at very high depth → Gini near 1."""
    depths = np.array([0, 1000])
    weights = np.array([999, 1])
    gini = compute_gini_rle(depths, weights)
    assert gini > 0.9


def test_gini_all_zero():
    """All depths are zero → Gini should be 0."""
    depths = np.array([0, 0, 0])
    weights = np.array([100, 200, 300])
    assert compute_gini_rle(depths, weights) == pytest.approx(0.0, abs=0.01)


def test_gini_empty():
    depths = np.array([], dtype=np.int64)
    weights = np.array([], dtype=np.int64)
    assert compute_gini_rle(depths, weights) == pytest.approx(0.0, abs=0.01)


# ---- compute_lorenz_data_rle tests ----

def test_lorenz_shape():
    depths = np.array([5, 10, 15])
    weights = np.array([100, 100, 100])
    x, y = compute_lorenz_data_rle(depths, weights)
    # Should have len(depths)+1 points (including origin)
    assert len(x) == 4
    assert len(y) == 4


def test_lorenz_endpoints():
    depths = np.array([5, 10])
    weights = np.array([50, 50])
    x, y = compute_lorenz_data_rle(depths, weights)
    assert x[0] == 0.0
    assert y[0] == 0.0
    assert x[-1] == pytest.approx(1.0)
    assert y[-1] == pytest.approx(1.0)


def test_lorenz_monotonic():
    depths = np.array([1, 5, 10, 20])
    weights = np.array([100, 200, 300, 400])
    x, y = compute_lorenz_data_rle(depths, weights)
    assert all(x[i] <= x[i + 1] for i in range(len(x) - 1))
    assert all(y[i] <= y[i + 1] for i in range(len(y) - 1))


# ---- coverage_for_value tests ----

def test_coverage_for_value_returns_dataframe(single_ref_df):
    result = coverage_for_value(single_ref_df, coverage=5)
    assert isinstance(result, pl.DataFrame)


def test_coverage_for_value_schema(single_ref_df):
    result = coverage_for_value(single_ref_df, coverage=5)
    assert set(result.columns) == {"coverage", "percent", "id"}


def test_coverage_for_value_zero_threshold(single_ref_df):
    result = coverage_for_value(single_ref_df, coverage=0)
    # refA: 100bp at depth 0, 100bp at depth 5, 100bp at depth 20
    # depth > 0: 200 bases out of 300 = 66.67%
    assert result["percent"][0] == pytest.approx(200 / 300 * 100, abs=0.1)


def test_coverage_for_value_high_threshold(single_ref_df):
    result = coverage_for_value(single_ref_df, coverage=100)
    # No bases have depth > 100
    assert result["percent"][0] == pytest.approx(0.0)


def test_coverage_for_value_id_truncation():
    """Ref name longer than 20 chars gets truncated."""
    long_name = "a" * 30
    df = pl.DataFrame(
        {
            "ref": [long_name],
            "start": [0],
            "end": [100],
            "depth": [10],
            "n_bases": [100],
            "total_bases": [100],
        }
    )
    result = coverage_for_value(df, coverage=5)
    assert len(result["id"][0]) == 20


def test_coverage_for_many_values_length(single_ref_df):
    coverage_values = [0, 5, 10, 15, 20]
    result = coverage_for_many_values(single_ref_df, coverage_values)
    assert result.height == len(coverage_values)


def test_coverage_for_many_values_all_ids_same(single_ref_df):
    coverage_values = [0, 5, 10]
    result = coverage_for_many_values(single_ref_df, coverage_values)
    ids = result["id"].unique().to_list()
    assert len(ids) == 1
