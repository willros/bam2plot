import pytest
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from bam2plot.main import (
    return_ref_for_plotting,
    coverage_plot,
    plot_gc,
    plot_from_reads,
    plot_depth_histogram,
    plot_depth_histogram_global,
    plot_lorenz_curves,
    plot_insert_size_distribution,
    COLOR_ZERO,
    COLOR_HIGH,
    COLOR_LOW,
)


def test_return_ref_for_plotting_columns(mosdepth_df):
    result = return_ref_for_plotting(mosdepth_df, "refA", thresh=10)
    assert "pos" in result.columns
    assert "rolling" in result.columns
    assert "color" in result.columns
    assert "index" in result.columns


def test_return_ref_for_plotting_color_zero(mosdepth_df):
    result = return_ref_for_plotting(mosdepth_df, "refA", thresh=10)
    # refA: first 100bp have depth=0, so rolling==0 → COLOR_ZERO
    zero_rows = result.filter(pl.col("rolling") == 0)
    colors = zero_rows["color"].unique().to_list()
    assert colors == [COLOR_ZERO]


def test_return_ref_for_plotting_color_high(mosdepth_df):
    result = return_ref_for_plotting(mosdepth_df, "refA", thresh=10)
    # refA: last 100bp have depth=20, rolling>10 → COLOR_HIGH
    high_rows = result.filter(pl.col("rolling") > 10)
    colors = high_rows["color"].unique().to_list()
    assert colors == [COLOR_HIGH]


def test_return_ref_for_plotting_color_low(mosdepth_df):
    result = return_ref_for_plotting(mosdepth_df, "refA", thresh=10)
    # refA: middle 100bp have depth=5, 0 < rolling <= 10 → COLOR_LOW
    low_rows = result.filter(
        (pl.col("rolling") > 0) & (pl.col("rolling") <= 10)
    )
    colors = low_rows["color"].unique().to_list()
    assert colors == [COLOR_LOW]


def test_return_ref_for_plotting_small_ref_no_sampling(mosdepth_df):
    # refA is 300bp (< 10,000), so modulo=1 → no rows dropped
    result = return_ref_for_plotting(mosdepth_df, "refA", thresh=10)
    # Total should be 300 rows (one per base position)
    assert result.height == 300


def test_coverage_plot_returns_figure(plotting_ready_df):
    fig = coverage_plot(
        plotting_ready_df,
        x_col="pos",
        y_col="rolling",
        color_col="color",
        thresh=10,
        rolling_window=1,
        sample_name="test_sample",
    )
    assert isinstance(fig, Figure)


def test_plot_gc_returns_figure():
    df = pl.DataFrame(
        {
            "position": list(range(1, 101)),
            "rolling_gc": [0.5] * 100,
        }
    )
    fig = plot_gc(df, "Test GC Plot")
    assert isinstance(fig, Figure)


def test_return_ref_for_plotting_mid_size_ref():
    """Refs with 10k-99k bases had rolling_window=0 before the fix."""
    n_bases = 50_000
    rows = [("midref", 0, n_bases, 10)]
    raw = pl.DataFrame(rows, schema=["ref", "start", "end", "depth"], orient="row")
    from bam2plot.main import enrich_coverage_df
    df = enrich_coverage_df(raw, thresh=5)
    # Should not raise — rolling_window is clamped to 1
    result = return_ref_for_plotting(df, "midref", thresh=5)
    assert result.height > 0
    assert "color" in result.columns


def test_plot_from_reads_returns_figure(reads_coverage_df):
    fig = plot_from_reads(
        reads_coverage_df,
        sample_name="test_sample",
        ref_name="test_ref",
        window=10,
    )
    assert isinstance(fig, Figure)


def test_plot_depth_histogram_returns_figure(mosdepth_df):
    fig = plot_depth_histogram(mosdepth_df, "refA")
    assert isinstance(fig, Figure)


def test_plot_depth_histogram_global_returns_figure(mosdepth_df):
    fig = plot_depth_histogram_global(mosdepth_df, ["refA", "refB"])
    assert isinstance(fig, Figure)


def test_plot_lorenz_curves_returns_figure(mosdepth_df):
    fig = plot_lorenz_curves(mosdepth_df, ["refA", "refB"])
    assert isinstance(fig, Figure)


def test_plot_insert_size_distribution_returns_figure():
    sizes = np.array([200, 250, 300, 350, 400, 300, 280, 320])
    fig = plot_insert_size_distribution(sizes, "test_sample")
    assert isinstance(fig, Figure)
