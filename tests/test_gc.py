import polars as pl
from bam2plot.main import add_gc, add_rolling_mean


def test_add_gc_marks_gc_as_1(gc_sequence_df):
    result = add_gc(gc_sequence_df)
    gc_bases = result.filter(pl.col("sequence").str.contains("G|C"))
    assert gc_bases["gc"].to_list() == [1] * gc_bases.height


def test_add_gc_marks_at_as_0(gc_sequence_df):
    result = add_gc(gc_sequence_df)
    at_bases = result.filter(~pl.col("sequence").str.contains("G|C"))
    assert at_bases["gc"].to_list() == [0] * at_bases.height


def test_add_gc_preserves_columns(gc_sequence_df):
    original_cols = set(gc_sequence_df.columns)
    result = add_gc(gc_sequence_df)
    assert original_cols.issubset(set(result.columns))


def test_add_rolling_mean_adds_column(gc_sequence_df):
    df = add_gc(gc_sequence_df)
    result = add_rolling_mean(df, window=3)
    assert "rolling_gc" in result.columns


def test_add_rolling_mean_values(gc_sequence_df):
    # ATGCGCATAT → gc: 0,0,1,1,1,1,0,0,0,0
    df = add_gc(gc_sequence_df)
    result = add_rolling_mean(df, window=3)
    rolling = result["rolling_gc"].to_list()
    # First two values are null (window=3), third onward are rolling means
    import pytest

    assert rolling[0] is None
    assert rolling[1] is None
    # positions 3-onwards: mean of [0,0,1]=0.333, [0,1,1]=0.666, [1,1,1]=1.0, ...
    assert rolling[2] == pytest.approx(1 / 3, abs=1e-6)
    assert rolling[3] == pytest.approx(2 / 3, abs=1e-6)
    assert rolling[4] == pytest.approx(1.0, abs=1e-6)
