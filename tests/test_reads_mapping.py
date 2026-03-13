import polars as pl

from bam2plot.main import _add_alignment_to_depth


def test_add_alignment_to_depth_converts_half_open_interval_to_one_based_positions():
    df = pl.DataFrame({"pos": [1, 2, 3, 4, 5, 6], "depth": [0, 0, 0, 0, 0, 0]})

    result = _add_alignment_to_depth(df, start=2, end=5)

    assert result["depth"].to_list() == [0, 0, 1, 1, 1, 0]


def test_add_alignment_to_depth_handles_reference_start_correctly():
    df = pl.DataFrame({"pos": [1, 2, 3, 4], "depth": [0, 0, 0, 0]})

    result = _add_alignment_to_depth(df, start=0, end=2)

    assert result["depth"].to_list() == [1, 1, 0, 0]
