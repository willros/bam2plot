import pytest
from bam2plot.main import (
    check_range,
    check_threads,
    parse_whitelist,
    refs_with_most_coverage,
)


def test_refs_returns_list(mosdepth_df):
    result = refs_with_most_coverage(mosdepth_df)
    assert isinstance(result, list)


def test_refs_default_returns_all(mosdepth_df):
    result = refs_with_most_coverage(mosdepth_df)
    unique_refs = mosdepth_df.n_unique("ref")
    assert len(result) == unique_refs


def test_refs_top_1(mosdepth_df):
    result = refs_with_most_coverage(mosdepth_df, n=1)
    assert len(result) == 1
    # refB has pct_over_thresh = 40/80 = 0.5; refA has 100/300 = 0.333
    # So refB should come first
    assert result[0] == "refB"


def test_refs_order(mosdepth_df):
    result = refs_with_most_coverage(mosdepth_df, n=2)
    # refB (0.5) should come before refA (0.333)
    assert result == ["refB", "refA"]


def test_refs_over_100_exits(mosdepth_df):
    with pytest.raises(SystemExit):
        refs_with_most_coverage(mosdepth_df, n=101)


def test_check_range_valid():
    assert check_range(1) == 1
    assert check_range(50) == 50
    assert check_range(100) == 100


def test_check_range_invalid():
    with pytest.raises(SystemExit):
        check_range(0)
    with pytest.raises(SystemExit):
        check_range(101)
    with pytest.raises(SystemExit):
        check_range(-5)


def test_check_threads_valid():
    assert check_threads(1) == 1
    assert check_threads(8) == 8


def test_check_threads_invalid():
    with pytest.raises(SystemExit):
        check_threads(0)
    with pytest.raises(SystemExit):
        check_threads(-1)


def test_whitelist_filter_empty_exits(mosdepth_df):
    """Filtering by a whitelist that matches nothing should sys.exit(1)."""
    from bam2plot.main import main_from_bam

    with pytest.raises(SystemExit):
        main_from_bam(
            bam="test/barcode87.bam",
            outpath="/tmp/test_whitelist_empty",
            whitelist=["NONEXISTENT_REF"],
            rolling_window=100,
            threshold=10,
            index=False,
            sort_and_index=False,
            zoom=False,
            cum_plot=False,
            plot_type="png",
            number_of_refs=10,
            threads=None,
        )


def test_parse_whitelist_splits_comma_separated_values():
    assert parse_whitelist("chr1, chr2,chr3") == ["chr1", "chr2", "chr3"]


def test_parse_whitelist_normalizes_mixed_inputs():
    assert parse_whitelist(["chr1, chr2", "chr3", ""]) == ["chr1", "chr2", "chr3"]
