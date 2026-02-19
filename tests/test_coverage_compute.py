import numpy as np
import pysam
import polars as pl
import pytest

from bam2plot.main import bam_to_raw_df, enrich_coverage_df


@pytest.fixture
def synthetic_bam(tmp_path):
    """Create a minimal BAM with two references and known read positions.

    refA (200 bp): two reads at [10,60) and [50,100) → overlap at [50,60)
    refB (100 bp): one read at [0,40)
    """
    bam_path = str(tmp_path / "test.bam")

    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.0", "SO": "coordinate"},
        "SQ": [
            {"SN": "refA", "LN": 200},
            {"SN": "refB", "LN": 100},
        ],
    })

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        # refA read 1: [10, 60)
        a = pysam.AlignedSegment(outf.header)
        a.query_name = "read1"
        a.query_sequence = "A" * 50
        a.flag = 0
        a.reference_id = 0  # refA
        a.reference_start = 10
        a.mapping_quality = 60
        a.cigar = [(0, 50)]  # 50M
        a.query_qualities = pysam.qualitystring_to_array("I" * 50)
        outf.write(a)

        # refA read 2: [50, 100)
        a = pysam.AlignedSegment(outf.header)
        a.query_name = "read2"
        a.query_sequence = "A" * 50
        a.flag = 0
        a.reference_id = 0  # refA
        a.reference_start = 50
        a.mapping_quality = 60
        a.cigar = [(0, 50)]  # 50M
        a.query_qualities = pysam.qualitystring_to_array("I" * 50)
        outf.write(a)

        # refB read 1: [0, 40)
        a = pysam.AlignedSegment(outf.header)
        a.query_name = "read3"
        a.query_sequence = "A" * 40
        a.flag = 0
        a.reference_id = 1  # refB
        a.reference_start = 0
        a.mapping_quality = 60
        a.cigar = [(0, 40)]  # 40M
        a.query_qualities = pysam.qualitystring_to_array("I" * 40)
        outf.write(a)

    return bam_path


# ---- bam_to_raw_df tests ----

def test_bam_to_raw_df_schema(synthetic_bam):
    df = bam_to_raw_df(synthetic_bam)
    assert set(df.columns) == {"ref", "start", "end", "depth"}


def test_bam_to_raw_df_both_refs_present(synthetic_bam):
    df = bam_to_raw_df(synthetic_bam)
    refs = df["ref"].unique().sort().to_list()
    assert "refA" in refs
    assert "refB" in refs


def test_bam_to_raw_df_intervals_cover_full_ref(synthetic_bam):
    """Intervals for each ref should cover from 0 to ref_length with no gaps."""
    df = bam_to_raw_df(synthetic_bam)

    for ref_name, ref_len in [("refA", 200), ("refB", 100)]:
        ref_df = df.filter(pl.col("ref") == ref_name).sort("start")
        assert ref_df["start"][0] == 0, f"{ref_name} should start at 0"
        assert ref_df["end"][-1] == ref_len, f"{ref_name} should end at {ref_len}"
        # Check no gaps: each row's start == previous row's end
        starts = ref_df["start"].to_list()
        ends = ref_df["end"].to_list()
        for i in range(1, len(starts)):
            assert starts[i] == ends[i - 1], f"Gap in {ref_name} at position {ends[i-1]}"


def test_bam_to_raw_df_depth_values(synthetic_bam):
    """refA: 0-depth at [0,10), depth 1 at [10,50), depth 2 at [50,60), depth 1 at [60,100), 0-depth at [100,200)."""
    df = bam_to_raw_df(synthetic_bam)
    ref_a = df.filter(pl.col("ref") == "refA").sort("start")

    # Build a lookup: position -> depth
    depth_at = {}
    for row in ref_a.iter_rows(named=True):
        for pos in range(row["start"], row["end"]):
            depth_at[pos] = row["depth"]

    # Check expected depths
    assert depth_at[0] == 0
    assert depth_at[9] == 0
    assert depth_at[10] == 1
    assert depth_at[49] == 1
    assert depth_at[50] == 2  # overlap region
    assert depth_at[59] == 2
    assert depth_at[60] == 1
    assert depth_at[99] == 1
    assert depth_at[100] == 0
    assert depth_at[199] == 0


def test_bam_to_raw_df_zero_depth_gaps(synthetic_bam):
    """Regions with no reads should have depth 0."""
    df = bam_to_raw_df(synthetic_bam)
    ref_b = df.filter(pl.col("ref") == "refB").sort("start")

    # refB: [0,40) depth 1, [40,100) depth 0
    depth_at = {}
    for row in ref_b.iter_rows(named=True):
        for pos in range(row["start"], row["end"]):
            depth_at[pos] = row["depth"]

    assert depth_at[0] == 1
    assert depth_at[39] == 1
    assert depth_at[40] == 0
    assert depth_at[99] == 0


# ---- enrich_coverage_df tests ----

@pytest.fixture
def raw_df():
    """Simple raw DataFrame for testing enrich_coverage_df."""
    return pl.DataFrame({
        "ref": ["refA", "refA", "refA", "refB", "refB"],
        "start": [0, 100, 200, 0, 40],
        "end": [100, 200, 300, 40, 80],
        "depth": [0, 5, 20, 15, 0],
    })


def test_enrich_coverage_df_schema(raw_df):
    df = enrich_coverage_df(raw_df, thresh=10)
    expected_cols = {
        "ref", "start", "end", "depth", "n_bases", "total_bases",
        "mean_coverage", "mean_coverage_total", "pct_over_zero",
        "pct_over_thresh", "pct_total_over_zero", "pct_total_over_thresh",
        "median_coverage", "median_coverage_total", "gini_coefficient",
    }
    assert set(df.columns) == expected_cols


def test_enrich_n_bases(raw_df):
    df = enrich_coverage_df(raw_df, thresh=10)
    n_bases = df["n_bases"].to_list()
    assert n_bases == [100, 100, 100, 40, 40]


def test_enrich_total_bases_per_ref(raw_df):
    df = enrich_coverage_df(raw_df, thresh=10)
    # refA rows should all have total_bases=300, refB should all have 80
    ref_a = df.filter(pl.col("ref") == "refA")
    ref_b = df.filter(pl.col("ref") == "refB")
    assert ref_a["total_bases"].unique().to_list() == [300]
    assert ref_b["total_bases"].unique().to_list() == [80]


def test_enrich_pct_ranges(raw_df):
    df = enrich_coverage_df(raw_df, thresh=10)
    assert (df["pct_over_zero"] >= 0).all()
    assert (df["pct_over_zero"] <= 1).all()
    assert (df["pct_over_thresh"] >= 0).all()
    assert (df["pct_over_thresh"] <= 1).all()


# ---- indexed BAM tests (same sweep-line, but with an indexed BAM) ----

@pytest.fixture
def synthetic_indexed_bam(synthetic_bam):
    """Sort and index the synthetic BAM."""
    sorted_path = synthetic_bam.replace(".bam", ".sorted.bam")
    pysam.sort("-o", sorted_path, synthetic_bam)
    pysam.index(sorted_path)
    return sorted_path


def test_bam_to_raw_df_indexed_schema(synthetic_indexed_bam):
    df = bam_to_raw_df(synthetic_indexed_bam)
    assert set(df.columns) == {"ref", "start", "end", "depth"}


def test_bam_to_raw_df_indexed_depth_values(synthetic_indexed_bam):
    """Same depth expectations as the unindexed test."""
    df = bam_to_raw_df(synthetic_indexed_bam)
    ref_a = df.filter(pl.col("ref") == "refA").sort("start")

    depth_at = {}
    for row in ref_a.iter_rows(named=True):
        for pos in range(row["start"], row["end"]):
            depth_at[pos] = row["depth"]

    assert depth_at[0] == 0
    assert depth_at[9] == 0
    assert depth_at[10] == 1
    assert depth_at[49] == 1
    assert depth_at[50] == 2
    assert depth_at[59] == 2
    assert depth_at[60] == 1
    assert depth_at[99] == 1
    assert depth_at[100] == 0
    assert depth_at[199] == 0


def test_bam_to_raw_df_indexed_intervals_cover_full_ref(synthetic_indexed_bam):
    """Intervals for each ref should cover from 0 to ref_length with no gaps."""
    df = bam_to_raw_df(synthetic_indexed_bam)

    for ref_name, ref_len in [("refA", 200), ("refB", 100)]:
        ref_df = df.filter(pl.col("ref") == ref_name).sort("start")
        assert ref_df["start"][0] == 0, f"{ref_name} should start at 0"
        assert ref_df["end"][-1] == ref_len, f"{ref_name} should end at {ref_len}"
        starts = ref_df["start"].to_list()
        ends = ref_df["end"].to_list()
        for i in range(1, len(starts)):
            assert starts[i] == ends[i - 1], f"Gap in {ref_name} at position {ends[i-1]}"


def test_indexed_and_unindexed_match(synthetic_bam, synthetic_indexed_bam):
    """Indexed and unindexed BAMs with same data produce identical DataFrames."""
    df_unindexed = bam_to_raw_df(synthetic_bam).sort("ref", "start")
    df_indexed = bam_to_raw_df(synthetic_indexed_bam).sort("ref", "start")

    assert df_unindexed.shape == df_indexed.shape
    assert df_unindexed.equals(df_indexed)
