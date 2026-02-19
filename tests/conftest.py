import matplotlib

matplotlib.use("Agg")

import polars as pl
import numpy as np
import pytest

from bam2plot.main import enrich_coverage_df


@pytest.fixture
def mosdepth_df():
    """Two refs: refA (300bp) and refB (80bp) with varied depths.

    Uses enrich_coverage_df() to produce the enriched schema so the fixture
    always matches the current column set.
    """
    rows_a = [
        ("refA", 0, 100, 0),
        ("refA", 100, 200, 5),
        ("refA", 200, 300, 20),
    ]
    rows_b = [
        ("refB", 0, 40, 15),
        ("refB", 40, 80, 0),
    ]
    rows = rows_a + rows_b

    raw = pl.DataFrame(
        rows, schema=["ref", "start", "end", "depth"], orient="row"
    )

    return enrich_coverage_df(raw, thresh=10)


@pytest.fixture
def single_ref_df(mosdepth_df):
    """mosdepth_df filtered to refA only."""
    return mosdepth_df.filter(pl.col("ref") == "refA")


@pytest.fixture
def gc_sequence_df():
    """10-base sequence ATGCGCATAT with name/sequence/position columns."""
    seq = list("ATGCGCATAT")
    return pl.DataFrame(
        {
            "name": ["test"] * len(seq),
            "sequence": seq,
            "position": list(range(1, len(seq) + 1)),
        }
    )


@pytest.fixture
def plotting_ready_df():
    """Small DataFrame suitable for coverage_plot()."""
    n = 50
    pos = list(range(n))
    rolling = [0.0] * 10 + [5.0] * 20 + [15.0] * 20
    from bam2plot.main import COLOR_ZERO, COLOR_LOW, COLOR_HIGH

    color = (
        [COLOR_ZERO] * 10 + [COLOR_LOW] * 20 + [COLOR_HIGH] * 20
    )
    return pl.DataFrame(
        {
            "pos": pos,
            "rolling": rolling,
            "color": color,
            "ref": ["testref"] * n,
            "pct_over_thresh": [0.4] * n,
        }
    )


@pytest.fixture
def reads_coverage_df():
    """100 rows with pos/depth/rolling for plot_from_reads()."""
    n = 100
    return pl.DataFrame(
        {
            "pos": list(range(1, n + 1)),
            "depth": [float(i % 20) for i in range(n)],
            "rolling": [float(i % 20) for i in range(n)],
        }
    )
