import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import polars as pl
import pytest
import re

from bam2plot.main import generate_html_report, _fig_to_base64


@pytest.fixture
def sample_fig():
    """A minimal matplotlib figure for testing."""
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    plt.close(fig)
    return fig


@pytest.fixture
def html_report(tmp_path, mosdepth_df, sample_fig):
    """Generate an HTML report and return its path."""
    coverage_figures = [("refA", sample_fig), ("refB", sample_fig)]
    generate_html_report(
        sample_name="test_sample",
        bam="/path/to/test.bam",
        threshold=10,
        rolling_window=100,
        df=mosdepth_df,
        top_refs=["refA", "refB"],
        coverage_figures=coverage_figures,
        cum_fig=sample_fig,
        outpath=str(tmp_path),
    )
    return tmp_path / "test_sample_report.html"


def test_html_report_created(html_report):
    assert html_report.exists()
    assert html_report.stat().st_size > 0


def test_html_report_is_standalone(html_report):
    content = html_report.read_text()
    # No external stylesheet links
    assert "<link " not in content.lower()
    # No external script tags
    assert not re.search(r'<script\s+src=', content, re.IGNORECASE)


def test_html_report_contains_base64_images(html_report):
    content = html_report.read_text()
    matches = re.findall(r"data:image/png;base64,", content)
    # 2 coverage plots + 1 cumulative = 3
    assert len(matches) == 3


def test_html_report_contains_stats(html_report):
    content = html_report.read_text()
    assert "Mean coverage" in content
    assert "&gt; 0X" in content or "> 0X" in content
    assert "refA" in content
    assert "refB" in content
    assert "test_sample" in content


def test_html_report_no_cum_fig(tmp_path, mosdepth_df, sample_fig):
    """When cum_fig is None, no cumulative section appears."""
    coverage_figures = [("refA", sample_fig)]
    generate_html_report(
        sample_name="nocum",
        bam="/path/to/test.bam",
        threshold=10,
        rolling_window=100,
        df=mosdepth_df,
        top_refs=["refA"],
        coverage_figures=coverage_figures,
        cum_fig=None,
        outpath=str(tmp_path),
    )
    report = tmp_path / "nocum_report.html"
    content = report.read_text()
    matches = re.findall(r"data:image/png;base64,", content)
    # Only 1 coverage plot, no cumulative
    assert len(matches) == 1
    assert "Cumulative Coverage" not in content


def test_fig_to_base64_returns_valid_string(sample_fig):
    result = _fig_to_base64(sample_fig)
    assert isinstance(result, str)
    assert len(result) > 100
    # Should be valid base64 (no whitespace, only base64 chars)
    assert re.match(r'^[A-Za-z0-9+/=]+$', result)
