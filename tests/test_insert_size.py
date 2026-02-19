import matplotlib

matplotlib.use("Agg")

import numpy as np
import pysam
import pytest
from matplotlib.figure import Figure

from bam2plot.main import extract_insert_sizes, plot_insert_size_distribution


@pytest.fixture
def paired_end_bam(tmp_path):
    """BAM with 2 proper pairs (template_length 300 and 400)."""
    bam_path = str(tmp_path / "paired.bam")

    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.0", "SO": "coordinate"},
        "SQ": [{"SN": "ref1", "LN": 1000}],
    })

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        # Pair 1 — read1
        a = pysam.AlignedSegment(outf.header)
        a.query_name = "pair1"
        a.query_sequence = "A" * 50
        a.flag = 0x1 | 0x2 | 0x40  # paired, proper pair, read1
        a.reference_id = 0
        a.reference_start = 0
        a.mapping_quality = 60
        a.cigar = [(0, 50)]
        a.query_qualities = pysam.qualitystring_to_array("I" * 50)
        a.template_length = 300
        a.next_reference_id = 0
        a.next_reference_start = 250
        outf.write(a)

        # Pair 1 — read2
        a = pysam.AlignedSegment(outf.header)
        a.query_name = "pair1"
        a.query_sequence = "A" * 50
        a.flag = 0x1 | 0x2 | 0x80  # paired, proper pair, read2
        a.reference_id = 0
        a.reference_start = 250
        a.mapping_quality = 60
        a.cigar = [(0, 50)]
        a.query_qualities = pysam.qualitystring_to_array("I" * 50)
        a.template_length = -300
        a.next_reference_id = 0
        a.next_reference_start = 0
        outf.write(a)

        # Pair 2 — read1
        a = pysam.AlignedSegment(outf.header)
        a.query_name = "pair2"
        a.query_sequence = "A" * 50
        a.flag = 0x1 | 0x2 | 0x40  # paired, proper pair, read1
        a.reference_id = 0
        a.reference_start = 500
        a.mapping_quality = 60
        a.cigar = [(0, 50)]
        a.query_qualities = pysam.qualitystring_to_array("I" * 50)
        a.template_length = 400
        a.next_reference_id = 0
        a.next_reference_start = 850
        outf.write(a)

        # Pair 2 — read2
        a = pysam.AlignedSegment(outf.header)
        a.query_name = "pair2"
        a.query_sequence = "A" * 50
        a.flag = 0x1 | 0x2 | 0x80  # paired, proper pair, read2
        a.reference_id = 0
        a.reference_start = 850
        a.mapping_quality = 60
        a.cigar = [(0, 50)]
        a.query_qualities = pysam.qualitystring_to_array("I" * 50)
        a.template_length = -400
        a.next_reference_id = 0
        a.next_reference_start = 500
        outf.write(a)

    return bam_path


@pytest.fixture
def single_end_bam(tmp_path):
    """BAM with only single-end reads (no paired-end data)."""
    bam_path = str(tmp_path / "single.bam")

    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.0", "SO": "coordinate"},
        "SQ": [{"SN": "ref1", "LN": 500}],
    })

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        a = pysam.AlignedSegment(outf.header)
        a.query_name = "read1"
        a.query_sequence = "A" * 100
        a.flag = 0  # single-end
        a.reference_id = 0
        a.reference_start = 0
        a.mapping_quality = 60
        a.cigar = [(0, 100)]
        a.query_qualities = pysam.qualitystring_to_array("I" * 100)
        outf.write(a)

    return bam_path


def test_extract_insert_sizes_paired(paired_end_bam):
    sizes = extract_insert_sizes(paired_end_bam)
    assert sizes is not None
    assert sorted(sizes.tolist()) == [300, 400]


def test_extract_insert_sizes_no_double_counting(paired_end_bam):
    """Only read1 of each pair should be counted (not read2)."""
    sizes = extract_insert_sizes(paired_end_bam)
    # 2 pairs → exactly 2 insert sizes
    assert len(sizes) == 2


def test_extract_insert_sizes_single_end(single_end_bam):
    result = extract_insert_sizes(single_end_bam)
    assert result is None


def test_plot_insert_size_returns_figure():
    sizes = np.array([200, 250, 300, 350, 400])
    fig = plot_insert_size_distribution(sizes, "test")
    assert isinstance(fig, Figure)
