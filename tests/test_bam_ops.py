import pytest
from unittest.mock import patch, MagicMock
from bam2plot.main import sort_bam, index_bam


def test_sort_bam_calls_pysam():
    with patch("bam2plot.main.pysam") as mock_pysam:
        sort_bam("input.bam", "output.sorted.bam")
        mock_pysam.sort.assert_called_once_with(
            "-o", "output.sorted.bam", "input.bam"
        )


def test_sort_bam_failure_exits():
    with patch("bam2plot.main.pysam") as mock_pysam:
        mock_pysam.sort.side_effect = Exception("broken bam")
        with pytest.raises(SystemExit):
            sort_bam("input.bam", "output.sorted.bam")


def test_index_bam_calls_pysam():
    with patch("bam2plot.main.pysam") as mock_pysam:
        index_bam("input.bam", "input.bam.bai")
        mock_pysam.index.assert_called_once_with(
            "input.bam", "input.bam.bai"
        )


def test_index_bam_failure_exits():
    with patch("bam2plot.main.pysam") as mock_pysam:
        mock_pysam.index.side_effect = Exception("not sorted")
        with pytest.raises(SystemExit):
            index_bam("input.bam", "input.bam.bai")
