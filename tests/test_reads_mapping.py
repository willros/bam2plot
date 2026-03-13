import numpy as np

from bam2plot.main import _add_alignment_to_array


def test_add_alignment_to_array_converts_half_open_interval_to_one_based_positions():
    diff = np.zeros(7, dtype=np.int64)

    _add_alignment_to_array(diff, start=2, end=5)

    assert np.cumsum(diff[:-1]).tolist() == [0, 0, 1, 1, 1, 0]


def test_add_alignment_to_array_handles_reference_start_correctly():
    diff = np.zeros(5, dtype=np.int64)

    _add_alignment_to_array(diff, start=0, end=2)

    assert np.cumsum(diff[:-1]).tolist() == [1, 1, 0, 0]
