from pathlib import Path

import pytest

from bam2plot.main import (
    PE_reads_to_df,
    _load_single_reference_sequence,
    _refresh_coverage_stats,
    main_guci,
)


def _write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    content = "".join(
        f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n"
        for name, sequence in records
    )
    path.write_text(content)


def test_pe_reads_to_df_pairs_reads_by_name_when_files_are_out_of_order(tmp_path):
    read_1 = tmp_path / "reads_1.fastq"
    read_2 = tmp_path / "reads_2.fastq"

    _write_fastq(read_1, [("pair2/1", "CCCC"), ("pair1/1", "AAAA")])
    _write_fastq(read_2, [("pair1/2", "TTTT"), ("pair2/2", "GGGG")])

    result = PE_reads_to_df(str(read_1), str(read_2))

    assert result.rows() == [("AAAA", "TTTT"), ("CCCC", "GGGG")]


def test_pe_reads_to_df_rejects_unmatched_pairs(tmp_path):
    read_1 = tmp_path / "reads_1.fastq"
    read_2 = tmp_path / "reads_2.fastq"

    _write_fastq(read_1, [("pair1/1", "AAAA"), ("pair2/1", "CCCC")])
    _write_fastq(read_2, [("pair1/2", "TTTT")])

    with pytest.raises(SystemExit):
        PE_reads_to_df(str(read_1), str(read_2))


def test_load_single_reference_sequence_rejects_multi_record_reference(tmp_path):
    reference = tmp_path / "reference.fa"
    reference.write_text(">ref1\nAAAA\n>ref2\nCCCC\n")

    with pytest.raises(SystemExit):
        _load_single_reference_sequence(str(reference))


def test_refresh_coverage_stats_recomputes_global_totals_for_filtered_subset(single_ref_df):
    refreshed = _refresh_coverage_stats(single_ref_df, threshold=10)

    assert refreshed["mean_coverage_total"][0] == pytest.approx(2500 / 300)
    assert refreshed["median_coverage_total"][0] == pytest.approx(5.0)
    assert refreshed["pct_total_over_zero"][0] == pytest.approx(200 / 300)
    assert refreshed["pct_total_over_thresh"][0] == pytest.approx(100 / 300)


def test_main_guci_plot_type_both_writes_png_and_svg(tmp_path):
    reference = tmp_path / "reference.fa"
    output_dir = tmp_path / "plots"
    reference.write_text(">ref1\nACGTACGT\n")

    with pytest.raises(SystemExit) as exc_info:
        main_guci(str(reference), window=2, out_folder=str(output_dir), plot_type="both")

    assert exc_info.value.code == 0
    assert (output_dir / "gc_reference.png").exists()
    assert (output_dir / "gc_reference.svg").exists()
