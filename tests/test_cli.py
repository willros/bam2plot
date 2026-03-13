from pathlib import Path
import json
import subprocess
import sys

import pysam
from bam2plot import __version__

REPO_ROOT = Path(__file__).resolve().parents[1]


def _run_cli(*args: str):
    result = subprocess.run(
        [sys.executable, "-m", "bam2plot", *args],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    return result


def test_package_cli_supports_version_flag():
    result = subprocess.run(
        [sys.executable, "-m", "bam2plot", "--version"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert result.stdout.strip() == __version__


def _write_fastq(path: Path, name: str, sequence: str) -> None:
    path.write_text(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def _write_reference(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))


def _write_bam(path: Path) -> None:
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.0", "SO": "coordinate"},
            "SQ": [{"SN": "ref1", "LN": 200}],
        }
    )

    with pysam.AlignmentFile(str(path), "wb", header=header) as outf:
        for query_name, start in [("read1", 10), ("read2", 50)]:
            segment = pysam.AlignedSegment(outf.header)
            segment.query_name = query_name
            segment.query_sequence = "A" * 50
            segment.flag = 0
            segment.reference_id = 0
            segment.reference_start = start
            segment.mapping_quality = 60
            segment.cigar = [(0, 50)]
            segment.query_qualities = pysam.qualitystring_to_array("I" * 50)
            outf.write(segment)


def test_from_bam_cli_writes_summary_and_zoom_aware_report(tmp_path):
    bam_path = tmp_path / "cli.bam"
    out_dir = tmp_path / "from_bam_out"
    _write_bam(bam_path)

    _run_cli(
        "from_bam",
        "-b",
        str(bam_path),
        "-o",
        str(out_dir),
        "--threads",
        "1",
        "-n",
        "1",
        "-z",
        "50 60",
    )

    summary = json.loads((out_dir / "cli_summary.json").read_text())
    report = (out_dir / "cli_report.html").read_text()

    assert summary["analysis_scope"] == "zoom"
    assert summary["zoom"] == {"start": 50, "end": 60}
    assert summary["references"][0]["total_bases"] == 11
    assert "Region:</strong> 50-60" in report
    assert (out_dir / "cli_summary.tsv").exists()


def test_from_reads_cli_generates_plot(tmp_path):
    reference = tmp_path / "reference.fa"
    reads = tmp_path / "reads.fastq"
    out_dir = tmp_path / "from_reads_out"

    reference_sequence = "ACGTTGCATGCA" * 50
    _write_reference(reference, [("ref1", reference_sequence)])
    _write_fastq(reads, "read1", reference_sequence[100:300])

    _run_cli(
        "from_reads",
        "-r1",
        str(reads),
        "-ref",
        str(reference),
        "-o",
        str(out_dir),
    )

    assert (out_dir / "reads_bam2plot_ref1.png").exists()


def test_guci_cli_supports_multi_contig_reference(tmp_path):
    reference = tmp_path / "reference.fa"
    out_dir = tmp_path / "guci_out"
    _write_reference(
        reference,
        [("ref1", "ACGTACGTACGT"), ("ref2", "GGGGCCCCAAAA")],
    )

    _run_cli(
        "guci",
        "-ref",
        str(reference),
        "-w",
        "3",
        "-o",
        str(out_dir),
    )

    output = out_dir / "gc_reference.png"
    assert output.exists()
    assert output.stat().st_size > 0
