from pathlib import Path
import subprocess
import sys

from bam2plot import __version__


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_package_cli_supports_version_flag():
    result = subprocess.run(
        [sys.executable, "-m", "bam2plot", "--version"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert result.stdout.strip() == __version__
