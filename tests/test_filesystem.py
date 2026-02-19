import pytest
from pathlib import Path
from bam2plot.main import make_dir, files_not_exists


def test_make_dir_creates_directory(tmp_path):
    new_dir = tmp_path / "new_output"
    make_dir(str(new_dir))
    assert new_dir.exists()
    assert new_dir.is_dir()


def test_make_dir_existing_no_error(tmp_path):
    existing = tmp_path / "existing"
    existing.mkdir()
    # Should not raise
    make_dir(str(existing))
    assert existing.exists()


def test_make_dir_nested(tmp_path):
    nested = tmp_path / "a" / "b" / "c"
    make_dir(str(nested))
    assert nested.exists()
    assert nested.is_dir()


def test_files_not_exists_existing_file(tmp_path):
    f = tmp_path / "real_file.txt"
    f.write_text("hello")
    # Should not raise
    files_not_exists(str(f))


def test_files_not_exists_missing_file():
    with pytest.raises(SystemExit):
        files_not_exists("/nonexistent/path/to/file.bam")
