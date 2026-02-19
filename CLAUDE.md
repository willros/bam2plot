# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

bam2plot is a bioinformatics CLI tool for generating coverage plots from BAM files. It uses mosdepth for coverage extraction, Polars for data manipulation, and matplotlib/seaborn for visualization. Published on PyPI and citable via Zenodo.

## Build & Development Commands

```bash
# Install from source (editable)
pip install -e .

# Build distribution
make build          # runs python -m build

# Upload to PyPI
make upload         # runs twine upload dist/*

# Format code
make format         # runs black bam2plot/*.py

# Clean build artifacts
make clean          # removes dist/, build/, *.egg-info
```

**External runtime dependency:** mosdepth must be installed separately:
```bash
conda install -c bioconda mosdepth
```

There is no test suite. Manual testing can be done with the sample BAM file:
```bash
bam2plot from_bam -b test/barcode87.bam -o /tmp/test_output -s -c
```

## Architecture

The entire implementation lives in a single file: `bam2plot/main.py` (~1000 lines). Entry point is `cli()` which is registered as the `bam2plot` console script in setup.py.

### Three subcommands with distinct pipelines

1. **`from_bam`** — BAM → mosdepth → Polars DataFrame → coverage plot
   - `bam2plot_from_bam()` → `main_from_bam()` parses args, optionally sorts/indexes BAM with pysam, runs mosdepth subprocess, converts output via `mosdepth_to_df()`, plots top N references

2. **`from_reads`** — FASTQ + reference FASTA → mappy alignment → coverage DataFrame → plot
   - `bam2plot_from_reads()` → `main_from_reads()` uses mappy (minimap2 wrapper) for alignment, supports both single-end (`map_fastq_to_ref_long_read()`) and paired-end (`map_fastq_to_ref_PE_read()`)

3. **`guci`** — Reference FASTA → per-base GC content → rolling mean → plot
   - `bam2plot_guci()` → `main_guci()` uses pyfastx to read sequences, calculates GC via `add_gc()`

### Key data flow patterns

- Coverage data is always represented as Polars DataFrames with columns for position and depth
- `mosdepth_to_df()` is the bridge between mosdepth output and the plotting pipeline
- `coverage_plot()` uses matplotlib LineCollection with a three-color scheme: red (0X), yellow (below threshold), blue (above threshold)
- `return_ref_for_plotting()` applies rolling mean smoothing before visualization
- `plot_cumulative_coverage_for_all()` uses seaborn FacetGrid for multi-reference subplots

### Dependencies (pinned in setup.py)

pysam, seaborn, polars, mappy, pyfastx, pyarrow, numpy — all pinned to specific versions. Python >= 3.10 required.

### Version note

`__init__.py` contains `__version__` which may lag behind the version in `setup.py` (currently 0.3.7). Both should be updated together on version bumps.
