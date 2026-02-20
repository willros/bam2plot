# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

bam2plot is a bioinformatics CLI tool for generating coverage plots and QC reports from BAM files. It uses a parallel pysam sweep-line algorithm for coverage computation, Polars for data manipulation, and matplotlib/seaborn for visualization. Published on PyPI and citable via Zenodo.

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

No external runtime dependencies required — all coverage computation is done via pysam.

Tests can be run with pytest:
```bash
pytest tests/
```

Manual testing can be done with the sample BAM file:
```bash
bam2plot from_bam -b test/barcode87.bam -o /tmp/test_output -s -c
```

## Architecture

The entire implementation lives in a single file: `bam2plot/main.py` (~1500 lines). Entry point is `cli()` which is registered as the `bam2plot` console script in setup.py.

### Three subcommands with distinct pipelines

1. **`from_bam`** — BAM → parallel pysam sweep-line → RLE DataFrame → coverage plot + QC reports + HTML report
   - `bam2plot_from_bam()` → `main_from_bam()` parses args, optionally sorts/indexes BAM with pysam, computes coverage via `bam_to_raw_df()`, enriches with `enrich_coverage_df()`, generates all plots and HTML report

2. **`from_reads`** — FASTQ + reference FASTA → mappy alignment → coverage DataFrame → plot
   - `bam2plot_from_reads()` → `main_from_reads()` uses mappy (minimap2 wrapper) for alignment, supports both single-end (`map_fastq_to_ref_long_read()`) and paired-end (`map_fastq_to_ref_PE_read()`)

3. **`guci`** — Reference FASTA → per-base GC content → rolling mean → plot
   - `bam2plot_guci()` → `main_guci()` uses pyfastx to read sequences, calculates GC via `add_gc()`

### Key data flow patterns

- Coverage data is always represented as Polars DataFrames with RLE-encoded columns (ref, start, end, depth)
- `bam_to_raw_df()` computes coverage via parallel sweep-line, outputs RLE DataFrame
- `enrich_coverage_df()` adds per-reference statistics: mean, median, Gini, percent above thresholds
- `coverage_plot()` uses matplotlib LineCollection with a three-color scheme: red (0X), yellow (below threshold), blue (above threshold)
- `return_ref_for_plotting()` applies rolling mean smoothing before visualization
- `plot_cumulative_coverage_for_all()` uses seaborn FacetGrid for multi-reference subplots
- `plot_depth_histogram()` and `plot_depth_histogram_global()` generate weighted histograms with mean/median lines
- `plot_lorenz_curves()` generates FacetGrid with Gini annotations per reference
- `generate_html_report()` produces a self-contained HTML file with all plots embedded as base64

### Dependencies (pinned in setup.py)

pysam, seaborn, polars, mappy, pyfastx, pyarrow, numpy, pandas — all pinned to specific versions. Python >= 3.10 required.

### Version note

`__init__.py` and `setup.py` both contain the version (currently 0.4.0). Both should be updated together on version bumps.
