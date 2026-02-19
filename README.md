# bam2plot

Generate coverage plots from BAM files. No external tools required.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15052225.svg)](https://zenodo.org/records/15052225)

## Installation

```bash
pip install bam2plot
```

Python >= 3.10 required. All dependencies (pysam, polars, matplotlib, etc.) are installed automatically.

## Quick start

```bash
# Plot coverage from a BAM file (sort + index automatically)
bam2plot from_bam -b input.bam -o output_folder -s

# With cumulative coverage plots and custom threshold
bam2plot from_bam -b input.bam -o output_folder -s -c -t 20

# Plot coverage directly from reads (no BAM needed)
bam2plot from_reads -r1 reads.fastq -ref reference.fasta -o output_folder

# Plot GC content of a reference
bam2plot guci -ref reference.fasta -w 1000 -o output_folder
```

## Output examples

Coverage plot with three-color depth visualization (red = 0X, yellow = below threshold, blue = above threshold):

![coverage plot](example/normal.png)

Cumulative coverage plots per reference (with `-c`):

![cumulative plot](example/cumplot.png)

Terminal output while running:

![terminal](example/running.png)

## Subcommands

### `from_bam` -- BAM to coverage plot

```
bam2plot from_bam -b BAM -o OUTPATH [options]

  -b, --bam              BAM file (required)
  -o, --outpath          Output directory (required)
  -s, --sort_and_index   Sort and index the BAM before plotting
  -i, --index            Index only (BAM must already be sorted)
  -t, --threshold        Coverage depth threshold (default: 10)
  -r, --rolling_window   Rolling window size for smoothing (default: 100)
  -n, --number_of_refs   How many references to plot (default: 10, max: 100)
  -w, --whitelist        Only plot these references
  -z, --zoom             Zoom into a region, e.g. -z='1000 5000'
  -c, --cum_plot         Also generate cumulative coverage plots
  -p, --plot_type        Output format: png, svg, or both (default: png)
```

### `from_reads` -- Reads + reference to coverage plot

Aligns reads to a reference using minimap2 (via mappy) and plots coverage. Supports both long reads and paired-end short reads.

```
bam2plot from_reads -r1 READ_1 -ref REFERENCE -o OUT_FOLDER [options]

  -r1, --read_1          FASTQ file (required)
  -r2, --read_2          Second FASTQ for paired-end reads
  -ref, --reference      Reference FASTA (required)
  -o, --out_folder       Output directory (required)
  -gc, --guci            Overlay GC content on the coverage plot
  -r, --rolling_window   Rolling window size (default: 50)
  -p, --plot_type        Output format: png, svg, or both (default: png)
```

### `guci` -- GC content plot

Computes per-base GC content of a reference FASTA with a rolling mean window.

```
bam2plot guci -ref REFERENCE -w WINDOW -o OUT_FOLDER [options]

  -ref, --reference      Reference FASTA (required)
  -w, --window           Rolling window size (required)
  -o, --out_folder       Output directory (required)
  -p, --plot_type        Output format: png, svg, or both (default: png)
```

## How it works

### Coverage computation (`from_bam`)

bam2plot computes per-base coverage depth directly from BAM files using pysam, with no external dependencies like mosdepth or samtools.

The pipeline has three stages:

**1. Sweep-line depth computation**

For each reference sequence, a coverage array is allocated (one element per base position plus a sentinel). For every aligned read, `+1` is added at the read's start position and `-1` at its end position. A cumulative sum over the array then yields the exact per-position depth. This is the classic sweep-line algorithm, running in O(reads + reference_length) time.

Secondary, QC-failed, and duplicate reads are filtered out (matching samtools/mosdepth defaults).

**2. Parallel dispatch for indexed BAMs**

When the BAM has an index (`.bai` file), bam2plot parallelizes across references using Python's `multiprocessing.Pool`. Each worker process opens the BAM independently, seeks to its assigned reference via `pysam.fetch(contig)`, and computes coverage for that reference alone. This eliminates the Python per-read iteration bottleneck by distributing it across CPUs -- benchmarked at 1.3-1.6x faster than mosdepth ([details](BENCHMARKS.md)).

For unindexed BAMs, a single-pass sequential sweep is used instead, iterating all reads once and dispatching to per-reference arrays by reference ID.

**3. Run-length encoding and enrichment**

The per-position depth array is compressed into run-length encoded (RLE) intervals: consecutive positions with identical depth are merged into `(ref, start, end, depth)` rows. This typically reduces millions of positions to hundreds of thousands of intervals, making downstream operations efficient.

The RLE DataFrame is then enriched with per-reference statistics: mean coverage, percentage of bases above zero, percentage above the user-specified threshold, and genome-wide totals.

**4. Visualization**

For each reference, the RLE intervals are exploded back to per-position depth, smoothed with a rolling mean, and rendered as a colored line plot using matplotlib's `LineCollection`. Three colors indicate depth status: red (0X), yellow (below threshold), blue (above threshold). Optionally, seaborn `FacetGrid` cumulative coverage plots are generated across references.

### Read alignment (`from_reads`)

The `from_reads` subcommand skips BAM files entirely. It aligns FASTQ reads directly to a reference FASTA using mappy (the Python binding for minimap2). Long reads use the `map-ont` preset; paired-end reads use the `sr` preset. Coverage is accumulated in a Polars DataFrame and plotted with the same visualization pipeline.

### GC content (`guci`)

The `guci` subcommand reads a reference FASTA with pyfastx, marks each base as GC (1) or AT (0), computes a rolling mean over the specified window, and plots the result.

## Citation

If you use bam2plot, please cite via [Zenodo](https://zenodo.org/records/15052225).

## License

MIT
