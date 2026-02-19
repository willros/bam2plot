# Benchmarks

Coverage computation performance of `bam_to_raw_df()` (parallel pysam sweep-line) compared to the previous mosdepth-based pipeline.

## Setup

- **Platform**: Linux 6.8.0, 72 CPUs
- **Python**: 3.10.14
- **pysam**: 0.22.0, **mosdepth**: 0.3.10
- **BAMs**: synthetic, sorted and indexed, uniform read distribution (100-300 bp reads)

## Results

### 500k reads, 5 references x 500 kb (2.5 Mb total)

| Method | Median | Min | vs mosdepth |
|---|---|---|---|
| **bam_to_raw_df (parallel)** | **0.552s** | 0.534s | **1.3x faster** |
| mosdepth end-to-end | 0.723s | 0.632s | baseline |
| sequential sweep-line | 0.833s | 0.809s | 1.15x slower |

### 2M reads, 10 references x 1 Mb (10 Mb total)

| Method | Median | Min | vs mosdepth |
|---|---|---|---|
| **bam_to_raw_df (parallel)** | **1.684s** | 1.631s | **1.6x faster** |
| mosdepth end-to-end | 2.748s | 2.721s | baseline |
| sequential sweep-line | 3.307s | 3.248s | 1.21x slower |

Output is bit-identical to mosdepth across all tests.

## Approaches evaluated

Several strategies were benchmarked before settling on the parallel sweep-line:

| Approach | 500k reads | Notes |
|---|---|---|
| `pysam.count_coverage()` | 7.93s | 10x slower; does per-nucleotide (A/C/G/T) counting, overkill for depth |
| `pysam.depth()` + Polars parse | 0.80s | Tied with mosdepth, but adds string/file I/O overhead |
| `np.add.at()` batch sweep | 1.07s | List-building negates batch benefit vs inline indexing |
| Sequential sweep (flag bitmask) | 0.83s | 15% slower than mosdepth; Python iteration is the bottleneck |
| **Parallel sweep (multiprocessing)** | **0.55s** | Each worker processes one reference via `pysam.fetch(contig)` |

## Why parallel wins

The bottleneck in pure-Python coverage computation is iterating reads via pysam (each read creates a Python `AlignedSegment` wrapper). Parallelizing across references with `multiprocessing.Pool` distributes this cost across CPUs. Each worker independently opens the BAM, seeks to its assigned reference via the index, and builds a coverage array. The per-worker overhead (process creation, numpy serialization) is small relative to the per-reference read iteration time.

## Reproducibility

Benchmarks can be reproduced by creating synthetic BAMs with known read distributions and timing `bam_to_raw_df()` from `bam2plot.main`. See the git history for the exact benchmark scripts used during development.
