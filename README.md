# BAM2PLOT

Plot your bam files!

bam2plot generates coverage plots:
![plot](example/normal.png)

... and if `-c` is added, cumulative coverage plots for each reference (*e.g.* chromosomes) for each sample:
![plot](example/cumplot.png)

If the flag `--highlight` is given, the regions with a coverage below the `--treshold` are highlighted:
![plot](example/highlight.png)

Below is an example of how bam2plot looks when runned in the terminal:
![plot](example/running.png)

## Dependencies
`bam2plot` depends on `perbase`, which you can install via:
```bash
cargo install perbase 
# or
conda install -c bioconda perbase
```
## Installation

You can install `bam2plot` using the following pip command:

```bash
pip install bam2plot
```

## Usage
Once installed, you can use the `bam2plot` command to perform coverage analysis on BAM files and generate coverage plots. Here's how to use it:

```bash
usage: bam2plot [-h] -b BAM [-o OUTPATH] [-w WHITELIST] [-t THRESHOLD] [-r ROLLING_WINDOW]
                [-i | --index | --no-index] [-s | --sort_and_index | --no-sort_and_index] [-z ZOOM]
                [-l | --log_scale | --no-log_scale] [-c | --cum_plot | --no-cum_plot]
                [-hl | --highlight | --no-highlight] [-p {png,svg,both}]

Plot your bam files!

options:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     bam file
  -o OUTPATH, --outpath OUTPATH
                        Where to save the plots.
  -w WHITELIST, --whitelist WHITELIST
                        Only include these references/chromosomes.
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold of mean coverage depth
  -r ROLLING_WINDOW, --rolling_window ROLLING_WINDOW
                        Rolling window size
  -i, --index, --no-index
                        Index bam file (default: False)
  -s, --sort_and_index, --no-sort_and_index
                        Index and sort bam file (default: False)
  -z ZOOM, --zoom ZOOM  Zoom into this region. Example: -z='100 2000'
  -l, --log_scale, --no-log_scale
                        Log scale of Y axis (default: False)
  -c, --cum_plot, --no-cum_plot
                        Generate cumulative plots of all chromosomes (default: False)
  -hl, --highlight, --no-highlight
                        Highlights regions where coverage is below treshold. (default: False)
  -p {png,svg,both}, --plot_type {png,svg,both}
                        How to save the plots
```

## Outputs

The BAM2Plot package generates coverage plots for each reference sequence in the BAM file. The plots are saved as SVG and PNG files in the specified output path or the default location.

## Examples

Here's an example of how to use the BAM2Plot package:

```bash
bam2plot --bam input.bam --outpath output_folder --rolling_window 50 --threshold 5 -s -c -hl
```
