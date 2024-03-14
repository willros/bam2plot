from pathlib import Path
import subprocess
import pysam
import pandas as pd
from io import StringIO
import _io
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import os
import platform
import pandas
import matplotlib.ticker as mtick
import numpy as np
import argparse

import polars as pl

# for windows users
if platform.system() == "Windows":
    matplotlib.use("Agg")

sns.set_theme()


SORTED_TEMP = "TEMP112233.sorted.bam"
SORTED_TEMP_INDEX = f"{SORTED_TEMP}.bai"


class bcolors:
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    UNDERLINE = "\033[4m"


def print_green(text):
    print(f"{bcolors.OKGREEN}{text}{bcolors.ENDC}")


def print_warning(text):
    print(f"{bcolors.WARNING}{bcolors.UNDERLINE}{text}{bcolors.ENDC}")


def print_fail(text):
    print(f"{bcolors.FAIL}{bcolors.UNDERLINE}{text}{bcolors.ENDC}")


def print_blue(text):
    print(f"{bcolors.OKBLUE}{text}{bcolors.ENDC}")


def sort_bam(bam: str, new_name: str) -> None:
    try:
        pysam.sort("-o", new_name, bam)
    except:
        print_fail("[ERROR]: The file is not a bam file!")
        exit(1)


def index_bam(bam: str, new_name: str) -> None:
    try:
        pysam.index(bam, new_name)
    except:
        print_fail("[ERROR]: The file is not sorted! Run 'bam2plot <bam> -s'")
        exit(1)


def run_perbase(bam: str) -> _io.StringIO:
    return StringIO(
        subprocess.check_output(
            f"perbase only-depth {bam}", shell=True, stderr=subprocess.DEVNULL
        ).decode(errors="ignore")
    )


def perbase_to_df(perbase: _io.StringIO) -> pl.DataFrame:
    return (
        pl.read_csv(
            perbase,
            separator="\t",
        )
        .with_columns(pos=pl.int_ranges("POS", "END"))
        .explode("pos")
        .rename({"DEPTH": "depth", "REF": "ref"})
        .select(["ref", "pos", "depth"])
    )


def print_coverage_info(df, threshold: int) -> None:
    name = df["ref"][0]
    over_zero = df.with_columns(over_zero=pl.col("depth") > 0)["over_zero"].mean() * 100
    over_treshold = (
        df.with_columns(over_t=pl.col("depth") > threshold)["over_t"].mean() * 100
    )
    print_blue(f"[SUMMARIZE]: Coverage information for: {name}")

    print_blue(f"[SUMMARIZE]: Percent bases with coverage above 0X: {over_zero: .1f}%")

    print_blue(
        f"[SUMMARIZE]: Percent bases with coverage above {threshold}X: {over_treshold: .1f}%"
    )

    print_blue(f"   [SUMMARIZE]: median coverage: {df['depth'].median(): .0f}X")
    print_blue(f"   [SUMMARIZE]: mean coverage: {df['depth'].mean(): .0f}X")


def print_total_reference_info(df, threshold: int) -> None:
    mean_coverage = df["depth"].mean()

    over_zero = df.with_columns(over_zero=pl.col("depth") > 0)["over_zero"].mean() * 100
    over_treshold = (
        df.with_columns(over_t=pl.col("depth") > threshold)["over_t"].mean() * 100
    )

    print_blue(f"[SUMMARIZE]: Mean coverage of all basepairs: {mean_coverage: .1f}X")

    print_blue(f"[SUMMARIZE]: Percent bases with coverage above 0X: {over_zero: .1f}%")

    print_blue(
        f"[SUMMARIZE]: Percent bases with coverage above {threshold}X: {over_treshold: .1f}%"
    )


def under_threshold(df, threshold):
    df = df.assign(zero=lambda x: x.depth < threshold)

    start = []
    stop = []

    start_value = False
    for row in df.itertuples():
        if row.zero and not start_value:
            start_value = True
            start.append(row.pos)

        if not row.zero and start_value:
            start_value = False
            stop.append(row.pos)

    if len(start) > len(stop):
        stop.append(df.pos.max())
    return start, stop


def plot_coverage(
    df,
    sample_name: str,
    threshold: int,
    rolling_window: int,
    log_scale: bool = False,
    highlight: bool = False,
) -> matplotlib.figure.Figure:
    if log_scale:
        df = df.with_columns(depth=(pl.col("depth") + 1).log10()).with_columns(
            rolling=(pl.col("depth") + 1).log10()
        )

        threshold = np.log10(threshold)

    mean_coverage = df["depth"].mean()

    over_treshold = (
        df.with_columns(over_t=pl.col("depth") > threshold)["over_t"].mean() * 100
    )

    df = df.to_pandas()

    coverage_plot = plt.figure(figsize=(15, 8))
    ax = sns.lineplot(data=df, x="pos", y="rolling")
    zero = plt.axhline(y=0, color="red")
    zero.set_label("Zero")
    mean = plt.axhline(y=mean_coverage, color="green")
    mean.set_label(f"Mean coverage: {mean_coverage: .1f}X")
    plt.legend(loc="upper right")
    plt.title(
        f"Percent bases with coverage above {threshold}X: {over_treshold: .1f}% | Rolling window: {rolling_window} nt"
    )
    plt.suptitle(f"Ref: {df.iloc[0].ref} | Sample: {sample_name}")
    ax.set(xlabel="Position", ylabel="Depth")

    if highlight:
        for a, b in zip(*under_threshold(df, threshold)):
            plt.fill_between([a, b], 0, mean_coverage, color="red", alpha=0.2)

    plt.close()

    return coverage_plot


def coverage_for_value(df, coverage: int):
    number_of_bases = df["pos"].max()
    _id = df["ref"][0][:20]

    percent = df.with_columns(over_t=pl.col("depth") > coverage)["over_t"].mean() * 100

    return pl.DataFrame(
        {
            "coverage": [coverage],
            "percent": [percent],
            "id": [_id],
        }
    )


def coverage_for_many_values(df, values):
    dfs = []
    for coverage in values:
        coverage_df = coverage_for_value(df, coverage)
        dfs.append(coverage_df)
    return pl.concat(dfs)


def plot_cumulative_coverage_for_all(df):
    max_cov = df["depth"].max()
    coverage_values = np.linspace(0, max_cov, 15)

    all_coverage = pl.concat(
        [
            coverage_for_many_values(df.filter(pl.col("ref") == ref), coverage_values)
            for ref in df["ref"].unique()
        ]
    )
    all_coverage = all_coverage.to_pandas()
    grid = sns.FacetGrid(all_coverage, col="id", height=2.5, col_wrap=5)
    grid.map_dataframe(sns.lineplot, x="coverage", y="percent")
    plt.close()
    return grid.fig


def make_dir(outpath: str) -> None:
    outpath = Path(outpath)
    if not outpath.exists():
        outpath.mkdir(parents=True)


def cli():
    parser = argparse.ArgumentParser(description="Plot your bam files!")
    parser.add_argument("-b", "--bam", required=True, help="bam file")
    parser.add_argument(
        "-o",
        "--outpath",
        required=False,
        default="bam2plots",
        help="Where to save the plots.",
    )
    parser.add_argument(
        "-w",
        "--whitelist",
        required=False,
        default=None,
        help="Only include these references/chromosomes.",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        required=False,
        default=3,
        help="Threshold of mean coverage depth",
        type=int,
    )
    parser.add_argument(
        "-r",
        "--rolling_window",
        required=False,
        default=100,
        help="Rolling window size",
        type=int,
    )
    parser.add_argument(
        "-i",
        "--index",
        required=False,
        default=False,
        help="Index bam file",
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "-s",
        "--sort_and_index",
        required=False,
        default=False,
        help="Index and sort bam file",
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "-z",
        "--zoom",
        required=False,
        default=False,
        help="Zoom into this region. Example: -z='100 2000'",
    )
    parser.add_argument(
        "-l",
        "--log_scale",
        required=False,
        default=False,
        help="Log scale of Y axis",
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "-c",
        "--cum_plot",
        required=False,
        default=False,
        help="Generate cumulative plots of all chromosomes",
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "-hl",
        "--highlight",
        required=False,
        default=False,
        help="Highlights regions where coverage is below treshold.",
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "-p",
        "--plot_type",
        required=False,
        default="png",
        choices=["png", "svg", "both"],
        help="How to save the plots",
    )

    args = parser.parse_args()
    command = "\nbam2plot \n" + "".join(f"{k}: {v}\n" for k, v in vars(args).items())
    print_green(command)

    main(
        bam=args.bam,
        outpath=args.outpath,
        whitelist=args.whitelist,
        rolling_window=args.rolling_window,
        threshold=args.threshold,
        index=args.index,
        sort_and_index=args.sort_and_index,
        zoom=args.zoom,
        log_scale=args.log_scale,
        cum_plot=args.cum_plot,
        highlight=args.highlight,
        plot_type=args.plot_type,
    )


def if_sort_and_index(sort_and_index, index, bam):
    if not sort_and_index and not index:
        return run_perbase(bam)

    if sort_and_index:
        print_green("[INFO]: Sorting bam file")
        sort_bam(bam, new_name=SORTED_TEMP)
        print_green("[INFO]: Indexing bam file")
        index_bam(SORTED_TEMP, new_name=SORTED_TEMP_INDEX)
    if index:
        print_green("[INFO]: Indexing bam file")
        index_name = f"{bam}.bai"
        index_bam(bam, new_name=index_name)

    try:
        if sort_and_index:
            perbase = run_perbase(SORTED_TEMP)
        if index:
            perbase = run_perbase(bam)
    except:
        print_fail("[ERROR]: Could not run perbase on bam file")
        exit(1)
    finally:
        if sort_and_index:
            os.remove(SORTED_TEMP)
            os.remove(SORTED_TEMP_INDEX)
        if index:
            os.remove(index_name)

        return perbase


def process_dataframe(perbase, sort_and_index, index):
    try:
        print_green("[INFO]: Processing dataframe")
        df = perbase_to_df(perbase)
    except:
        print_fail("[ERROR]: Could not process dataframe")
        if not sort_and_index:
            print_warning(
                "[WARNING]: Is the file indexed? If not, run 'bam2plot <file.bam> -i'"
            )
            print_warning(
                "[WARNING]: Is the file sorted? If not, run 'bam2plot <file.bam> -s'"
            )
            exit(1)
        if not index:
            print_warning(
                "[WARNING]: Is the file indexed? If not, run 'bam2plot <file.bam> -i'"
            )
            exit(1)

    return df


def save_plot_coverage(plot, outpath, sample_name, reference, plot_type):
    out_file = f"{outpath}/{sample_name}_bam2plot"
    name = f"{out_file}_{reference}"

    if plot_type == "png":
        plot.savefig(f"{name}.png")

    if plot_type == "svg":
        plot.savefig(f"{name}.svg")

    if plot_type == "both":
        plot.savefig(f"{name}.svg")
        plot.savefig(f"{name}.png")

    print_green(f"[INFO]: Plot for {reference} generated")


def save_plot_cum(cum_plot, outpath, bam, plot_type):
    cum_plot_name = f"{outpath}/{Path(bam).stem}_cumulative_coverage"
    print_green(f"[INFO]: Cumulative plot generated!")

    if plot_type == "png":
        cum_plot.savefig(f"{cum_plot_name}.png")

    if plot_type == "svg":
        cum_plot.savefig(f"{cum_plot_name}.svg")

    if plot_type == "both":
        cum_plot.savefig(f"{cum_plot_name}.png")
        cum_plot.savefig(f"{cum_plot_name}.svg")

    print_green(f"[INFO]: Cumulative plot generated!")


def main(
    bam,
    outpath,
    whitelist,
    rolling_window,
    threshold,
    index,
    sort_and_index,
    zoom,
    log_scale,
    cum_plot,
    highlight,
    plot_type,
) -> None:
    print_green(f"[INFO]: Running bam2plot on {bam}!")
    if zoom:
        start = int(zoom.split(" ")[0])
        end = int(zoom.split(" ")[1])
        if start >= end:
            print_fail("[ERROR]: Start value of zoom must be lower than end value.")
            exit(1)

    if not Path(bam).exists():
        print_fail(f"[ERROR]: The file {bam} does not exist")
        exit(1)

    make_dir(outpath)
    sample_name = Path(bam).stem

    perbase = if_sort_and_index(sort_and_index, index, bam)

    df = process_dataframe(perbase, sort_and_index, index)

    print_total_reference_info(df, threshold)

    if whitelist:
        whitelist = [whitelist] if type(whitelist) == str else whitelist
        print_green(
            f"[INFO]: Only looking for references in the whitelist: {whitelist}"
        )
        df = df.filter(pl.col("ref").is_in(whitelist))

    plot_number = df["ref"].n_unique()
    if plot_number == 0:
        print_fail("[ERROR]: No reference to plot against!")
        exit(1)

    plot_text = "plot" if plot_number == 1 else "plots"
    print_green(f"[INFO]: Generating {plot_number} {plot_text}:")

    for reference in df["ref"].unique():
        mpileup_df = df.filter(pl.col("ref") == reference).with_columns(
            rolling=pl.col("depth").rolling_mean(window_size=rolling_window)
        )

        if zoom:
            mpileup_df = mpileup_df.filter(pl.col("pos").is_between(start, end))
            if mpileup_df.shape[0] == 0:
                print_warning("[WARNING]: No positions to plot after zoom")
                continue

        if mpileup_df.shape[0] == 0:
            print_warning("[WARNING]: No positions to plot")
            continue

        print_coverage_info(mpileup_df, threshold)

        plot = plot_coverage(
            mpileup_df,
            sample_name,
            threshold=threshold,
            rolling_window=rolling_window,
            log_scale=log_scale,
            highlight=highlight,
        )

        save_plot_coverage(plot, outpath, sample_name, reference, plot_type)

    print_green("[INFO]: Coverage plots done!")

    if cum_plot:
        print_green("[INFO]: Generating cumulative coverage plots for each reference")
        cum_plot = plot_cumulative_coverage_for_all(df)
        save_plot_cum(cum_plot, outpath, bam, plot_type)

    print_green(f"[INFO]: Plots location: {Path(outpath).resolve()}")
    exit(0)


if __name__ == "__main__":
    cli()
