from pathlib import Path
import sys
import os
import platform
import argparse
import multiprocessing
import warnings
import base64
import io
from datetime import datetime

import polars as pl
import seaborn as sns
import numpy as np
import pysam
import mappy as mp
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import pyfastx


warnings.filterwarnings("ignore", category=FutureWarning)

# for windows users
if platform.system() == "Windows":
    matplotlib.use("Agg")

sns.set_theme()


SORTED_TEMP = "TEMP112233.sorted.bam"
SORTED_TEMP_INDEX = f"{SORTED_TEMP}.bai"

# Coverage plot color scheme
COLOR_ZERO = "#FF6F61"   # red: 0X coverage
COLOR_LOW = "#FFC107"    # yellow: below threshold
COLOR_HIGH = "#4A90E2"   # blue: above threshold


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
    except Exception as e:
        print_fail(f"[ERROR]: The bam file is broken or it is not a bam file! ({e})")
        sys.exit(1)


def index_bam(bam: str, new_name: str) -> None:
    try:
        pysam.index(bam, new_name)
    except Exception as e:
        print_fail(
            f"[ERROR]: The file is not sorted or it is not a bam file! Run 'bam2plot <bam> -s' ({e})"
        )
        sys.exit(1)


_FILTER_FLAGS = 0x4 | 0x100 | 0x200 | 0x400  # unmapped | secondary | qcfail | duplicate


def _sweep_one_ref(args):
    """Compute RLE coverage for a single reference (worker for parallel dispatch)."""
    bam, ref_name, ref_len = args
    af = pysam.AlignmentFile(bam, "rb")
    arr = np.zeros(ref_len + 1, dtype=np.int32)
    for read in af.fetch(ref_name):
        if read.flag & _FILTER_FLAGS:
            continue
        arr[read.reference_start] += 1
        arr[read.reference_end] -= 1
    af.close()
    depth = np.cumsum(arr[:-1])
    return _rle_encode(ref_name, depth)


def _rle_encode(ref_name, depth):
    """Run-length encode a per-position depth array."""
    if len(depth) == 0:
        return (ref_name, np.array([], dtype=np.int64), np.array([], dtype=np.int64), np.array([], dtype=np.int64))
    changes = np.diff(depth, prepend=depth[0] - 1)
    starts = np.where(changes != 0)[0]
    ends = np.append(starts[1:], len(depth))
    depths = depth[starts]
    return (ref_name, starts.astype(np.int64), ends.astype(np.int64), depths.astype(np.int64))


def _results_to_dataframe(results):
    """Concatenate per-reference RLE tuples into a single DataFrame."""
    all_rows = []
    for ref_name, starts, ends, depths in results:
        if len(starts) == 0:
            continue
        ref_col = np.full(len(starts), ref_name, dtype=object)
        all_rows.append((ref_col, starts, ends, depths))
    if not all_rows:
        return pl.DataFrame(
            schema={"ref": pl.Utf8, "start": pl.Int64, "end": pl.Int64, "depth": pl.Int64}
        )
    return pl.DataFrame({
        "ref": np.concatenate([r[0] for r in all_rows]),
        "start": np.concatenate([r[1] for r in all_rows]),
        "end": np.concatenate([r[2] for r in all_rows]),
        "depth": np.concatenate([r[3] for r in all_rows]),
    })


def bam_to_raw_df(bam: str) -> pl.DataFrame:
    """Compute per-base coverage from a BAM file.

    For indexed BAMs, parallelizes across references using multiprocessing.
    Falls back to sequential sweep-line for unindexed BAMs.
    Filters secondary, qcfail, and duplicate reads (matching mosdepth/samtools defaults).

    Returns a 4-column DataFrame: ref, start, end, depth (run-length encoded).
    """
    af = pysam.AlignmentFile(bam, "rb")
    ref_names = list(af.references)
    ref_lens = list(af.lengths)
    has_index = af.has_index()
    af.close()

    n_refs = len(ref_names)
    args = [(bam, ref_names[i], ref_lens[i]) for i in range(n_refs)]

    if has_index and n_refs > 1:
        n_workers = min(n_refs, os.cpu_count() or 1, 4)
        with multiprocessing.Pool(n_workers) as pool:
            results = pool.map(_sweep_one_ref, args)
    else:
        # Sequential: single pass for unindexed BAMs or single-ref BAMs
        if has_index:
            results = [_sweep_one_ref(a) for a in args]
        else:
            arrays = [np.zeros(ref_lens[i] + 1, dtype=np.int32) for i in range(n_refs)]
            af = pysam.AlignmentFile(bam, "rb")
            for read in af:
                if read.flag & _FILTER_FLAGS:
                    continue
                arrays[read.reference_id][read.reference_start] += 1
                arrays[read.reference_id][read.reference_end] -= 1
            af.close()
            results = [
                _rle_encode(ref_names[i], np.cumsum(arrays[i][:-1]))
                for i in range(n_refs)
            ]

    return _results_to_dataframe(results)


def enrich_coverage_df(raw_df: pl.DataFrame, thresh: int) -> pl.DataFrame:
    """Add coverage statistics columns to a raw (ref, start, end, depth) DataFrame.

    Returns the 12-column schema expected by all downstream code.
    """
    return (
        raw_df
        .with_columns(n_bases=pl.col("end") - pl.col("start"))
        .with_columns(cum_coverage=pl.col("n_bases") * pl.col("depth"))
        .with_columns(total_bases=pl.col("n_bases").sum().over("ref"))
        .with_columns(
            mean_coverage=pl.col("cum_coverage").sum().over("ref")
            / pl.col("total_bases").over("ref")
        )
        .with_columns(
            mean_coverage_total=pl.col("cum_coverage").sum()
            / (pl.col(["ref"]).is_first_distinct() * pl.col("total_bases")).sum()
        )
        .with_columns(
            over_zero=pl.when(pl.col("depth") > 0).then(pl.col("n_bases")).otherwise(0)
        )
        .with_columns(
            pct_over_zero=pl.col("over_zero").sum().over("ref")
            / pl.col("total_bases").over("ref")
        )
        .with_columns(
            over_thresh=pl.when(pl.col("depth") > thresh)
            .then(pl.col("n_bases"))
            .otherwise(0)
        )
        .with_columns(
            pct_over_thresh=pl.col("over_thresh").sum().over("ref")
            / pl.col("total_bases").over("ref")
        )
        .with_columns(
            pct_total_over_zero=pl.col("over_zero").sum()
            / (pl.col(["ref"]).is_first_distinct() * pl.col("total_bases")).sum()
        )
        .with_columns(
            pct_total_over_thresh=pl.col("over_thresh").sum()
            / (pl.col(["ref"]).is_first_distinct() * pl.col("total_bases")).sum()
        )
        .select(pl.col("*").exclude("over_zero", "over_thresh", "cum_coverage"))
    )


def print_coverage_info(df, threshold: int) -> None:
    print_blue(f'[SUMMARIZE]: Coverage information for: {df["ref"][0]}')
    print_blue(
        f'[SUMMARIZE]: Percent bases with coverage above 0X: {df["pct_over_zero"][0] * 100: .1f}%'
    )
    print_blue(
        f'[SUMMARIZE]: Percent bases with coverage above {threshold}X: {df["pct_over_thresh"][0] * 100: .1f}%'
    )
    print_blue(f'   [SUMMARIZE]: mean coverage: {df["mean_coverage"][0]: .1f}X')


def print_total_reference_info(df, threshold: int) -> None:
    if df["mean_coverage_total"].shape[0] > 0:
        print_blue(
            f'[SUMMARIZE]: Mean coverage of all basepairs: {df["mean_coverage_total"][0]: .1f}X'
        )
    print_blue(
        f'[SUMMARIZE]: Percent bases with coverage above 0X: {df["pct_total_over_zero"][0] * 100: .1f}%'
    )
    print_blue(
        f'[SUMMARIZE]: Percent bases with coverage above {threshold}X: {df["pct_total_over_thresh"][0] * 100: .1f}%'
    )


def refs_with_most_coverage(df, n=None) -> list:
    if n is None:
        n = df.n_unique("ref")

    if n > 100:
        print_warning(
            f"[WARNING]: Number of reference to plot is {n} – which is too many!"
        )
        print_warning("[WARNING]: Please choose a number < 100 with --number_to_plot")
        sys.exit(1)

    return (
        df.unique(subset="ref")
        .sort("pct_over_thresh", descending=True)[0:n, "ref"]
        .to_list()
    )


def return_ref_for_plotting(df, ref, thresh, rolling_window: int = None):

    if df.filter(pl.col("ref") == ref)["total_bases"][0] < 10_000:
        rolling_window = 1
        modulo = 1
    else:
        if rolling_window is None:
            rolling_window = (
                df.filter(pl.col("ref") == ref)["total_bases"][0] // 100_000
            )
        modulo = df.filter(pl.col("ref") == ref)["total_bases"][0] // 1000

    return (
        df.filter(pl.col("ref") == ref)
        .sort("start")
        .with_columns(pos=pl.int_ranges(start="start", end="end"))
        .explode(pl.col("pos"))
        .with_columns(rolling=pl.col("depth").rolling_mean(window_size=rolling_window))
        .with_columns(
            zero=pl.when(pl.col("rolling") == 0)
            .then(pl.lit(COLOR_ZERO))
            .when(pl.col("rolling") > thresh)
            .then(pl.lit(COLOR_HIGH))
            .otherwise(pl.lit(COLOR_LOW))
        )
        .fill_null(0)
        .with_row_index()
        .filter(pl.col("index") % modulo == 0)
    )


def coverage_plot(
    df: pl.DataFrame,
    x_col: str,
    y_col: str,
    color_col: str,
    thresh: int,
    rolling_window: int,
    sample_name: str,
):
    df_pd = df.to_pandas()

    x = df_pd[x_col].values
    y = df_pd[y_col].values
    colors = df_pd[color_col].values

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, colors=colors, linewidth=3)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.add_collection(lc)
    # ax.autoscale()
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(-1, y.max())

    legend_elements = [
        Line2D([0], [0], color=COLOR_ZERO, lw=1, label="0X"),
        Line2D([0], [0], color=COLOR_LOW, lw=1, label="Over 0X"),
        Line2D([0], [0], color=COLOR_HIGH, lw=1, label=f"Over {thresh}X"),
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    plt.title(
        f"Percent bases with coverage above {thresh}X: {df['pct_over_thresh'][0] * 100: .1f}% | Rolling window: {rolling_window} nt"
    )
    plt.suptitle(f"Ref: {df['ref'][0]} | Sample: {sample_name}")
    plt.close()

    return fig


def coverage_for_value(df, coverage: int):
    number_of_bases = df["total_bases"][0]
    _id = df["ref"][0][:20]

    percent_df = df.with_columns(
        over_t=pl.when(pl.col("depth") > coverage).then(pl.col("n_bases")).otherwise(0)
    )

    percent = percent_df["over_t"].sum() / number_of_bases * 100

    return pl.DataFrame(
        {
            "coverage": [coverage],
            "percent": [percent],
            "id": [_id],
        }
    )


def coverage_for_many_values(df, coverage_values):
    dfs = []
    for coverage in coverage_values:
        coverage_df = coverage_for_value(df, coverage)
        dfs.append(coverage_df)
    return pl.concat(dfs)


def plot_cumulative_coverage_for_all(df, n):
    max_cov = df["depth"].max()
    coverage_values = np.linspace(0, max_cov, 15)
    top_n_refs = refs_with_most_coverage(df, n=n)

    all_coverage = pl.concat(
        [
            coverage_for_many_values(df.filter(pl.col("ref") == ref), coverage_values)
            for ref in top_n_refs
        ]
    )
    grid = sns.FacetGrid(all_coverage, col="id", height=2.5, col_wrap=5)
    grid.map_dataframe(sns.lineplot, x="coverage", y="percent")
    plt.close()
    return grid.fig


def make_dir(outpath: str) -> None:
    outpath = Path(outpath)
    if not outpath.exists():
        outpath.mkdir(parents=True)


def cli():
    args = sys.argv
    valid_subcommand = ["from_bam", "from_reads", "guci"]

    if len(args) < 2:
        print_fail("You must call bam2plot with the following subcommands:")
        print_fail(f"   [1]: '{valid_subcommand[0]}'")
        print_fail(f"   [2]: '{valid_subcommand[1]}'")
        print_fail(f"   [3]: '{valid_subcommand[2]}'")
        sys.exit(1)

    sub_command = args[1]

    if not sub_command in valid_subcommand:
        print_fail(f"The following is not a valid sub command: {sub_command}")
        print_fail(f"Choose between:")
        print_fail(f"   [1]: '{valid_subcommand[0]}'")
        print_fail(f"   [2]: '{valid_subcommand[1]}'")
        print_fail(f"   [3]: '{valid_subcommand[2]}'")
        sys.exit(1)

    if sub_command == "from_bam":
        bam2plot_from_bam()
    elif sub_command == "from_reads":
        bam2plot_from_reads()
    elif sub_command == "guci":
        bam2plot_guci()

    else:
        sys.exit(0)


def bam2plot_guci():
    parser = argparse.ArgumentParser(
        description="Plot GC content of your reference fasta!"
    )
    parser.add_argument("sub_command")
    parser.add_argument("-ref", "--reference", required=True, help="Reference fasta")
    parser.add_argument(
        "-w",
        "--window",
        required=True,
        help="Rolling window size",
        type=int,
    )
    parser.add_argument(
        "-o",
        "--out_folder",
        required=True,
        help="Where to save the plots.",
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

    main_guci(
        ref=args.reference,
        window=args.window,
        out_folder=args.out_folder,
        plot_type=args.plot_type,
    )


def main_guci(
    ref,
    window,
    out_folder,
    plot_type,
):
    print_green(f"[INFO]: Running bam2plot guci!")

    files_not_exists(ref)
    make_dir(out_folder)

    print_green(f"[INFO]: Starting processing {Path(ref).stem}!")

    df = guci(ref, window)
    title = f"GC content of: {Path(ref).stem}. Rolling window: {window}"
    plot = plot_gc(df, title)

    out_name = f"{out_folder}/gc_{Path(ref).stem}.{plot_type}"
    plot.savefig(out_name)

    print_green(f"[INFO]: Guci plot done!")
    print_green(f"[INFO]: Plot location: {Path(out_folder).resolve()}")

    sys.exit(0)


def bam2plot_from_reads():
    parser = argparse.ArgumentParser(
        description="Align your reads and plot the coverage!"
    )
    parser.add_argument("sub_command")
    parser.add_argument(
        "-r1", "--read_1", required=True, help="Fastq file 1 (Required)"
    )
    parser.add_argument(
        "-r2", "--read_2", required=False, default=None, help="Fastq file 2 (Optional)"
    )
    parser.add_argument("-ref", "--reference", required=True, help="Reference fasta")
    parser.add_argument(
        "-gc",
        "--guci",
        required=False,
        default=False,
        help="Plot GC content?",
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "-o",
        "--out_folder",
        required=True,
        help="Where to save the plots.",
    )
    parser.add_argument(
        "-r",
        "--rolling_window",
        required=False,
        default=50,
        help="Rolling window size",
        type=int,
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

    main_from_reads(
        read_1=args.read_1,
        read_2=args.read_2,
        ref=args.reference,
        window=args.rolling_window,
        gc=args.guci,
        out_folder=args.out_folder,
        plot_type=args.plot_type,
    )


def map_fastq_to_ref_long_read(fastq, ref, preset="map-ont") -> pl.DataFrame:
    a = mp.Aligner(ref, preset=preset)

    # see if the ref genome is the adequate one
    test = 0
    ref_len = None
    for _, seq, _ in mp.fastx_read(fastq):
        for hit in a.map(seq):
            ref_len = hit.ctg_len

        if ref_len is not None:
            break

        test += 1
        if test > 1000:
            print_fail("[ERROR]: No alignment to this reference")
            return None

    df = pl.DataFrame({"pos": np.arange(1, ref_len + 1), "depth": 0})

    for _, seq, _ in mp.fastx_read(fastq):
        for hit in a.map(seq):
            df = df.with_columns(
                depth=pl.when(pl.col("pos").is_between(hit.r_st, hit.r_en))
                .then(pl.col("depth") + 1)
                .otherwise(pl.col("depth"))
            )
    return df


def PE_reads_to_df(read_1, read_2) -> pl.DataFrame:
    def parse_PE_read(fastx, read_name, comment_name, seq_name):

        reads = pyfastx.Fastx(fastx, comment=True)
        names = []
        seqs = []
        comments = []

        for name, seq, qual, comment in reads:
            names.append(name)
            seqs.append(seq)
            comments.append(comment)

        return pl.DataFrame({read_name: names, comment_name: comments, seq_name: seqs})

    read_1 = parse_PE_read(
        read_1, read_name="name_1", comment_name="comment_1", seq_name="seq_1"
    )
    read_2 = parse_PE_read(
        read_2, read_name="name_2", comment_name="comment_2", seq_name="seq_2"
    )

    try:
        return (
            pl.concat([read_1, read_2], how="horizontal")
            .sort(pl.col("comment_1"), pl.col("comment_2"))
            .with_columns(
                comment_1=pl.col("comment_1").str.replace("/.*", ""),
                comment_2=pl.col("comment_2").str.replace("/.*", ""),
            )
            .filter(pl.col("comment_1") == pl.col("comment_2"))
            .select(pl.col("seq_1"), pl.col("seq_2"))
        )

    except Exception as e:
        print_fail(f"[ERROR]: fastq files do not match up! ({e})")
        sys.exit(1)


def map_fastq_to_ref_PE_read(fastq_1, fastq_2, ref, preset="sr") -> pl.DataFrame:
    a = mp.Aligner(ref, preset=preset)

    # see if the ref genome is the adequate one
    test = 0
    ref_len = None
    PE_df = PE_reads_to_df(fastq_1, fastq_2)
    for reads in PE_df.iter_rows(named=True):
        for hit in a.map(reads["seq_1"], reads["seq_2"]):
            ref_len = hit.ctg_len

        if ref_len is not None:
            break

        test += 1
        if test > 1000:
            print_fail(f"No alignment to this reference: {ref}")
            sys.exit(1)

    df = pl.DataFrame({"pos": np.arange(1, ref_len + 1), "depth": 0})

    counter = 0
    print_green(f"Processing read: 1 of {PE_df.shape[0]}")
    for reads in PE_df.iter_rows(named=True):
        counter += 1
        if counter % 50_000 == 0:
            print_green(f"Processing read: {counter} of {PE_df.shape[0]}")
        for hit in a.map(reads["seq_1"], reads["seq_2"]):
            df = df.with_columns(
                depth=pl.when(pl.col("pos").is_between(hit.r_st, hit.r_en))
                .then(pl.col("depth") + 1)
                .otherwise(pl.col("depth"))
            )
    return df


def plot_from_reads(
    df,
    sample_name,
    ref_name,
    window,
):
    coverage_plot = plt.figure(figsize=(15, 8))
    ax = sns.lineplot(data=df, x="pos", y="rolling")
    plt.title(f" Sample: {sample_name}| Ref: {ref_name} | Rolling window: {window}")
    ax.set(xlabel="Position", ylabel="Depth")
    zero = plt.axhline(y=0, color="red")
    zero.set_label("Zero")
    plt.close()

    return coverage_plot


def ref_to_seq_df(fastx_file: str) -> pl.DataFrame:
    fastx = pyfastx.Fastx(fastx_file)
    reads = list(zip(*[[x[0], x[1]] for x in fastx]))

    df = (
        pl.DataFrame(
            {
                "name": reads[0],
                "sequence": reads[1],
            }
        )
        .select(pl.col("sequence").str.split("").explode(), pl.col("name"))
        .filter(pl.col("sequence") != "")
        .with_row_index(name="position", offset=1)
    )

    return df


def add_gc(df) -> pl.DataFrame:
    return df.with_columns(
        gc=pl.when(pl.col("sequence").str.contains("G|C")).then(1).otherwise(0)
    )


def add_rolling_mean(df, window):
    return df.with_columns(rolling_gc=pl.col("gc").rolling_mean(window_size=window))


def guci(fastx_file, window):
    df = ref_to_seq_df(fastx_file)
    df = add_gc(df)
    df = add_rolling_mean(df, window).select(pl.col("position"), pl.col("rolling_gc"))
    return df


def plot_gc(df, title):
    fig = plt.figure(figsize=(50, 15))
    plt.plot(df["position"], df["rolling_gc"])
    plt.gca().set_yticklabels([f"{x:.0%}" for x in plt.gca().get_yticks()])
    plt.gca().xaxis.set_major_formatter(
        matplotlib.ticker.StrMethodFormatter("{x:,.0f}")
    )
    plt.ylabel("% GC content")
    plt.xlabel("Position")
    plt.title(title)
    plt.close()
    return fig


def files_not_exists(*files):
    for file in files:
        if file is None:
            continue
        if not Path(file).exists():
            print_fail(f"[ERROR]: The file {file} does not exist")
            sys.exit(1)


def main_from_reads(
    read_1,
    read_2,
    ref,
    window,
    gc,
    out_folder,
    plot_type,
) -> None:
    print_green(f"[INFO]: Running bam2plot from_reads!")

    files_not_exists(read_1, read_2, ref)
    make_dir(out_folder)
    sample_name = Path(read_1).stem
    ref_name = Path(ref).stem

    if read_2 is None:
        print_green(f"[INFO]: Running bam2plot on reads from: {read_1}!")
        df = map_fastq_to_ref_long_read(read_1, ref)
    else:
        print_green(f"[INFO]: Running bam2plot on reads from: {read_1} and {read_2}!")
        df = map_fastq_to_ref_PE_read(read_1, read_2, ref)

    if gc:
        df = df.with_columns(
            rolling=pl.col("depth").rolling_mean(window_size=window)
        ).with_columns(rolling=pl.col("rolling") / pl.col("rolling").max())
    else:
        df = df.with_columns(rolling=pl.col("depth").rolling_mean(window_size=window))

    coverage_plot = plot_from_reads(
        df,
        sample_name=sample_name,
        ref_name=ref_name,
        window=window,
    )

    if gc:
        guci_df = guci(ref, window)
        sns.lineplot(
            x="position",
            y="rolling_gc",
            data=guci_df,
            color="orange",
            ax=coverage_plot.axes[0],
        )

    save_plot_coverage(coverage_plot, out_folder, sample_name, ref_name, plot_type)

    print_green(f"[INFO]: Coverage plot done!")
    print_green(f"[INFO]: Plot location: {Path(out_folder).resolve()}")

    sys.exit(0)


def check_range(value):
    ivalue = int(value)
    if ivalue < 1 or ivalue > 100:
        print_fail(
            f"   [ERROR]: {value} is an invalid number of references. Choose a number between 1 and 100."
        )
        sys.exit(1)
    return ivalue


def bam2plot_from_bam():
    parser = argparse.ArgumentParser(description="Plot your bam files!")
    parser.add_argument("sub_command")
    parser.add_argument("-b", "--bam", required=True, help="bam file")
    parser.add_argument(
        "-o",
        "--outpath",
        required=True,
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
        default=10,
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
        "-c",
        "--cum_plot",
        required=False,
        default=False,
        help="Generate cumulative plots of all chromosomes",
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
    parser.add_argument(
        "-n",
        "--number_of_refs",
        required=False,
        default=10,
        type=check_range,  # cap
        help="How many references (chromosomes) to plot",
    )

    args = parser.parse_args()
    command = "\nbam2plot \n" + "".join(f"{k}: {v}\n" for k, v in vars(args).items())
    print_green(command)

    main_from_bam(
        bam=args.bam,
        outpath=args.outpath,
        whitelist=args.whitelist,
        rolling_window=args.rolling_window,
        threshold=args.threshold,
        index=args.index,
        sort_and_index=args.sort_and_index,
        zoom=args.zoom,
        cum_plot=args.cum_plot,
        plot_type=args.plot_type,
        number_of_refs=args.number_of_refs,
    )


def if_sort_and_index(sort_and_index, index, bam) -> str:
    """Optionally sort/index the BAM and return the path to use for coverage."""
    if sort_and_index:
        print_green("[INFO]: Sorting bam file")
        sort_bam(bam, new_name=SORTED_TEMP)
        print_green("[INFO]: Indexing bam file")
        index_bam(SORTED_TEMP, new_name=SORTED_TEMP_INDEX)
        return SORTED_TEMP

    if index:
        print_green("[INFO]: Indexing bam file")
        index_bam(bam, new_name=f"{bam}.bai")

    return bam


def process_dataframe(bam, threshold):
    try:
        print_green("[INFO]: Processing dataframe")
        raw_df = bam_to_raw_df(bam)
        df = enrich_coverage_df(raw_df, threshold)
        return df
    except Exception as e:
        print_fail(f"[ERROR]: Could not process dataframe ({e})")
        print_warning(
            "[WARNING]: Is the file properly prepared? If not, consider running 'bam2plot <file.bam> -s' or 'bam2plot <file.bam> -i'"
        )
        sys.exit(1)


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

    if plot_type == "png":
        cum_plot.savefig(f"{cum_plot_name}.png")

    if plot_type == "svg":
        cum_plot.savefig(f"{cum_plot_name}.svg")

    if plot_type == "both":
        cum_plot.savefig(f"{cum_plot_name}.png")
        cum_plot.savefig(f"{cum_plot_name}.svg")

    print_green(f"[INFO]: Cumulative plot generated!")


def _fig_to_base64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def generate_html_report(
    sample_name: str,
    bam: str,
    threshold: int,
    rolling_window: int,
    df: pl.DataFrame,
    top_refs: list,
    coverage_figures: list,
    cum_fig,
    outpath: str,
) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Global summary stats
    first_row = df.row(0, named=True)
    global_mean_cov = first_row["mean_coverage_total"]
    global_pct_zero = first_row["pct_total_over_zero"] * 100
    global_pct_thresh = first_row["pct_total_over_thresh"] * 100

    # Per-reference stats
    per_ref_rows = []
    for ref in top_refs:
        ref_df = df.filter(pl.col("ref") == ref)
        row = ref_df.row(0, named=True)
        per_ref_rows.append(
            f"""<tr>
            <td>{ref}</td>
            <td>{row['total_bases']:,}</td>
            <td>{row['mean_coverage']:.1f}X</td>
            <td>{row['pct_over_zero'] * 100:.1f}%</td>
            <td>{row['pct_over_thresh'] * 100:.1f}%</td>
            </tr>"""
        )
    per_ref_html = "\n".join(per_ref_rows)

    # Coverage plot sections
    plot_sections = []
    for ref_name, fig in coverage_figures:
        b64 = _fig_to_base64(fig)
        plot_sections.append(
            f"""<div class="plot-section">
            <h3>{ref_name}</h3>
            <img src="data:image/png;base64,{b64}" alt="Coverage plot for {ref_name}">
            </div>"""
        )
    plots_html = "\n".join(plot_sections)

    # Cumulative plot
    cum_html = ""
    if cum_fig is not None:
        cum_b64 = _fig_to_base64(cum_fig)
        cum_html = f"""<div class="plot-section">
        <h2>Cumulative Coverage</h2>
        <img src="data:image/png;base64,{cum_b64}" alt="Cumulative coverage plot">
        </div>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>bam2plot Report — {sample_name}</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 0; padding: 20px; background: #f5f5f5; color: #333; }}
.container {{ max-width: 1200px; margin: 0 auto; background: #fff; padding: 30px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
h1 {{ color: #2c3e50; border-bottom: 2px solid #4A90E2; padding-bottom: 10px; }}
h2 {{ color: #34495e; margin-top: 30px; }}
h3 {{ color: #4A90E2; }}
table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
th, td {{ padding: 10px 14px; text-align: left; border: 1px solid #ddd; }}
th {{ background: #4A90E2; color: #fff; }}
tr:nth-child(even) {{ background: #f9f9f9; }}
.meta {{ color: #777; font-size: 0.9em; margin-bottom: 20px; }}
.plot-section {{ margin: 20px 0; }}
.plot-section img {{ max-width: 100%; height: auto; }}
</style>
</head>
<body>
<div class="container">
<h1>bam2plot Report</h1>
<div class="meta">
<p><strong>Sample:</strong> {sample_name} &nbsp;|&nbsp; <strong>BAM:</strong> {bam}</p>
<p><strong>Threshold:</strong> {threshold}X &nbsp;|&nbsp; <strong>Rolling window:</strong> {rolling_window} nt &nbsp;|&nbsp; <strong>Generated:</strong> {timestamp}</p>
</div>

<h2>Global Summary</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Mean coverage</td><td>{global_mean_cov:.1f}X</td></tr>
<tr><td>Bases with coverage &gt; 0X</td><td>{global_pct_zero:.1f}%</td></tr>
<tr><td>Bases with coverage &gt; {threshold}X</td><td>{global_pct_thresh:.1f}%</td></tr>
</table>

<h2>Per-Reference Statistics</h2>
<table>
<tr><th>Reference</th><th>Total bases</th><th>Mean coverage</th><th>&gt; 0X</th><th>&gt; {threshold}X</th></tr>
{per_ref_html}
</table>

<h2>Coverage Plots</h2>
{plots_html}

{cum_html}
</div>
</body>
</html>"""

    report_path = Path(outpath) / f"{sample_name}_report.html"
    report_path.write_text(html)
    print_green(f"[INFO]: HTML report saved to {report_path.resolve()}")


def main_from_bam(
    bam,
    outpath,
    whitelist,
    rolling_window,
    threshold,
    index,
    sort_and_index,
    zoom,
    cum_plot,
    plot_type,
    number_of_refs,
) -> None:
    print_green(f"[INFO]: Running bam2plot from_bam!")

    if zoom:
        try:
            parts = zoom.split(" ")
            start = int(parts[0])
            end = int(parts[1])
        except (IndexError, ValueError):
            print_fail("[ERROR]: Invalid zoom format. Use: -z='100 2000'")
            sys.exit(1)
        if start >= end:
            print_fail("[ERROR]: Start value of zoom must be lower than end value.")
            sys.exit(1)

    if not Path(bam).exists():
        print_fail(f"[ERROR]: The file {bam} does not exist")
        sys.exit(1)

    make_dir(outpath)
    sample_name = Path(bam).stem

    effective_bam = if_sort_and_index(sort_and_index, index, bam)
    try:
        df = process_dataframe(effective_bam, threshold)
    finally:
        if sort_and_index:
            Path(SORTED_TEMP).unlink(missing_ok=True)
            Path(SORTED_TEMP_INDEX).unlink(missing_ok=True)
        if index:
            Path(f"{bam}.bai").unlink(missing_ok=True)

    print_total_reference_info(df, threshold)

    if whitelist:
        whitelist = [whitelist] if type(whitelist) == str else whitelist
        print_green(
            f"[INFO]: Only looking for references in the whitelist: {whitelist}"
        )
        df = df.filter(pl.col("ref").is_in(whitelist))

    if number_of_refs == 0:
        print_fail("[ERROR]: No reference to plot against!")
        sys.exit(1)

    plot_text = "plot" if number_of_refs == 1 else "plots"
    print_green(f"[INFO]: Generating {number_of_refs} {plot_text}:")

    top_n_refs = refs_with_most_coverage(df, n=number_of_refs)
    coverage_figures = []

    for i, reference in enumerate(top_n_refs):
        df_to_plot = return_ref_for_plotting(df, reference, threshold, rolling_window)

        if zoom:
            df_to_plot = df_to_plot.filter(pl.col("pos").is_between(start, end))
            if df_to_plot.shape[0] == 0:
                print_warning("[WARNING]: No positions to plot after zoom")
                continue

        if df_to_plot.shape[0] == 0:
            print_warning("[WARNING]: No positions to plot")
            continue

        print_coverage_info(df_to_plot, threshold)

        plot = coverage_plot(
            df_to_plot,
            "pos",
            "rolling",
            "zero",
            threshold,
            rolling_window,
            sample_name,
        )
        coverage_figures.append((reference, plot))
        plot_name = f"{i}_{sample_name}"
        save_plot_coverage(plot, outpath, sample_name, reference, plot_type)

    print_green("[INFO]: Coverage plots done!")

    cum_fig = None
    if cum_plot:
        print_green("[INFO]: Generating cumulative coverage plots for each reference")
        cum_fig = plot_cumulative_coverage_for_all(df, n=number_of_refs)
        save_plot_cum(cum_fig, outpath, bam, plot_type)

    generate_html_report(
        sample_name=sample_name,
        bam=bam,
        threshold=threshold,
        rolling_window=rolling_window,
        df=df,
        top_refs=top_n_refs,
        coverage_figures=coverage_figures,
        cum_fig=cum_fig,
        outpath=outpath,
    )

    print_green(f"[INFO]: Plots location: {Path(outpath).resolve()}")
    sys.exit(0)


if __name__ == "__main__":
    cli()
