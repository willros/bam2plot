from pathlib import Path
import fire
import pysam
import pandas as pd
from io import StringIO
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import os
import polars as pl

# for windows users
matplotlib.use('Agg')


SORTED_TEMP = "TEMP112233.sorted"
MPILEUP_TEMP = "TEMP112233.mpileup"


def sort_bam(bam: str, name: str) -> None:
    pysam.sort("-o", name, bam)

    
def run_mpileup(sorted_bam: str) -> None:
    file = open(MPILEUP_TEMP, "w")
    pysam.mpileup("-a", sorted_bam, save_stdout="TEMP112233.mpileup", catch_stdout=False)
    file.close()


def create_mpileup_df(
    mpileup: str,
    rolling_window: int,
) -> pd.DataFrame:
    return (
        pl.read_csv(
            mpileup, 
            sep="\t", 
            has_header=False, 
            new_columns=["id", "Position", "base", "coverage", "x", "y"]
        )
        .drop(["base", "x", "y"])
        .with_columns(
            pl.col("coverage").rolling_mean(window_size=10).alias("Depth")
        )
        .to_pandas()
    )


def plot_coverage(
    mpileup_df: pd.DataFrame,
    sample_name: str,
    threshold: int,
    rolling_window,
) -> matplotlib.figure.Figure:
    mean_coverage = mpileup_df.coverage.mean()
    coverage = (
        sum(1 if x > threshold else 0 for x in mpileup_df.coverage)
        / mpileup_df.shape[0]
        * 100
    )

    sns.set_theme()
    coverage_plot = plt.figure(figsize=(15, 8))
    sns.lineplot(data=mpileup_df, x="Position", y="Depth")
    zero = plt.axhline(y=0, color="red")
    zero.set_label("Zero")
    mean = plt.axhline(y=mean_coverage, color="green")
    mean.set_label(f"Mean coverage: {mean_coverage: .1f}X")
    plt.legend(loc="upper right")
    plt.title(
        f"Percent bases with coverage above {threshold}X: {coverage: .1f}% | Rolling window: {rolling_window} nt"
    )
    plt.suptitle(f"Ref: {mpileup_df.iloc[0].id} | Sample: {sample_name}")
    plt.close()
    return coverage_plot


def make_dir(outpath: str) -> None:
    outpath = Path(outpath)
    if not outpath.exists():
        outpath.mkdir(parents=True)



def cli(
    bam: str,
    outpath: str = "",
    rolling_window: int = 100,
    threshold: int = 3,
) -> None:
    sample_name = Path(bam).stem
        
    if outpath == "":
        outpath = "bam2plot"
        make_dir(outpath)
    else:
        make_dir(outpath)
        
    out_file = f"{outpath}/{sample_name}_bam2plot"

    print("Sorting bam file")
    sort_bam(bam, name=SORTED_TEMP)
    print("Creating mpileup")
    run_mpileup(SORTED_TEMP)
    os.remove(SORTED_TEMP)
    print("Parsing bamfile")
    df = create_mpileup_df(MPILEUP_TEMP, rolling_window=rolling_window)
    os.remove(MPILEUP_TEMP)

    plot_number = df.id.nunique()
    if plot_number == 0:
        print("No reference to plot against!")
        exit(1)

    print(f"Generating {plot_number} plots:")
    for reference in df.id.unique():
        mpileup_df = df.loc[lambda x: x.id == reference]
        plot = plot_coverage(
            mpileup_df, sample_name, threshold=threshold, rolling_window=rolling_window
        )

        name = f"{out_file}_{reference}"
        plot.savefig(f"{name}.svg")
        plot.savefig(f"{name}.png")
        print(f"Plot for {reference} generated")

    print("Plots done!")
    print(f"Plots location: {Path(outpath).resolve()}")
    exit(0)


def run() -> None:
    fire.Fire(cli)
