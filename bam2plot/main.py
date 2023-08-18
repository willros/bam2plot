from pathlib import Path
import fire
import pysam
import pandas as pd
from io import StringIO
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import os
import _io


SORTED_TEMP = "temp112233.sorted"


def sort_bam(bam: str, name: str) -> None:
    pysam.sort("-o", name, bam)

def run_mpileup(sorted_bam: str) -> _io.StringIO:
    return StringIO(pysam.mpileup("-a", sorted_bam))
    
def create_mpileup_df(
    sorted_bam: str,
    rolling_window: int,
) -> pd.DataFrame:
    return (
        pd.read_csv(
            run_mpileup(sorted_bam),
            sep="\t",
            header=None,
            names=["id", "Position", "base", "coverage", "x", "y"],
        )
        .assign(Depth=lambda x: x.coverage.rolling(rolling_window).mean())
        .drop(columns=["base", "x", "y"])
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
    plt.title(f"Percent bases with coverage above {threshold}X: {coverage: .1f}% | Rolling window: {rolling_window} nt")
    plt.suptitle(f"Ref: {mpileup_df.iloc[0].id} | Sample: {sample_name}")
    plt.close()
    return coverage_plot


def cli(
    bam: str,
    sample_name: str = "",
    outpath: str = "",
    rolling_window: int = 100,
    threshold: int = 3,
) -> None:
    if sample_name == "":
        sample_name = Path(bam).stem
    if outpath == "":
        outpath = f"{sample_name}_bam2plot"
        
    print("Sorting bam file")
    sort_bam(bam, name=SORTED_TEMP)
    print("Parsing bamfile")
    df = create_mpileup_df(SORTED_TEMP, rolling_window=rolling_window)
    os.remove(SORTED_TEMP)
    
    plot_number = df.id.nunique()
    if plot_number == 0:
        print("No reference to plot against!")
        exit(1)
        
    print(f"Generating {plot_number} plots:")
    for reference in df.id.unique():
        mpileup_df = df.loc[lambda x: x.id == reference]
        plot = plot_coverage(mpileup_df, sample_name, threshold=threshold, rolling_window=rolling_window)
        
        name = f"{outpath}_{reference}"
        plot.savefig(f"{name}.svg")
        plot.savefig(f"{name}.png")
        print(f"Plot for {reference} generated")
    
        
    print("Plots done!")
    print(f"Plots location: {Path(outpath).resolve()}")
    exit(0)


def run() -> None:
    fire.Fire(
        cli
    )
