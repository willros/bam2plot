from pathlib import Path
import subprocess
import fire
import pysam
import pandas as pd
from io import StringIO
import _io
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import os
import pandas

# for windows users
matplotlib.use('Agg')


SORTED_TEMP = "TEMP112233.sorted.bam"
SORTED_TEMP_INDEX = f"{SORTED_TEMP}.bai"

def sort_bam(bam: str, new_name: str) -> None:
    pysam.sort("-o", new_name, bam)
    
def index_bam(bam: str, new_name: str) -> None:
    pysam.index(bam, new_name)
    
    
    
    
def run_perbase(bam: str) -> _io.StringIO:
    return StringIO(
        subprocess.check_output(
            f"perbase only-depth {bam}", 
            shell=True, 
            stderr=subprocess.DEVNULL
        )
        .decode()
    )

def perbase_to_df(perbase: _io.StringIO) -> pd.DataFrame:
    return (
        pd.read_csv(
            perbase, 
            sep="\t", 
        )
        .rename(columns={"REF": "id", "POS": "Position", "DEPTH": "coverage"})
    )



def interpolate_df_ranges(df, rolling_window: int):
    stop = df.END.max()
    df_to_interpolate = (
        pd.DataFrame()
        .assign(Position=range(1, stop + 1))
    )
    
    return (
        df
        .merge(df_to_interpolate, on="Position", how="right")
        .drop(columns="END")
        .assign(coverage=lambda x: x.coverage.ffill())
        .assign(id=lambda x: x.id.ffill())
        .assign(Depth=lambda x: x.coverage.rolling(rolling_window).mean())
    )


def wrangle_perbase_df(df, rolling_window: int):
    return pd.concat(
        [
            interpolate_df_ranges(df.loc[lambda x: x.id == ref], rolling_window) 
            for ref in df.id.unique()
        ],
        ignore_index=True
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
    index: bool = False,
    sort_and_index: bool = False
) -> None:
    sample_name = Path(bam).stem
        
    if outpath == "":
        outpath = "bam2plot"
        make_dir(outpath)
    else:
        make_dir(outpath)
        
    out_file = f"{outpath}/{sample_name}_bam2plot"
    
        
    if sort_and_index:
        print("Sorting bam file")
        sort_bam(bam, new_name=SORTED_TEMP)
        print("Indexing bam file")
        index_bam(SORTED_TEMP, new_name=SORTED_TEMP_INDEX)
        perbase = run_perbase(SORTED_TEMP)
        os.remove(SORTED_TEMP)
        os.remove(SORTED_TEMP_INDEX)
    else:
        if index:
            print("Indexing bam file")
            index_name = f"{bam}.bai"
            index_bam(bam, new_name=index_name)
        perbase = run_perbase(bam)
        if index:
            os.remove(index_name)

    try:
        print("Processing dataframe")
        df = perbase_to_df(perbase)
        wrangled_df = wrangle_perbase_df(df, rolling_window)
    except pandas.errors.EmptyDataError as e:
        print("Error while processing bam")
        if not sort_and_index:
            print("Try sorting and indexing the file (-s True)")
            exit(1)
        if not index:
            print("Try indexing the file (-i True)")
            exit(1)

    plot_number = wrangled_df.id.nunique()
    if plot_number == 0:
        print("No reference to plot against!")
        exit(1)

    print(f"Generating {plot_number} plots:")
    for reference in wrangled_df.id.unique():
        mpileup_df = wrangled_df.loc[lambda x: x.id == reference]
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
