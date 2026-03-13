"""Public Python API for bam2plot.

The package root exposes a small set of reusable functions intended for
notebooks and pipelines. CLI-only entrypoints remain in ``bam2plot.main``.
"""

from importlib import import_module
from typing import TYPE_CHECKING

__author__ = "William Rosenbaum"
__version__ = "0.4.0"

__all__ = [
    "__author__",
    "__version__",
    "bam_to_raw_df",
    "enrich_coverage_df",
    "extract_insert_sizes",
    "weighted_median_rle",
    "compute_lorenz_data_rle",
    "compute_gini_rle",
    "write_summary_outputs",
    "generate_html_report",
    "guci",
]

_PUBLIC_API = {
    "bam_to_raw_df": "bam2plot.main",
    "enrich_coverage_df": "bam2plot.main",
    "extract_insert_sizes": "bam2plot.main",
    "weighted_median_rle": "bam2plot.main",
    "compute_lorenz_data_rle": "bam2plot.main",
    "compute_gini_rle": "bam2plot.main",
    "write_summary_outputs": "bam2plot.main",
    "generate_html_report": "bam2plot.main",
    "guci": "bam2plot.main",
}


def __getattr__(name: str):
    if name in _PUBLIC_API:
        module = import_module(_PUBLIC_API[name])
        return getattr(module, name)
    raise AttributeError(f"module 'bam2plot' has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(__all__)


if TYPE_CHECKING:
    from bam2plot.main import (
        bam_to_raw_df,
        compute_gini_rle,
        compute_lorenz_data_rle,
        enrich_coverage_df,
        extract_insert_sizes,
        generate_html_report,
        guci,
        weighted_median_rle,
        write_summary_outputs,
    )
