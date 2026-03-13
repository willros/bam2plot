import bam2plot
import bam2plot.main as main_mod


def test_public_api_exports_expected_names():
    expected = {
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
    }

    assert expected.issubset(set(bam2plot.__all__))
    assert expected.issubset(set(dir(bam2plot)))


def test_public_api_reexports_selected_helpers_from_main():
    assert bam2plot.bam_to_raw_df is main_mod.bam_to_raw_df
    assert bam2plot.enrich_coverage_df is main_mod.enrich_coverage_df
    assert bam2plot.extract_insert_sizes is main_mod.extract_insert_sizes
    assert bam2plot.weighted_median_rle is main_mod.weighted_median_rle
    assert bam2plot.compute_lorenz_data_rle is main_mod.compute_lorenz_data_rle
    assert bam2plot.compute_gini_rle is main_mod.compute_gini_rle
    assert bam2plot.write_summary_outputs is main_mod.write_summary_outputs
    assert bam2plot.generate_html_report is main_mod.generate_html_report
    assert bam2plot.guci is main_mod.guci
