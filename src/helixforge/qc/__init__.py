"""Quality control for HelixForge.

This module provides tools for quality control of gene predictions:

- QC flags for various issues (confidence, splice, homology, structure)
- Aggregation of results from multiple modules
- Filtering based on QC criteria with preset profiles
- Report generation (HTML with Chart.js visualizations)

Example:
    >>> from helixforge.qc import QCAggregator, GeneFilter, QCReportGenerator
    >>> from helixforge.qc import Flags, GeneQC
    >>>
    >>> # Aggregate results from modules
    >>> aggregator = QCAggregator()
    >>> gene_qcs = aggregator.aggregate(
    ...     confidence_results=confidence_data,
    ...     splice_results=splice_data,
    ...     homology_results=homology_data,
    ... )
    >>>
    >>> # Filter genes
    >>> result = GeneFilter.publication_ready().apply(gene_qcs)
    >>> print(f"Publication-ready: {result.pass_count}")
    >>>
    >>> # Generate report
    >>> generator = QCReportGenerator()
    >>> generator.generate(gene_qcs, "qc_report.html")
"""

from helixforge.qc.flags import (
    FlagCategory,
    FlagSeverity,
    Flags,
    GeneQC,
    QCFlag,
    filter_by_category,
    filter_by_severity,
    max_severity,
    summarize_flags,
)
from helixforge.qc.aggregate import (
    QCAggregator,
    QCAggregatorConfig,
    export_qc_tsv,
    get_genes_by_tier,
    get_genes_with_flag,
    summarize_qc_results,
)
from helixforge.qc.filters import (
    FilterCriteria,
    FilterResult,
    GeneFilter,
    calculate_filter_statistics,
    split_by_tier,
    summarize_filter_results,
    write_filtered_gff,
    write_tiered_gene_lists,
)
from helixforge.qc.report import (
    QCReportGenerator,
    ReportConfig,
    ReportData,
    generate_json_report,
    generate_summary_report,
)

__all__ = [
    # Flags
    "QCFlag",
    "Flags",
    "GeneQC",
    "FlagSeverity",
    "FlagCategory",
    "summarize_flags",
    "max_severity",
    "filter_by_severity",
    "filter_by_category",
    # Aggregation
    "QCAggregator",
    "QCAggregatorConfig",
    "summarize_qc_results",
    "get_genes_by_tier",
    "get_genes_with_flag",
    "export_qc_tsv",
    # Filtering
    "GeneFilter",
    "FilterCriteria",
    "FilterResult",
    "split_by_tier",
    "write_tiered_gene_lists",
    "write_filtered_gff",
    "calculate_filter_statistics",
    "summarize_filter_results",
    # Report
    "QCReportGenerator",
    "ReportConfig",
    "ReportData",
    "generate_summary_report",
    "generate_json_report",
]
