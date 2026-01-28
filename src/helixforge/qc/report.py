"""QC report generation.

This module provides tools for generating comprehensive QC reports
in HTML and other formats.

Example:
    >>> from helixforge.qc import QCReportGenerator
    >>> generator = QCReportGenerator()
    >>> generator.generate(
    ...     gene_qcs=gene_qcs,
    ...     output_path="qc_report.html",
    ...     title="Genome Annotation QC Report",
    ... )
"""

import json
from collections.abc import Iterable
from datetime import datetime
from pathlib import Path
from typing import Any

import attrs

from helixforge.qc.aggregate import summarize_qc_results
from helixforge.qc.flags import (
    FlagCategory,
    FlagSeverity,
    Flags,
    GeneQC,
)


# =============================================================================
# Report Configuration
# =============================================================================


@attrs.define
class ReportConfig:
    """Configuration for QC report generation.

    Attributes:
        title: Report title.
        description: Report description.
        include_charts: Whether to include Chart.js visualizations.
        include_gene_table: Whether to include detailed gene table.
        max_genes_in_table: Maximum genes to show in table.
        include_summary: Whether to include summary statistics.
        include_flag_details: Whether to include flag details section.
        custom_css: Optional custom CSS to include.
    """

    title: str = "HelixForge QC Report"
    description: str = ""
    include_charts: bool = True
    include_gene_table: bool = True
    max_genes_in_table: int = 1000
    include_summary: bool = True
    include_flag_details: bool = True
    custom_css: str = ""


# =============================================================================
# Report Data
# =============================================================================


@attrs.define
class ReportData:
    """Data prepared for report rendering.

    Attributes:
        title: Report title.
        description: Report description.
        generated_at: Report generation timestamp.
        summary: Summary statistics dict.
        tier_data: Tier distribution data for charts.
        severity_data: Severity distribution data for charts.
        category_data: Category distribution data for charts.
        flag_counts_data: Per-flag count data for charts.
        aed_distribution_data: AED score distribution data for histogram.
        confidence_distribution_data: Confidence score distribution data.
        flag_details: List of flag detail dicts.
        gene_table: List of gene data dicts for table.
        config: Report configuration.
    """

    title: str
    description: str
    generated_at: str
    summary: dict[str, Any]
    tier_data: dict[str, Any]
    severity_data: dict[str, Any]
    category_data: dict[str, Any]
    flag_counts_data: dict[str, Any]
    aed_distribution_data: dict[str, Any]
    confidence_distribution_data: dict[str, Any]
    flag_details: list[dict[str, Any]]
    gene_table: list[dict[str, Any]]
    config: ReportConfig


# =============================================================================
# Report Generator
# =============================================================================


@attrs.define
class QCReportGenerator:
    """Generate HTML QC reports with visualizations.

    Uses Jinja2 templates and Chart.js for interactive visualizations.

    Example:
        >>> generator = QCReportGenerator()
        >>> generator.generate(gene_qcs, "report.html")
    """

    config: ReportConfig = attrs.Factory(ReportConfig)

    def generate(
        self,
        gene_qcs: dict[str, GeneQC] | Iterable[GeneQC],
        output_path: Path | str,
        title: str | None = None,
        description: str | None = None,
    ) -> Path:
        """Generate HTML QC report.

        Args:
            gene_qcs: Dict or iterable of GeneQC objects.
            output_path: Path to write HTML report.
            title: Override report title.
            description: Override report description.

        Returns:
            Path to generated report.
        """
        # Convert to dict if needed
        if isinstance(gene_qcs, dict):
            qc_dict = gene_qcs
        else:
            qc_dict = {qc.gene_id: qc for qc in gene_qcs}

        # Prepare report data
        report_data = self._prepare_data(qc_dict, title, description)

        # Render HTML
        html_content = self._render_html(report_data)

        # Write output
        output_path = Path(output_path)
        output_path.write_text(html_content)

        return output_path

    def _prepare_data(
        self,
        gene_qcs: dict[str, GeneQC],
        title: str | None,
        description: str | None,
    ) -> ReportData:
        """Prepare data for report rendering.

        Args:
            gene_qcs: Dict of GeneQC objects.
            title: Override title.
            description: Override description.

        Returns:
            ReportData object.
        """
        # Get summary statistics
        summary = summarize_qc_results(gene_qcs)

        # Prepare chart data
        tier_data = self._prepare_tier_data(summary)
        severity_data = self._prepare_severity_data(summary)
        category_data = self._prepare_category_data(summary)
        flag_counts_data = self._prepare_flag_counts_data(summary)
        aed_distribution_data = self._prepare_aed_distribution(gene_qcs)
        confidence_distribution_data = self._prepare_confidence_distribution(gene_qcs)

        # Prepare flag details
        flag_details = self._prepare_flag_details()

        # Prepare gene table (limited to avoid bloating HTML)
        gene_table = self._prepare_gene_table(gene_qcs)

        return ReportData(
            title=title or self.config.title,
            description=description or self.config.description,
            generated_at=datetime.now().isoformat(),
            summary=summary,
            tier_data=tier_data,
            severity_data=severity_data,
            category_data=category_data,
            flag_counts_data=flag_counts_data,
            aed_distribution_data=aed_distribution_data,
            confidence_distribution_data=confidence_distribution_data,
            flag_details=flag_details,
            gene_table=gene_table,
            config=self.config,
        )

    def _prepare_tier_data(self, summary: dict[str, Any]) -> dict[str, Any]:
        """Prepare tier distribution chart data."""
        tier_counts = summary.get("tier_counts", {})
        tiers = ["high", "medium", "low", "reject", "unclassified"]
        colors = {
            "high": "#28a745",
            "medium": "#ffc107",
            "low": "#fd7e14",
            "reject": "#dc3545",
            "unclassified": "#6c757d",
        }

        labels = []
        values = []
        chart_colors = []

        for tier in tiers:
            count = tier_counts.get(tier, 0)
            if count > 0:
                labels.append(tier.capitalize())
                values.append(count)
                chart_colors.append(colors.get(tier, "#6c757d"))

        return {
            "labels": labels,
            "values": values,
            "colors": chart_colors,
        }

    def _prepare_severity_data(self, summary: dict[str, Any]) -> dict[str, Any]:
        """Prepare severity distribution chart data."""
        severity_counts = summary.get("severity_counts", {})
        colors = {
            "info": "#17a2b8",
            "warning": "#ffc107",
            "error": "#fd7e14",
            "critical": "#dc3545",
        }

        labels = []
        values = []
        chart_colors = []

        for severity in ["info", "warning", "error", "critical"]:
            count = severity_counts.get(severity, 0)
            labels.append(severity.capitalize())
            values.append(count)
            chart_colors.append(colors.get(severity, "#6c757d"))

        return {
            "labels": labels,
            "values": values,
            "colors": chart_colors,
        }

    def _prepare_category_data(self, summary: dict[str, Any]) -> dict[str, Any]:
        """Prepare category distribution chart data."""
        category_counts = summary.get("category_counts", {})
        colors = {
            "confidence": "#007bff",
            "splice": "#6610f2",
            "homology": "#20c997",
            "structure": "#e83e8c",
            "annotation": "#fd7e14",
        }

        labels = []
        values = []
        chart_colors = []

        for category in ["confidence", "splice", "homology", "structure", "annotation"]:
            count = category_counts.get(category, 0)
            if count > 0:
                labels.append(category.capitalize())
                values.append(count)
                chart_colors.append(colors.get(category, "#6c757d"))

        return {
            "labels": labels,
            "values": values,
            "colors": chart_colors,
        }

    def _prepare_flag_counts_data(self, summary: dict[str, Any]) -> dict[str, Any]:
        """Prepare per-flag count chart data (top 15 flags)."""
        flag_counts = summary.get("flag_counts", {})

        # Sort by count descending, take top 15
        sorted_flags = sorted(flag_counts.items(), key=lambda x: x[1], reverse=True)[:15]

        labels = [f[0] for f in sorted_flags]
        values = [f[1] for f in sorted_flags]

        return {
            "labels": labels,
            "values": values,
            "total_flags": len(flag_counts),
        }

    def _prepare_aed_distribution(self, gene_qcs: dict[str, GeneQC]) -> dict[str, Any]:
        """Prepare AED score distribution data for histogram."""
        # Collect RNA-seq AED values
        rnaseq_aeds = []
        combined_aeds = []

        for qc in gene_qcs.values():
            if qc.rnaseq_aed is not None:
                rnaseq_aeds.append(qc.rnaseq_aed)
            if qc.combined_aed is not None:
                combined_aeds.append(qc.combined_aed)

        # Create histogram bins (0.0-0.1, 0.1-0.2, ..., 0.9-1.0)
        bins = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        bin_labels = ["0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
                      "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0"]

        def histogram(values: list[float], bins: list[float]) -> list[int]:
            counts = [0] * (len(bins) - 1)
            for v in values:
                for i in range(len(bins) - 1):
                    if bins[i] <= v < bins[i + 1]:
                        counts[i] += 1
                        break
                    elif v == 1.0 and i == len(bins) - 2:
                        counts[i] += 1
                        break
            return counts

        rnaseq_hist = histogram(rnaseq_aeds, bins) if rnaseq_aeds else [0] * 10
        combined_hist = histogram(combined_aeds, bins) if combined_aeds else [0] * 10

        return {
            "labels": bin_labels,
            "rnaseq_values": rnaseq_hist,
            "combined_values": combined_hist,
            "rnaseq_mean": sum(rnaseq_aeds) / len(rnaseq_aeds) if rnaseq_aeds else None,
            "combined_mean": sum(combined_aeds) / len(combined_aeds) if combined_aeds else None,
            "rnaseq_count": len(rnaseq_aeds),
            "combined_count": len(combined_aeds),
        }

    def _prepare_confidence_distribution(self, gene_qcs: dict[str, GeneQC]) -> dict[str, Any]:
        """Prepare confidence score distribution data for histogram."""
        scores = [qc.confidence_score for qc in gene_qcs.values() if qc.confidence_score is not None]

        # Create histogram bins
        bins = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        bin_labels = ["0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
                      "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0"]

        def histogram(values: list[float], bins: list[float]) -> list[int]:
            counts = [0] * (len(bins) - 1)
            for v in values:
                for i in range(len(bins) - 1):
                    if bins[i] <= v < bins[i + 1]:
                        counts[i] += 1
                        break
                    elif v == 1.0 and i == len(bins) - 2:
                        counts[i] += 1
                        break
            return counts

        hist = histogram(scores, bins) if scores else [0] * 10

        return {
            "labels": bin_labels,
            "values": hist,
            "mean": sum(scores) / len(scores) if scores else None,
            "count": len(scores),
        }

    def _prepare_flag_details(self) -> list[dict[str, Any]]:
        """Prepare flag details for display."""
        details = []
        for flag in Flags.get_all():
            details.append({
                "code": flag.code,
                "name": flag.name,
                "description": flag.description,
                "category": flag.category.value,
                "severity": flag.severity.value,
            })
        return sorted(details, key=lambda x: (x["category"], x["severity"], x["code"]))

    def _prepare_gene_table(
        self,
        gene_qcs: dict[str, GeneQC],
    ) -> list[dict[str, Any]]:
        """Prepare gene table data."""
        genes = sorted(gene_qcs.values(), key=lambda x: x.gene_id)

        # Limit if configured
        if self.config.max_genes_in_table:
            genes = genes[:self.config.max_genes_in_table]

        return [qc.to_dict() for qc in genes]

    def _render_html(self, data: ReportData) -> str:
        """Render HTML from report data.

        Uses embedded template to avoid external dependencies.
        """
        # Convert data to JSON for JavaScript
        tier_json = json.dumps(data.tier_data)
        severity_json = json.dumps(data.severity_data)
        category_json = json.dumps(data.category_data)
        flag_counts_json = json.dumps(data.flag_counts_data)
        aed_dist_json = json.dumps(data.aed_distribution_data)
        confidence_dist_json = json.dumps(data.confidence_distribution_data)
        gene_table_json = json.dumps(data.gene_table)

        # Build HTML
        html = self._get_html_template()

        # Replace placeholders
        html = html.replace("{{title}}", data.title)
        html = html.replace("{{description}}", data.description)
        html = html.replace("{{generated_at}}", data.generated_at)
        html = html.replace("{{total_genes}}", str(data.summary.get("total_genes", 0)))

        # Tier counts
        tier_counts = data.summary.get("tier_counts", {})
        html = html.replace("{{high_count}}", str(tier_counts.get("high", 0)))
        html = html.replace("{{medium_count}}", str(tier_counts.get("medium", 0)))
        html = html.replace("{{low_count}}", str(tier_counts.get("low", 0)))
        html = html.replace("{{reject_count}}", str(tier_counts.get("reject", 0)))

        # Percentages
        total = data.summary.get("total_genes", 1)
        if total == 0:
            total = 1
        html = html.replace("{{high_pct}}", f"{tier_counts.get('high', 0) / total * 100:.1f}")
        html = html.replace("{{medium_pct}}", f"{tier_counts.get('medium', 0) / total * 100:.1f}")
        html = html.replace("{{low_pct}}", f"{tier_counts.get('low', 0) / total * 100:.1f}")
        html = html.replace("{{reject_pct}}", f"{tier_counts.get('reject', 0) / total * 100:.1f}")

        # Average scores
        avg_conf = data.summary.get("avg_confidence")
        avg_splice = data.summary.get("avg_splice")
        avg_homology = data.summary.get("avg_homology")
        html = html.replace("{{avg_confidence}}", f"{avg_conf:.3f}" if avg_conf else "N/A")
        html = html.replace("{{avg_splice}}", f"{avg_splice:.3f}" if avg_splice else "N/A")
        html = html.replace("{{avg_homology}}", f"{avg_homology:.3f}" if avg_homology else "N/A")

        # Chart data
        html = html.replace("{{tier_data}}", tier_json)
        html = html.replace("{{severity_data}}", severity_json)
        html = html.replace("{{category_data}}", category_json)
        html = html.replace("{{flag_counts_data}}", flag_counts_json)
        html = html.replace("{{aed_dist_data}}", aed_dist_json)
        html = html.replace("{{confidence_dist_data}}", confidence_dist_json)
        html = html.replace("{{gene_table_data}}", gene_table_json)

        # AED summary stats
        aed_data = data.aed_distribution_data
        html = html.replace("{{rnaseq_aed_mean}}", f"{aed_data['rnaseq_mean']:.4f}" if aed_data['rnaseq_mean'] else "N/A")
        html = html.replace("{{combined_aed_mean}}", f"{aed_data['combined_mean']:.4f}" if aed_data['combined_mean'] else "N/A")

        # Gene table info
        html = html.replace("{{genes_in_table}}", str(len(data.gene_table)))
        html = html.replace("{{max_genes}}", str(data.config.max_genes_in_table))

        # Flag details table
        flag_rows = []
        for flag in data.flag_details:
            severity_class = f"severity-{flag['severity']}"
            flag_rows.append(f"""
                <tr class="{severity_class}">
                    <td><code>{flag['code']}</code></td>
                    <td>{flag['name']}</td>
                    <td>{flag['description']}</td>
                    <td>{flag['category']}</td>
                    <td><span class="badge {severity_class}">{flag['severity']}</span></td>
                </tr>
            """)
        html = html.replace("{{flag_rows}}", "\n".join(flag_rows))

        # Custom CSS
        html = html.replace("{{custom_css}}", data.config.custom_css)

        return html

    def _get_html_template(self) -> str:
        """Get the HTML template string."""
        return '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{title}}</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        :root {
            --primary: #007bff;
            --success: #28a745;
            --warning: #ffc107;
            --danger: #dc3545;
            --info: #17a2b8;
            --secondary: #6c757d;
        }

        * {
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            background-color: #f8f9fa;
            padding: 20px;
        }

        .container {
            max-width: 1400px;
            margin: 0 auto;
        }

        header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
        }

        header h1 {
            font-size: 2.5rem;
            margin-bottom: 10px;
        }

        header p {
            opacity: 0.9;
        }

        .meta {
            font-size: 0.9rem;
            opacity: 0.8;
            margin-top: 15px;
        }

        .section {
            background: white;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            padding: 25px;
            margin-bottom: 25px;
        }

        .section h2 {
            color: #333;
            border-bottom: 2px solid #eee;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }

        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
        }

        .stat-card {
            background: #f8f9fa;
            border-radius: 8px;
            padding: 20px;
            text-align: center;
        }

        .stat-card.high { border-left: 4px solid var(--success); }
        .stat-card.medium { border-left: 4px solid var(--warning); }
        .stat-card.low { border-left: 4px solid #fd7e14; }
        .stat-card.reject { border-left: 4px solid var(--danger); }

        .stat-value {
            font-size: 2rem;
            font-weight: bold;
            color: #333;
        }

        .stat-label {
            color: #666;
            font-size: 0.9rem;
        }

        .stat-pct {
            font-size: 0.85rem;
            color: #888;
        }

        .charts-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(350px, 1fr));
            gap: 25px;
        }

        .chart-container {
            background: #f8f9fa;
            border-radius: 8px;
            padding: 20px;
        }

        .chart-container h3 {
            text-align: center;
            margin-bottom: 15px;
            color: #555;
        }

        table {
            width: 100%;
            border-collapse: collapse;
        }

        th, td {
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #eee;
        }

        th {
            background: #f8f9fa;
            font-weight: 600;
            color: #333;
        }

        tr:hover {
            background-color: #f5f5f5;
        }

        .badge {
            display: inline-block;
            padding: 3px 8px;
            border-radius: 4px;
            font-size: 0.8rem;
            font-weight: 500;
        }

        .severity-info { background-color: #d1ecf1; color: #0c5460; }
        .severity-warning { background-color: #fff3cd; color: #856404; }
        .severity-error { background-color: #ffe5d0; color: #8a4500; }
        .severity-critical { background-color: #f8d7da; color: #721c24; }

        code {
            background: #e9ecef;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: monospace;
        }

        .scores-summary {
            display: flex;
            gap: 30px;
            justify-content: center;
            margin-top: 20px;
        }

        .score-item {
            text-align: center;
        }

        .score-value {
            font-size: 1.5rem;
            font-weight: bold;
            color: var(--primary);
        }

        .score-label {
            font-size: 0.85rem;
            color: #666;
        }

        .table-container {
            max-height: 500px;
            overflow-y: auto;
        }

        #geneFilter {
            width: 100%;
            padding: 10px 15px;
            border: 1px solid #ddd;
            border-radius: 5px;
            margin-bottom: 15px;
            font-size: 1rem;
        }

        {{custom_css}}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>{{title}}</h1>
            <p>{{description}}</p>
            <div class="meta">
                Generated: {{generated_at}} | Total Genes: {{total_genes}}
            </div>
        </header>

        <section class="section">
            <h2>Summary</h2>
            <div class="stats-grid">
                <div class="stat-card high">
                    <div class="stat-value">{{high_count}}</div>
                    <div class="stat-label">High Confidence</div>
                    <div class="stat-pct">{{high_pct}}%</div>
                </div>
                <div class="stat-card medium">
                    <div class="stat-value">{{medium_count}}</div>
                    <div class="stat-label">Medium Confidence</div>
                    <div class="stat-pct">{{medium_pct}}%</div>
                </div>
                <div class="stat-card low">
                    <div class="stat-value">{{low_count}}</div>
                    <div class="stat-label">Low Confidence</div>
                    <div class="stat-pct">{{low_pct}}%</div>
                </div>
                <div class="stat-card reject">
                    <div class="stat-value">{{reject_count}}</div>
                    <div class="stat-label">Rejected</div>
                    <div class="stat-pct">{{reject_pct}}%</div>
                </div>
            </div>

            <div class="scores-summary">
                <div class="score-item">
                    <div class="score-value">{{avg_confidence}}</div>
                    <div class="score-label">Avg Confidence</div>
                </div>
                <div class="score-item">
                    <div class="score-value">{{avg_splice}}</div>
                    <div class="score-label">Avg Splice Score</div>
                </div>
                <div class="score-item">
                    <div class="score-value">{{avg_homology}}</div>
                    <div class="score-label">Avg Homology Score</div>
                </div>
                <div class="score-item">
                    <div class="score-value">{{rnaseq_aed_mean}}</div>
                    <div class="score-label">Mean RNA-seq AED</div>
                </div>
                <div class="score-item">
                    <div class="score-value">{{combined_aed_mean}}</div>
                    <div class="score-label">Mean Combined AED</div>
                </div>
            </div>
        </section>

        <section class="section">
            <h2>Distribution Charts</h2>
            <div class="charts-grid">
                <div class="chart-container">
                    <h3>Tier Distribution</h3>
                    <canvas id="tierChart"></canvas>
                </div>
                <div class="chart-container">
                    <h3>Confidence Score Distribution</h3>
                    <canvas id="confidenceChart"></canvas>
                </div>
                <div class="chart-container">
                    <h3>AED Score Distribution</h3>
                    <canvas id="aedChart"></canvas>
                </div>
            </div>
        </section>

        <section class="section">
            <h2>Flag Analysis</h2>
            <div class="charts-grid">
                <div class="chart-container">
                    <h3>Top 15 Flags by Occurrence</h3>
                    <canvas id="flagCountsChart"></canvas>
                </div>
                <div class="chart-container">
                    <h3>Flag Severity Distribution</h3>
                    <canvas id="severityChart"></canvas>
                </div>
                <div class="chart-container">
                    <h3>Flag Category Distribution</h3>
                    <canvas id="categoryChart"></canvas>
                </div>
            </div>
        </section>

        <section class="section">
            <h2>QC Flags Reference</h2>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Code</th>
                            <th>Name</th>
                            <th>Description</th>
                            <th>Category</th>
                            <th>Severity</th>
                        </tr>
                    </thead>
                    <tbody>
                        {{flag_rows}}
                    </tbody>
                </table>
            </div>
        </section>

        <section class="section">
            <h2>Gene Details</h2>
            <p style="color: #666; margin-bottom: 15px;">
                Showing {{genes_in_table}} of {{total_genes}} genes (limited to {{max_genes}} for performance).
                Use the TSV output for complete data.
            </p>
            <input type="text" id="geneFilter" placeholder="Filter genes by ID, tier, or flag...">
            <div class="table-container">
                <table id="geneTable">
                    <thead>
                        <tr>
                            <th>Gene ID</th>
                            <th>Tier</th>
                            <th>Flags</th>
                            <th>Confidence</th>
                            <th>RNA-seq AED</th>
                            <th>Combined AED</th>
                        </tr>
                    </thead>
                    <tbody id="geneTableBody">
                    </tbody>
                </table>
            </div>
        </section>

        <footer style="text-align: center; color: #888; padding: 20px;">
            Generated by HelixForge QC System
        </footer>
    </div>

    <script>
        // Chart data
        const tierData = {{tier_data}};
        const severityData = {{severity_data}};
        const categoryData = {{category_data}};
        const flagCountsData = {{flag_counts_data}};
        const aedDistData = {{aed_dist_data}};
        const confidenceDistData = {{confidence_dist_data}};
        const geneTableData = {{gene_table_data}};

        // Tier Chart
        new Chart(document.getElementById('tierChart'), {
            type: 'doughnut',
            data: {
                labels: tierData.labels,
                datasets: [{
                    data: tierData.values,
                    backgroundColor: tierData.colors,
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    legend: { position: 'bottom' }
                }
            }
        });

        // Confidence Distribution Chart
        new Chart(document.getElementById('confidenceChart'), {
            type: 'bar',
            data: {
                labels: confidenceDistData.labels,
                datasets: [{
                    label: 'Gene Count',
                    data: confidenceDistData.values,
                    backgroundColor: '#007bff',
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    legend: { display: false },
                    title: { display: true, text: `Mean: ${confidenceDistData.mean ? confidenceDistData.mean.toFixed(3) : 'N/A'}` }
                },
                scales: { y: { beginAtZero: true } }
            }
        });

        // AED Distribution Chart
        new Chart(document.getElementById('aedChart'), {
            type: 'bar',
            data: {
                labels: aedDistData.labels,
                datasets: [
                    {
                        label: 'RNA-seq AED',
                        data: aedDistData.rnaseq_values,
                        backgroundColor: 'rgba(255, 99, 132, 0.7)',
                    },
                    {
                        label: 'Combined AED',
                        data: aedDistData.combined_values,
                        backgroundColor: 'rgba(54, 162, 235, 0.7)',
                    }
                ]
            },
            options: {
                responsive: true,
                plugins: {
                    legend: { position: 'top' },
                    title: { display: true, text: 'Lower AED = Better Support' }
                },
                scales: { y: { beginAtZero: true } }
            }
        });

        // Flag Counts Chart (horizontal bar for readability)
        new Chart(document.getElementById('flagCountsChart'), {
            type: 'bar',
            data: {
                labels: flagCountsData.labels,
                datasets: [{
                    label: 'Gene Count',
                    data: flagCountsData.values,
                    backgroundColor: '#6c757d',
                }]
            },
            options: {
                indexAxis: 'y',
                responsive: true,
                plugins: {
                    legend: { display: false },
                    title: { display: true, text: `Total unique flags: ${flagCountsData.total_flags}` }
                },
                scales: { x: { beginAtZero: true } }
            }
        });

        // Severity Chart
        new Chart(document.getElementById('severityChart'), {
            type: 'bar',
            data: {
                labels: severityData.labels,
                datasets: [{
                    label: 'Flag Count',
                    data: severityData.values,
                    backgroundColor: severityData.colors,
                }]
            },
            options: {
                responsive: true,
                plugins: { legend: { display: false } },
                scales: { y: { beginAtZero: true } }
            }
        });

        // Category Chart
        new Chart(document.getElementById('categoryChart'), {
            type: 'pie',
            data: {
                labels: categoryData.labels,
                datasets: [{
                    data: categoryData.values,
                    backgroundColor: categoryData.colors,
                }]
            },
            options: {
                responsive: true,
                plugins: { legend: { position: 'bottom' } }
            }
        });

        // Gene Table
        function renderGeneTable(data) {
            const tbody = document.getElementById('geneTableBody');
            tbody.innerHTML = '';

            data.forEach(gene => {
                const tr = document.createElement('tr');
                const tierClass = gene.tier === 'high' ? 'info' : gene.tier === 'medium' ? 'warning' : gene.tier === 'low' ? 'error' : 'critical';
                const flags = gene.flag_codes ? gene.flag_codes.slice(0, 3).join(', ') + (gene.flag_codes.length > 3 ? '...' : '') : '';
                tr.innerHTML = `
                    <td><code>${gene.gene_id}</code></td>
                    <td><span class="badge severity-${tierClass}">${gene.tier || 'N/A'}</span></td>
                    <td title="${gene.flag_codes ? gene.flag_codes.join(', ') : ''}">${flags}</td>
                    <td>${gene.confidence_score != null ? gene.confidence_score.toFixed(3) : 'N/A'}</td>
                    <td>${gene.rnaseq_aed != null ? gene.rnaseq_aed.toFixed(3) : 'N/A'}</td>
                    <td>${gene.combined_aed != null ? gene.combined_aed.toFixed(3) : 'N/A'}</td>
                `;
                tbody.appendChild(tr);
            });
        }

        // Initial render
        renderGeneTable(geneTableData);

        // Filter functionality
        document.getElementById('geneFilter').addEventListener('input', function(e) {
            const filter = e.target.value.toLowerCase();
            const filtered = geneTableData.filter(gene => {
                return gene.gene_id.toLowerCase().includes(filter) ||
                       (gene.tier && gene.tier.toLowerCase().includes(filter)) ||
                       (gene.flag_codes && gene.flag_codes.some(f => f.toLowerCase().includes(filter)));
            });
            renderGeneTable(filtered);
        });
    </script>
</body>
</html>'''


def generate_summary_report(
    gene_qcs: dict[str, GeneQC],
    output_path: Path | str,
) -> Path:
    """Generate a simple text summary report.

    Args:
        gene_qcs: Dict of GeneQC objects.
        output_path: Path to write report.

    Returns:
        Path to generated report.
    """
    summary = summarize_qc_results(gene_qcs)
    output_path = Path(output_path)

    lines = [
        "=" * 60,
        "HelixForge QC Summary Report",
        "=" * 60,
        "",
        f"Generated: {datetime.now().isoformat()}",
        f"Total Genes: {summary['total_genes']}",
        "",
        "Tier Distribution:",
        "-" * 30,
    ]

    tier_counts = summary.get("tier_counts", {})
    total = summary.get("total_genes", 1)
    if total == 0:
        total = 1
    for tier in ["high", "medium", "low", "reject", "unclassified"]:
        count = tier_counts.get(tier, 0)
        pct = count / total * 100 if total > 0 else 0
        lines.append(f"  {tier.capitalize():15} {count:6} ({pct:5.1f}%)")

    lines.extend([
        "",
        "Average Scores:",
        "-" * 30,
        f"  Confidence:    {summary.get('avg_confidence', 'N/A')}",
        f"  Splice:        {summary.get('avg_splice', 'N/A')}",
        f"  Homology:      {summary.get('avg_homology', 'N/A')}",
        "",
        "Flag Counts by Category:",
        "-" * 30,
    ])

    for cat in ["confidence", "splice", "homology", "structure", "annotation"]:
        count = summary.get("category_counts", {}).get(cat, 0)
        lines.append(f"  {cat.capitalize():15} {count:6}")

    lines.extend([
        "",
        "Flag Counts by Severity:",
        "-" * 30,
    ])

    for sev in ["info", "warning", "error", "critical"]:
        count = summary.get("severity_counts", {}).get(sev, 0)
        lines.append(f"  {sev.capitalize():15} {count:6}")

    lines.append("")
    lines.append("=" * 60)

    output_path.write_text("\n".join(lines))
    return output_path


def generate_json_report(
    gene_qcs: dict[str, GeneQC],
    output_path: Path | str,
) -> Path:
    """Generate a JSON QC report for programmatic access.

    Args:
        gene_qcs: Dict of GeneQC objects.
        output_path: Path to write report.

    Returns:
        Path to generated report.
    """
    summary = summarize_qc_results(gene_qcs)
    output_path = Path(output_path)

    report = {
        "generated_at": datetime.now().isoformat(),
        "summary": summary,
        "genes": {
            gene_id: qc.to_dict()
            for gene_id, qc in gene_qcs.items()
        },
    }

    output_path.write_text(json.dumps(report, indent=2))
    return output_path
