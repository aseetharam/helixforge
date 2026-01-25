"""Interactive visualization using Dash/Panel.

This module provides interactive genome browsers using Dash or Panel
for exploring gene predictions and evidence.

Features:
    - Interactive gene browser
    - Zoom and pan navigation
    - Evidence track toggling
    - Gene filtering and search
    - Export capabilities

Example:
    >>> from helixforge.viz.interactive import launch_browser
    >>> launch_browser(genes, evidence, port=8050)

TODO:
    - Implement Dash-based browser
    - Add Panel alternative
    - Add search functionality
    - Add filtering controls
    - Support for large genomes
"""

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from helixforge.core.models import GeneModel
    from helixforge.isoforms.evidence import EvidenceCollection

# =============================================================================
# Browser Launch
# =============================================================================


def launch_browser(
    genes: list["GeneModel"],
    evidence: "EvidenceCollection | None" = None,
    genome_path: Path | str | None = None,
    host: str = "127.0.0.1",
    port: int = 8050,
    debug: bool = False,
) -> None:
    """Launch interactive genome browser.

    Args:
        genes: Gene models to display.
        evidence: Optional evidence to display.
        genome_path: Path to genome FASTA for sequence display.
        host: Host address to bind to.
        port: Port number.
        debug: Enable debug mode.
    """
    # TODO: Implement browser launch
    raise NotImplementedError("launch_browser not yet implemented")


# =============================================================================
# Dash App
# =============================================================================


def create_dash_app(
    genes: list["GeneModel"],
    evidence: "EvidenceCollection | None" = None,
) -> "dash.Dash":  # type: ignore
    """Create a Dash application for genome browsing.

    Args:
        genes: Gene models to display.
        evidence: Optional evidence to display.

    Returns:
        Dash application instance.
    """
    # TODO: Implement Dash app
    raise NotImplementedError("create_dash_app not yet implemented")


def create_gene_view_component(gene: "GeneModel") -> "dash.html.Div":  # type: ignore
    """Create a Dash component for displaying a gene.

    Args:
        gene: Gene model to display.

    Returns:
        Dash HTML component.
    """
    # TODO: Implement gene view component
    raise NotImplementedError("create_gene_view_component not yet implemented")


def create_coverage_component(
    evidence: "EvidenceCollection",
) -> "dash.dcc.Graph":  # type: ignore
    """Create a Plotly coverage track component.

    Args:
        evidence: Evidence with coverage data.

    Returns:
        Dash Graph component.
    """
    # TODO: Implement coverage component
    raise NotImplementedError("create_coverage_component not yet implemented")


# =============================================================================
# Panel App
# =============================================================================


def create_panel_app(
    genes: list["GeneModel"],
    evidence: "EvidenceCollection | None" = None,
) -> "panel.viewable.Viewable":  # type: ignore
    """Create a Panel application for genome browsing.

    Args:
        genes: Gene models to display.
        evidence: Optional evidence to display.

    Returns:
        Panel viewable instance.
    """
    # TODO: Implement Panel app
    raise NotImplementedError("create_panel_app not yet implemented")


# =============================================================================
# Plotly Figures
# =============================================================================


def create_locus_figure(
    gene: "GeneModel",
    evidence: "EvidenceCollection | None" = None,
) -> "plotly.graph_objects.Figure":  # type: ignore
    """Create a Plotly figure for a gene locus.

    Args:
        gene: Gene model to display.
        evidence: Optional evidence to display.

    Returns:
        Plotly Figure object.
    """
    # TODO: Implement Plotly locus figure
    raise NotImplementedError("create_locus_figure not yet implemented")


def create_junction_figure(
    evidence: "EvidenceCollection",
) -> "plotly.graph_objects.Figure":  # type: ignore
    """Create a Plotly figure for splice junctions.

    Args:
        evidence: Evidence with junction data.

    Returns:
        Plotly Figure object.
    """
    # TODO: Implement Plotly junction figure
    raise NotImplementedError("create_junction_figure not yet implemented")
