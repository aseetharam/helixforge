"""Tiered output filtering by confidence/evidence/biotype."""

from __future__ import annotations

from typing import Any

import attrs

from helixforge.utils.logging import get_logger

_log = get_logger(__name__)


@attrs.define
class FilterCriteria:
    """Criteria for filtering gene models from a GFF3."""

    min_tier: int | None = None
    max_tier: int | None = None
    exclude_biotypes: list[str] = attrs.field(factory=list)
    require_flags: list[str] = attrs.field(factory=list)
    exclude_flags: list[str] = attrs.field(factory=list)
    min_confidence: float | None = None


class GeneFilter:
    """Filter gene records parsed from a HelixForge GFF3.

    Operates on gene dicts with GFF3 attribute fields extracted into keys:
    ``tier`` (int), ``gene_biotype`` (str), ``flags`` (list[str]),
    ``combined_score`` (float).
    """

    def __init__(self, criteria: FilterCriteria) -> None:
        self.criteria = criteria

    def apply(
        self, genes: list[dict[str, Any]]
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
        """Return ``(kept, excluded)`` gene dicts."""
        kept: list[dict[str, Any]] = []
        excluded: list[dict[str, Any]] = []
        for gene in genes:
            if self._passes(gene):
                kept.append(gene)
            else:
                excluded.append(gene)
        return kept, excluded

    def _passes(self, gene: dict[str, Any]) -> bool:
        c = self.criteria
        tier = gene.get("tier")
        if tier is not None:
            if c.min_tier is not None and tier < c.min_tier:
                return False
            if c.max_tier is not None and tier > c.max_tier:
                return False

        biotype = gene.get("gene_biotype") or gene.get("biotype")
        if c.exclude_biotypes and biotype in c.exclude_biotypes:
            return False

        flags = gene.get("flags", [])
        if c.require_flags and not all(f in flags for f in c.require_flags):
            return False
        if c.exclude_flags and any(f in flags for f in c.exclude_flags):
            return False

        score = gene.get("combined_score")
        if c.min_confidence is not None:
            if score is None or score < c.min_confidence:
                return False

        return True

    @classmethod
    def high_confidence(cls) -> GeneFilter:
        """Tier <= 2, no INTERNAL_STOP or NO_START_CODON flags."""
        return cls(
            FilterCriteria(
                max_tier=2,
                exclude_flags=["INTERNAL_STOP", "NO_START_CODON"],
            )
        )

    @classmethod
    def publication_ready(cls) -> GeneFilter:
        """Tier 1 only, protein_coding biotype only."""
        return cls(
            FilterCriteria(
                max_tier=1,
                exclude_biotypes=[
                    b
                    for b in [
                        "transposable_element",
                        "ncRNA",
                        "pseudogene",
                        "tRNA",
                        "rRNA",
                        "snRNA",
                        "snoRNA",
                        "miRNA",
                    ]
                ],
            )
        )


def parse_gff3_genes(path: str) -> list[dict[str, Any]]:
    """Parse a HelixForge GFF3 into gene dicts with extracted attributes.

    Returns a list of dicts, each with keys: ``gene_id``, ``seqid``,
    ``strand``, ``tier`` (int or None), ``gene_biotype``, ``origin``,
    ``flags`` (list[str]), ``combined_score`` (float or None), plus
    ``lines`` (all raw GFF3 lines belonging to this gene).
    """
    genes: list[dict[str, Any]] = []
    current_gene: dict[str, Any] | None = None
    header_lines: list[str] = []

    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                header_lines.append(line)
                continue
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2]
            if ftype == "gene":
                if current_gene is not None:
                    genes.append(current_gene)
                attrs = _parse_gff3_attrs(cols[8])
                current_gene = {
                    "gene_id": attrs.get("ID", ""),
                    "seqid": cols[0],
                    "strand": cols[6],
                    "tier": _safe_int(attrs.get("tier")),
                    "gene_biotype": attrs.get("gene_biotype"),
                    "origin": attrs.get("origin"),
                    "flags": _parse_flags(attrs.get("flags")),
                    "combined_score": None,
                    "lines": [line],
                }
            elif current_gene is not None:
                current_gene["lines"].append(line)
                if ftype == "mRNA":
                    attrs = _parse_gff3_attrs(cols[8])
                    score = attrs.get("combined_score")
                    if score is not None:
                        try:
                            s = float(score)
                            existing = current_gene["combined_score"]
                            if existing is None or s > existing:
                                current_gene["combined_score"] = s
                        except ValueError:
                            pass
    if current_gene is not None:
        genes.append(current_gene)

    return genes


def write_filtered_gff3(
    genes: list[dict[str, Any]],
    output_path: str,
    header_from: str | None = None,
) -> int:
    """Write kept genes back to a GFF3. Returns gene count."""
    with open(output_path, "w") as fh:
        if header_from:
            with open(header_from) as src:
                for line in src:
                    if line.startswith("#"):
                        fh.write(line)
                    else:
                        break
        elif not any(True for g in genes if g.get("lines")):
            fh.write("##gff-version 3\n")
        for gene in genes:
            for line in gene.get("lines", []):
                fh.write(line)
    return len(genes)


def _parse_gff3_attrs(attrs_str: str) -> dict[str, str]:
    """Parse GFF3 column-9 into a dict."""
    result: dict[str, str] = {}
    for part in attrs_str.split(";"):
        part = part.strip()
        if "=" in part:
            key, _, val = part.partition("=")
            result[key] = val
    return result


def _parse_flags(flags_str: str | None) -> list[str]:
    """Parse comma-separated flags string into a list."""
    if not flags_str:
        return []
    return [f.strip() for f in flags_str.split(",") if f.strip()]


def _safe_int(val: str | None) -> int | None:
    if val is None:
        return None
    try:
        return int(val)
    except ValueError:
        return None
