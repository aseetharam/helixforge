"""Phase 23 D1/D2 — nested sub-configs, from_flat, YAML round-trip. Floor: 5.

The grouping must not break any flat construction path (CLAUDE.md §15 / the
phase's "What NOT to do" #1): ``PipelineConfig(min_cds_overlap=…)``,
``dataclasses.replace(cfg, pad=…)`` and the CLI's flat kwargs all keep working,
while the new sub-config validators reject out-of-range knobs and the nested
structure round-trips through ``to_nested_dict`` / YAML.
"""

import dataclasses

import pytest

from helixforge.reconcile.pipeline import (
    ASConfig,
    ClassificationConfig,
    PipelineConfig,
    ReconcileConfig,
    ResourceConfig,
    ValidationConfig,
)


def _cfg(**kw):
    return PipelineConfig(genome_fasta="g.fa", helixer_gff3="h.gff3", **kw)


# --- from_flat reproduces the flat config -----------------------------------

def test_from_flat_equals_flat_constructor():
    flat = _cfg(min_cds_overlap=0.7, pad=False, max_isoforms=8, procs=3,
                min_tpm=1.5, junction_min_reads=5)
    built = PipelineConfig.from_flat(
        genome_fasta="g.fa", helixer_gff3="h.gff3", min_cds_overlap=0.7,
        pad=False, max_isoforms=8, procs=3, min_tpm=1.5, junction_min_reads=5,
    )
    assert built == flat


def test_from_flat_minimal_quality_gate():
    # The phase's quality-gate snippet must work.
    cfg = PipelineConfig.from_flat(genome_fasta="g", helixer_gff3="h")
    assert cfg.genome_fasta == "g" and cfg.min_cds_overlap == 0.6


def test_from_flat_expands_nested_groups():
    cfg = PipelineConfig.from_flat(
        genome_fasta="g.fa", helixer_gff3="h.gff3",
        reconciliation={"min_cds_overlap": 0.9},
        resources=ResourceConfig(procs=4, threads=9),
    )
    assert cfg.min_cds_overlap == 0.9
    assert cfg.procs == 4 and cfg.threads == 9


# --- grouped views project the flat fields ----------------------------------

def test_grouped_views_match_flat_fields():
    cfg = _cfg(min_cds_overlap=0.7, pad=False, max_isoforms=8, procs=3,
               near_zero_coverage=0.2, junction_tolerance=4)
    assert cfg.classification == ClassificationConfig(0.5, 1, 2.0, 0.2)
    assert cfg.as_config == ASConfig(True, True, 8, False, False, True, 200)
    assert cfg.reconciliation == ReconcileConfig(0.5, 0.7, 0.6, False, None)
    assert cfg.validation == ValidationConfig(300, 10, 100_000, 4, 3)
    assert cfg.resources == ResourceConfig(3, 4)


def test_dataclasses_replace_still_works_on_flat_fields():
    # The ablation runner + chunk builder rely on dataclasses.replace(cfg, …).
    cfg = _cfg()
    replaced = dataclasses.replace(cfg, pad=False, scoring_profile="permissive",
                                  helixer_support_weight=0.0, min_cds_overlap=0.8)
    assert replaced.pad is False
    assert replaced.scoring_profile == "permissive"
    assert replaced.as_config.pad is False
    assert replaced.reconciliation.min_cds_overlap == 0.8


# --- validators reject out-of-range -----------------------------------------

@pytest.mark.parametrize("kw", [
    {"min_cds_overlap": 2.0},
    {"min_cdna_overlap": -0.1},
    {"reciprocal_overlap": 1.5},
    {"near_zero_coverage": 5.0},
    {"max_isoforms": 0},
    {"procs": 0},
    {"flank": -1},
    {"min_tpm": -2.0},
])
def test_out_of_range_knob_rejected(kw):
    with pytest.raises(ValueError):
        _cfg(**kw)


def test_subconfig_validator_direct():
    with pytest.raises(ValueError):
        ReconcileConfig(min_cds_overlap=1.2)
    with pytest.raises(ValueError):
        ASConfig(max_isoforms=0)
    with pytest.raises(ValueError):
        ClassificationConfig(near_zero_coverage=2.0)


# --- nested-dict + YAML round-trip ------------------------------------------

def test_nested_dict_round_trip():
    cfg = _cfg(min_cds_overlap=0.7, pad=False, max_isoforms=8, procs=3,
               region="Chr1", stringtie_list=["a.gtf", "b.gtf"],
               novel_evidence_floor=1.5)
    nested = cfg.to_nested_dict()
    # the five groups are present and nested
    assert set(nested["reconciliation"]) == {
        "reciprocal_overlap", "min_cds_overlap", "min_cdna_overlap",
        "admit_novel", "novel_evidence_floor",
        # Phase 29 D1 paralog/tandem-array merge guards.
        "merge_min_gap_reads", "merge_require_canonical",
        "paralog_identity_threshold", "paralog_kmer_k",
    }
    assert PipelineConfig.from_flat(**nested) == cfg


def test_yaml_round_trip(tmp_path):
    yaml = pytest.importorskip("yaml")  # noqa: F841 - degrade if PyYAML absent
    cfg = _cfg(min_cds_overlap=0.7, pad=False, max_isoforms=8, procs=3,
               region="Chr1:100-200", junction_min_reads=7)
    path = tmp_path / "run.yaml"
    assert cfg.to_yaml(path) is not None
    assert PipelineConfig.from_yaml(path) == cfg


def test_region_tuple_round_trips_through_yaml_shape():
    # An explicit-span region tuple survives the list-ification a YAML reload does.
    cfg = _cfg(region=("Chr1", 100, 200))
    reloaded = PipelineConfig.from_flat(**cfg.to_nested_dict())
    assert reloaded.region == ("Chr1", 100, 200)
    assert reloaded == cfg
