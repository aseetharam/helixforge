#!/usr/bin/env python3
"""Araport11 benchmark driver (Phase 17 D4) — the M7 accuracy/ablation run.

Thin CLI over ``helixforge.bench.driver.run_araport11_benchmark``: reconcile the
A. thaliana set, then benchmark + ablate against Araport11 + SwissProt and write
``docs/BENCHMARK.md`` (or a chosen path). The pipeline itself is the standard
Phase-8 ``run_pipeline(PipelineConfig)``; this script only wires the inputs.

Example (resume from the M2 Mikado artifacts, like ``scripts/m3_run.py``)::

    /usr/bin/python3 scripts/benchmark_araport11.py \
        --genome helixforge_testing/genome/athaliana.fasta \
        --helixer helixforge_testing/helixer_output/Arabidopsis-thaliana_helixer.gff3 \
        --helixer-h5 helixforge_testing/helixer_output/Arabidopsis-thaliana_input.h5 \
        --reconciled m3_out/helixforge.gff3 \
        --reference helixforge_testing/araport11/Araport11.gff3 \
        --proteins helixforge_testing/swissprot/uniprot_sprot.fasta \
        --lineage brassicales_odb10 \
        --out bench_out --doc docs/BENCHMARK.md

The actual numbers are produced by running this; the repo ships the code + the
``docs/BENCHMARK.md`` template with the metric rows defined.
"""

from __future__ import annotations

import argparse

from helixforge.bench.driver import run_araport11_benchmark
from helixforge.reconcile.pipeline import PipelineConfig


def main(argv=None):
    p = argparse.ArgumentParser(description="HelixForge Araport11 benchmark driver")
    p.add_argument("--genome", required=True)
    p.add_argument("--helixer", required=True, help="Helixer GFF3 (before set)")
    p.add_argument("--helixer-h5", default=None, help="Helixer HDF5 (Helixer-support column)")
    p.add_argument("--reconciled", default=None,
                   help="Existing HelixForge GFF3; if omitted the pipeline is run")
    p.add_argument("--reference", required=True, help="Araport11 reference GFF3")
    p.add_argument("--proteins", required=True, help="SwissProt protein FASTA")
    p.add_argument("--lineage", default=None, help="compleasm/BUSCO lineage")
    p.add_argument("--omadb", default=None, help="OMArk OMA database")
    p.add_argument("--out", default="bench_out", help="Output directory")
    p.add_argument("--doc", default="docs/BENCHMARK.md", help="Markdown output path")
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args(argv)

    config = PipelineConfig(
        genome_fasta=args.genome,
        helixer_gff3=args.helixer,
        helixer_h5=args.helixer_h5,
        protein_db=args.proteins,
        output_prefix=str(args.out) + "/helixforge",
    )
    before_after, benchmark, ablation = run_araport11_benchmark(
        config, args.reference, args.proteins, args.out,
        reconciled_gff3=args.reconciled, lineage=args.lineage, omadb=args.omadb,
        threads=args.threads, doc_path=args.doc,
    )
    print(f"wrote {args.doc}")
    print(f"before/after rows: {len(before_after)}; benchmark rows: {len(benchmark)}; "
          f"ablation variants: {len(ablation)}")


if __name__ == "__main__":
    main()
