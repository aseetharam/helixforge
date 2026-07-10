# HelixForge containers — reproducibility image

`helixforge.def` is an [Apptainer](https://apptainer.org/) definition that locks
the **exact** CPU toolchain HelixForge v3 was validated on (M2/M3 on a real
*A. thaliana* genome). The pinned versions live in
[`docs/EXTERNAL_TOOLS.md`](../docs/EXTERNAL_TOOLS.md) — that file is the single
source of truth; this image and the `helixforge doctor` matrix track it.

An unpinned image defeats the purpose: the reconciliation stack is
**version-fragile** (Mikado 2.x needs pandas 1.5.3 / SQLAlchemy 1.4.49 — pandas
≥ 2 breaks `mikado serialise` with `'Engine' has no attribute 'cursor'`). Every
version in the `.def` is pinned for exactly this reason.

## What's inside (four isolated conda envs)

| Env | Tools (pinned) | Why isolated |
|---|---|---|
| `helixforge-env` | the `helixforge` package + `cli,stats,viz` extras (python 3.10) | provides the `helixforge` CLI + python; pandas ≥ 2 is fine here (never runs Mikado in-process) |
| `mikado-env` | mikado 2.3.4, pandas 1.5.3, sqlalchemy 1.4.49, transdecoder 5.5.0, diamond 2.1.16, portcullis | the fragile reconcile stack — its old pandas/SQLAlchemy must not leak into the package env |
| `prep-env` | STAR 2.7.11b, HISAT2 2.2.1, samtools 1.21, StringTie 3.0.3, miniprot 0.18 | raw-evidence staging (FASTQ→BAM→GTF, proteome→GFF) |
| `bench-env` | compleasm 0.2.6, BUSCO 5.7.1, OMArk 0.3.0, AGAT 1.4.1, gffcompare 0.12.6 | benchmarking; heavy/conflicting deps kept apart |

Each env's `bin/` is prepended to `PATH` in `%environment`, so every tool
resolves by its bare command name — exactly the names `helixforge doctor`
probes.

### The GPU-Helixer split (deliberate)

Helixer is **GPU-only** and **off the default path** (most users supply Helixer
output), so it is **not** in this CPU image — bundling its CUDA/TensorFlow stack
would bloat the image and pin it to a host driver. Run Helixer separately:

- a dedicated GPU image (`Helixer.py`, model `land_plant_v0.3`, pinned **0.3.5**), or
- a host module / existing GPU environment,

then feed its `*.gff3` + split HDF5 (`*_input.h5` + `*_predictions.h5`) into this
image. `io/hdf5.py` reads the native split layout directly (Phase 13 §C4).

## Build

```bash
# from the repo root (the %files section copies the whole source tree in)
apptainer build helixforge.sif containers/helixforge.def
```

Building Mikado + the benchmark tools pulls a lot from bioconda; expect a long
first build. Use `--fakeroot` if you lack root:

```bash
apptainer build --fakeroot helixforge.sif containers/helixforge.def
```

## Run

```bash
# the runscript is `helixforge`, so the image behaves like the CLI:
apptainer run helixforge.sif --help
apptainer run helixforge.sif doctor                       # preflight (see below)

# or exec a specific tool / a full pipeline:
apptainer exec helixforge.sif mikado --version
apptainer exec helixforge.sif helixforge run --help
```

Bind your data and a **fast local scratch** for the Mikado SQLite DB (never a
network FS — `serialise` is single-threaded and DB-insertion bound):

```bash
apptainer exec \
  --bind /path/to/data:/data \
  --bind /fast/scratch:/scratch \
  helixforge.sif helixforge run --workdir /scratch ...
```

## Verifying versions *inside* the image — `helixforge doctor`

`helixforge doctor` (Phase 13) is the regression-gate preflight. It resolves
every external tool on `PATH` (or via a `--*-bin` override), probes each
`--version`, and prints a found / missing / mismatch table **against the pinned
matrix** in `src/helixforge/prep/doctor.py` (which mirrors
`docs/EXTERNAL_TOOLS.md`):

```bash
apptainer exec helixforge.sif helixforge doctor
```

Expected inside this image: every reconcile + prep + benchmark tool `found` at
its pinned version; `helixer` shows `missing` (the intended GPU split — it is
not a required tool, so the report stays OK). A `mismatch` means the image
drifted from the matrix and must be rebuilt; a missing **required** tool
(mikado / diamond / transdecoder) makes `doctor` exit non-zero.

The `%test` section runs `helixforge doctor` and asserts the pandas 1.5.3 pin at
build time, so a broken image fails `apptainer test helixforge.sif`.

## Keeping pins in sync

When a pin changes, update **`docs/EXTERNAL_TOOLS.md` first** (source of truth),
then this `.def` and the `TOOL_MATRIX` in `src/helixforge/prep/doctor.py` in the
same commit. `helixforge doctor` will flag any drift between the three.
