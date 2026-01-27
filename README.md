# HelixForge

A modular toolkit for refining Helixer gene predictions into publication-quality annotations.

## Overview

HelixForge takes raw gene predictions from [Helixer](https://github.com/weberlab-hhu/Helixer) and refines them using:

- **RNA-seq evidence** for splice site validation and coverage-based confidence
- **Homology validation** against protein databases
- **Confidence scoring** with transparent, multi-factor metrics from HDF5 predictions
- **Quality control** with comprehensive reports and tiered filtering

## Installation

### From Source (Development)

```bash
git clone https://github.com/aseetharam/helixforge.git
cd helixforge
pip install -e ".[dev]"
```

### Dependencies

HelixForge requires Python 3.10 or later and depends on:
- h5py, pysam, pyfaidx for file I/O
- numpy, pandas, scipy for data processing
- matplotlib, plotly for visualization
- click, rich for CLI

## Quick Start

```bash
# Main workflow: refine Helixer predictions with RNA-seq evidence
helixforge refine \
    -p helixer_predictions.h5 \
    -g helixer_predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    -o refined.gff3 \
    -r refine_report.tsv

# Validate with protein homology
helixforge homology extract-proteins --gff refined.gff3 --genome genome.fa -o proteins.fa
helixforge homology search --proteins proteins.fa --database swissprot.dmnd -o hits.tsv
helixforge homology validate --search-results hits.tsv --gff refined.gff3 -o validation.tsv

# Aggregate QC results and generate report
helixforge qc aggregate --refine-tsv refine_report.tsv --homology-tsv validation.tsv -o qc.tsv
helixforge qc report --qc-tsv qc.tsv -o qc_report.html

# Create tiered output GFF files
helixforge qc tiered-output --qc-results qc.tsv --gff refined.gff3 -o tiered/
```

## Available Commands

```
helixforge
├── refine        # Main pipeline - splice correction, boundaries, confidence, evidence
├── confidence    # Standalone confidence scoring from HDF5
├── evidence      # Standalone RNA-seq evidence scoring
├── homology      # Protein homology validation
│   ├── extract-proteins
│   ├── search
│   └── validate
├── validate      # Quick protein validation (all-in-one)
├── qc            # Quality control and reporting
│   ├── aggregate
│   ├── report
│   ├── filter
│   ├── tiered-output
│   └── list-flags
├── parallel      # Large genome parallelization
│   ├── plan
│   ├── tasks
│   ├── aggregate
│   └── suggest
└── viz           # Visualization (optional)
```

## Pipeline Diagram

```
                          ┌─────────────────────┐
                          │   Helixer Output    │
                          │  predictions.h5     │
                          │  predictions.gff3   │
                          └──────────┬──────────┘
                                     │
                                     │  + RNA-seq BAM(s)
                                     ▼
                          ┌─────────────────────┐
                          │       refine        │
                          │  (main pipeline)    │
                          │                     │
                          │ • Splice correction │
                          │ • Boundary adjust   │
                          │ • Confidence score  │
                          │ • Evidence score    │
                          └──────────┬──────────┘
                                     │
                                     ▼
                          ┌─────────────────────┐
                          │   homology validate │
                          │ (protein homology)  │
                          └──────────┬──────────┘
                                     │
                                     ▼
                          ┌─────────────────────┐
                          │    qc aggregate     │
                          │    qc report        │
                          └──────────┬──────────┘
                                     │
                                     ▼
                          ┌─────────────────────┐
                          │   qc tiered-output  │
                          │   (final GFFs)      │
                          └─────────────────────┘
```

## Documentation

- [Workflow Guide](docs/workflow.md) - Complete step-by-step workflow
- [Refine Pipeline](docs/refine.md) - Main refinement pipeline details
- [Evidence Scoring](docs/evidence.md) - RNA-seq evidence scoring
- [Parallel Execution](docs/parallel.md) - HPC parallelization guide

For development guidelines, see [GUIDE.md](GUIDE.md).

## Development Status

| Phase | Description | Status |
|-------|-------------|--------|
| 0 | Project scaffold | COMPLETE |
| 1 | I/O Foundation | COMPLETE |
| 2 | Confidence Scoring | COMPLETE |
| 3 | Splice Refinement | COMPLETE |
| 4 | Parallel Execution | COMPLETE |
| 5 | Homology Pipeline | COMPLETE |
| 6 | QC System | COMPLETE |
| 7 | Visualization | PENDING |
| 8 | CLI Integration | COMPLETE |
| 9 | Testing/Documentation | IN PROGRESS |
| 10 | Isoform Support (v0.2) | PENDING |

## Contributing

We welcome contributions! Please see our contributing guidelines (coming soon).

### Development Setup

```bash
# Install with development dependencies
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install

# Run tests
pytest -v

# Run linter
ruff check src/

# Type checking
mypy src/helixforge
```

### Code Style

- Type hints required on all public functions
- Google-style docstrings
- 100 character line length
- Ruff for linting and formatting

## Citation

If you use HelixForge in your research, please cite:

> HelixForge: A toolkit for refining Helixer gene predictions. (2024). https://github.com/aseetharam/helixforge

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

HelixForge builds upon the excellent work of the Helixer team and the broader genomics community.
