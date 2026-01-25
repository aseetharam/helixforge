# HelixForge

A modular toolkit for refining Helixer gene predictions into publication-quality annotations.

## Overview

HelixForge takes raw gene predictions from [Helixer](https://github.com/weberlab-hhu/Helixer) and refines them using:

- **RNA-seq evidence** for splice site validation and coverage-based confidence
- **Homology validation** against protein databases
- **Confidence scoring** with transparent, multi-factor metrics
- **Isoform reconstruction** from empirical splice junction data
- **Quality control** with comprehensive reports and filtering

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
# Refine Helixer predictions with RNA-seq evidence
helixforge refine \
    --predictions helixer_output.h5 \
    --genome genome.fa \
    --bam rnaseq.bam \
    --output refined_genes.gff3

# Generate QC report
helixforge qc \
    --gff refined_genes.gff3 \
    --output qc_report.html
```

> **Note**: Full functionality is under active development. See the roadmap below.

## Feature Roadmap

### v0.1.0 (Current Development)
- [ ] Core data structures for gene models
- [ ] GFF3/FASTA/HDF5 I/O
- [ ] Basic confidence scoring
- [ ] CLI skeleton
- [ ] QC flag definitions

### v0.2.0 (Planned)
- [ ] RNA-seq evidence integration
- [ ] Splice junction extraction
- [ ] Gene boundary refinement
- [ ] Merge/split detection
- [ ] HTML QC reports

### v0.3.0 (Future)
- [ ] Isoform reconstruction
- [ ] Homology validation
- [ ] Interactive visualization
- [ ] SLURM cluster support

## Documentation

Full documentation is available at [docs/](docs/) (coming soon).

For development guidelines, see [CLAUDE.md](CLAUDE.md).

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
