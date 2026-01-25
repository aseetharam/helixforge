"""Configuration management for HelixForge.

This module handles loading, validating, and providing access to
HelixForge configuration settings. Configuration can come from:
- Default values
- Configuration files (YAML/TOML)
- Environment variables
- Command-line arguments

Example:
    >>> from helixforge.config import Config
    >>> config = Config.load("helixforge.yaml")
    >>> config.parallel.chunk_size
    1000000

TODO:
    - Implement Config class with attrs
    - Add configuration file loading (YAML/TOML)
    - Add environment variable overrides
    - Add validation for all settings
    - Add preset configurations for common use cases
"""

from pathlib import Path
from typing import Any

import attrs

# =============================================================================
# Default Configuration Values
# =============================================================================

# Parallel processing defaults
DEFAULT_CHUNK_SIZE = 1_000_000  # Base pairs per chunk
DEFAULT_MAX_WORKERS = 1
DEFAULT_MEMORY_LIMIT_MB = 8192

# Confidence scoring defaults
DEFAULT_MIN_CONFIDENCE = 0.5
DEFAULT_CONFIDENCE_WEIGHTS = {
    "sequence": 0.3,
    "splice": 0.25,
    "coverage": 0.25,
    "homology": 0.2,
}

# QC thresholds
DEFAULT_MIN_EXON_LENGTH = 3
DEFAULT_MAX_INTRON_LENGTH = 500_000
DEFAULT_MIN_CDS_LENGTH = 30

# Evidence thresholds
DEFAULT_MIN_JUNCTION_READS = 3
DEFAULT_MIN_COVERAGE = 1.0


# =============================================================================
# Configuration Classes
# =============================================================================


@attrs.define
class ParallelConfig:
    """Configuration for parallel processing.

    Attributes:
        chunk_size: Size of genomic chunks in base pairs.
        max_workers: Maximum number of parallel workers.
        memory_limit_mb: Memory limit per worker in MB.
        use_slurm: Whether to use SLURM for cluster execution.
    """

    chunk_size: int = DEFAULT_CHUNK_SIZE
    max_workers: int = DEFAULT_MAX_WORKERS
    memory_limit_mb: int = DEFAULT_MEMORY_LIMIT_MB
    use_slurm: bool = False

    # TODO: Add validation
    # TODO: Add SLURM-specific settings


@attrs.define
class ConfidenceConfig:
    """Configuration for confidence scoring.

    Attributes:
        min_confidence: Minimum confidence to retain a gene.
        weights: Component weights for confidence calculation.
    """

    min_confidence: float = DEFAULT_MIN_CONFIDENCE
    weights: dict[str, float] = attrs.Factory(lambda: DEFAULT_CONFIDENCE_WEIGHTS.copy())

    # TODO: Add validation that weights sum to 1.0


@attrs.define
class QCConfig:
    """Configuration for quality control.

    Attributes:
        min_exon_length: Minimum exon length in base pairs.
        max_intron_length: Maximum intron length in base pairs.
        min_cds_length: Minimum CDS length in base pairs.
    """

    min_exon_length: int = DEFAULT_MIN_EXON_LENGTH
    max_intron_length: int = DEFAULT_MAX_INTRON_LENGTH
    min_cds_length: int = DEFAULT_MIN_CDS_LENGTH

    # TODO: Add additional QC thresholds


@attrs.define
class EvidenceConfig:
    """Configuration for evidence integration.

    Attributes:
        min_junction_reads: Minimum reads to support a splice junction.
        min_coverage: Minimum coverage to consider a region expressed.
    """

    min_junction_reads: int = DEFAULT_MIN_JUNCTION_READS
    min_coverage: float = DEFAULT_MIN_COVERAGE

    # TODO: Add strand-specific settings


@attrs.define
class Config:
    """Main configuration container for HelixForge.

    Attributes:
        parallel: Parallel processing configuration.
        confidence: Confidence scoring configuration.
        qc: Quality control configuration.
        evidence: Evidence integration configuration.
    """

    parallel: ParallelConfig = attrs.Factory(ParallelConfig)
    confidence: ConfidenceConfig = attrs.Factory(ConfidenceConfig)
    qc: QCConfig = attrs.Factory(QCConfig)
    evidence: EvidenceConfig = attrs.Factory(EvidenceConfig)

    @classmethod
    def load(cls, path: Path | str | None = None) -> "Config":
        """Load configuration from file.

        Args:
            path: Path to configuration file (YAML or TOML).
                  If None, returns default configuration.

        Returns:
            Loaded configuration object.

        Raises:
            FileNotFoundError: If configuration file doesn't exist.
            ValueError: If configuration file is invalid.
        """
        # TODO: Implement file loading
        if path is None:
            return cls()

        # Placeholder for file loading
        raise NotImplementedError("Configuration file loading not yet implemented")

    def to_dict(self) -> dict[str, Any]:
        """Convert configuration to dictionary.

        Returns:
            Dictionary representation of configuration.
        """
        return attrs.asdict(self)

    def save(self, path: Path | str) -> None:
        """Save configuration to file.

        Args:
            path: Path to save configuration file.
        """
        # TODO: Implement file saving
        raise NotImplementedError("Configuration file saving not yet implemented")
