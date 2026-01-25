"""Placeholder tests to verify test infrastructure is working.

These tests will be replaced with actual unit tests as modules are implemented.
"""

import pytest


class TestInfrastructure:
    """Tests to verify the test infrastructure is working correctly."""

    def test_pytest_runs(self) -> None:
        """Verify pytest is configured and running."""
        assert True

    def test_fixtures_available(self, test_data_dir) -> None:
        """Verify fixtures from conftest.py are available."""
        assert test_data_dir.name == "data"

    @pytest.mark.skip(reason="Placeholder for future implementation")
    def test_helixforge_import(self) -> None:
        """Verify helixforge package can be imported.

        TODO: Enable this test once package structure is complete.
        """
        import helixforge

        assert helixforge.__version__ == "0.1.0-alpha"


class TestDataStructures:
    """Placeholder tests for core data structures."""

    @pytest.mark.skip(reason="GeneModel not yet implemented")
    def test_gene_model_creation(self) -> None:
        """Test creating a GeneModel instance."""
        # TODO: Implement when GeneModel is defined
        pass

    @pytest.mark.skip(reason="Transcript not yet implemented")
    def test_transcript_creation(self) -> None:
        """Test creating a Transcript instance."""
        # TODO: Implement when Transcript is defined
        pass

    @pytest.mark.skip(reason="GenomicRegion not yet implemented")
    def test_genomic_region_creation(self) -> None:
        """Test creating a GenomicRegion instance."""
        # TODO: Implement when GenomicRegion is defined
        pass
