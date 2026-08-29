"""
Configuration for the tests.
"""

from pathlib import Path

import pytest

# TODO: Make more tests for dictionaries in the serialisations?


@pytest.fixture(scope="session")
def package_dir() -> Path:
    """Return the global data directory."""
    return Path(__file__).parent.parent


@pytest.fixture(scope="session")
def test_data_dir() -> Path:
    """Return the global data directory."""
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def global_data_dir(package_dir: Path) -> Path:
    """Return the global data directory."""
    return package_dir / "data"


@pytest.fixture(scope="session")
def real_data_dir(global_data_dir: Path) -> Path:
    """Return the global data directory."""
    return global_data_dir / "RealData"
