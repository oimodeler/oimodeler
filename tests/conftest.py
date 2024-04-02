from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def global_data_dir():
    """Return the global data directory."""
    return Path(__file__).parent.parent / "data"


@pytest.fixture(scope="session")
def real_data_dir(global_data_dir: Path):
    """Return the global data directory."""
    return global_data_dir / "RealData"
