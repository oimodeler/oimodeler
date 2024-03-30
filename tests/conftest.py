from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def global_datadir():
    """Return the global data directory."""
    return Path(__file__).parent / "data"
