from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def global_datadir():
    return Path(__file__).parent / "data"
