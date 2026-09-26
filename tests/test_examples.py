"""
Tests for the example scripts in the examples/ directory. Very slow.

Warnings
--------
Work in progress.
"""

import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES_DIR = Path(__file__).parent.parent / "examples"


@pytest.mark.slow
@pytest.mark.skip(reason="Test not yet finished.")
# @pytest.mark.parametrize(
#     "script", sorted([s for s in EXAMPLES_DIR.rglob("*.py")])
# )
def test_scripts(script: Path) -> None:
    """Tests all example scripts contained in examples/."""
    result = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True,
        text=True,
        check=False,
        cwd=script.parent,
    )
    assert result.returncode == 0
