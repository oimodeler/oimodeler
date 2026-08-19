"""
Tests for the example scripts in the examples/ directory. Very slow.

Warnings
--------
Work in progress.
"""

import subprocess
import sys
from pathlib import Path

import papermill as pm
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


@pytest.mark.slow
# @pytest.mark.parametrize(
#     "notebook",
#     sorted([nb for nb in (EXAMPLES_DIR / "notebooks").rglob("*.ipynb")]),
# )
@pytest.mark.skip(reason="Test not yet finished.")
def test_notebooks(notebook: Path) -> None:
    """Tests all example Jupyter notebooks contained in examples/."""
    output_notebook = notebook.with_suffix(".out.ipynb")
    try:
        pm.execute_notebook(str(notebook), str(output_notebook))
        assert output_notebook.exists()
    finally:
        output_notebook.unlink(missing_ok=True)
