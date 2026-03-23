"""Helpers for running example scripts directly from a repo checkout."""

import sys
from pathlib import Path


def bootstrap_example_paths(file_path: str) -> None:
    """Ensure example scripts import from the local checkout's `src/` tree."""
    examples_dir = Path(file_path).resolve().parent
    src_dir = examples_dir.parent / "src"

    for path in (examples_dir, src_dir):
        path_str = str(path)
        if path_str not in sys.path:
            sys.path.insert(0, path_str)
