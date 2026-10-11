"""Compatibility CLI for the real-Seurat style and numeric regression suite."""
from pathlib import Path
import runpy
runpy.run_path(str(Path(__file__).with_name("verify_style.py")), run_name="__main__")
