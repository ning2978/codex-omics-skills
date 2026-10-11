#!/usr/bin/env python3
"""Compatibility entrypoint; shared engine executes in this module namespace."""
from pathlib import Path
import sys
MODALITY = "xenium"
_SHARED = Path(__file__).resolve().parents[2] / "_shared" / "dotplot"
sys.path.insert(0, str(_SHARED))
exec(compile((_SHARED / "engine.py").read_bytes(), str(_SHARED / "engine.py"), "exec"), globals())
