#!/usr/bin/env python3
"""Compatibility entry point for VCF-Polar v2.

The original v1.1 workflow acknowledged Ari's pipeline and
https://www.biostars.org/p/173338/. See README.md for migration notes.
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

try:
    from vcf_polar import main
except ModuleNotFoundError as exc:
    if exc.name != "pysam":
        raise
    sys.exit("VCF-Polar requires pysam. Install from the repository root: python -m pip install .")

if __name__ == "__main__":
    sys.exit(main())
