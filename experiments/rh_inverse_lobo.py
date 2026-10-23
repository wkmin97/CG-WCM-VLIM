#!/usr/bin/env python3
"""Separate inverse-calibration RH experiment; does not feed the conditional AGBD run."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common import DATA, RESULTS
from src.lobo_rh_inverse import run


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(DATA))
    parser.add_argument("--outdir", default=str(RESULTS/"rh_inverse_LOBO"))
    args = parser.parse_args()
    run(args.input, args.outdir)
