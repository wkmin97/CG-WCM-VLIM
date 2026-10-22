#!/usr/bin/env python3
"""Run the documented processed-data reproduction, then verify its outputs."""
import argparse
import os
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-inverse-rh", action="store_true", help="Also run the separate inverse-calibration RH experiment.")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    scripts = ["01_check_data.py", "02_calibrate_rh.py", "03_validate_agbd_lobo.py",
               "04_export_figures_and_tables.py", "05_regional_comparison.py", "verify_reproduction.py"]
    if args.include_inverse_rh:
        scripts.insert(2, "experiments/rh_inverse_lobo.py")
    env = os.environ.copy()
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        env[name] = "1"
    for script in scripts:
        print(f"\nRunning {script}", flush=True)
        subprocess.run([sys.executable, str(root/script)], cwd=root, env=env, check=True)
    print("\nProcessed-data workflow completed. Read docs/RELEASE_STATUS.md before public release.")


if __name__ == "__main__":
    main()
