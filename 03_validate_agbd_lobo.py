#!/usr/bin/env python3
"""Re-run the supplied conditional AGBD-model LOBO experiment.

The default layer model corresponds to the VLIM rows discussed in Table 3.
Pred_RH is reference-assisted. These scores test AGBD fitting conditional on
these RH inputs; they do not validate the full correction-to-biomass pipeline.
"""
import argparse
from src.common import DATA, RESULTS
from src.lobo_agbd import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(DATA))
    parser.add_argument("--outdir", default=str(RESULTS/"agbd_conditional_LOBO"))
    parser.add_argument("--models", default="layer", help="layer,integral,power,stepwise; extra families are source-code experiments, not verified Table 2 equivalents")
    parser.add_argument("--stepwise-max-features", type=int, default=4)
    args = parser.parse_args()
    args.block_col, args.target_col, args.expected_blocks = "site", "ALS_AGB", 12
    run(args)


if __name__ == "__main__":
    main()
