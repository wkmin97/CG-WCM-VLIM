# GEDI waveform correction and forest biomass estimation

Analysis code and processed inputs accompanying the manuscript:

**Restoring forest structural consistency for GEDI footprint-level aboveground biomass estimation through physically based waveform correction**

Wankun Min, Wenli Huang, Liwei Yao, and Yumin Chen.

The study examines terrain effects on simulated GEDI relative-height (RH) metrics and estimates aboveground biomass density (AGBD) using vertical profile information from 12 UAV lidar blocks.

**Release status: repository preparation draft.** The runnable components reproduce the supplied processed-data calculations. The complete raw-data-to-manuscript workflow is not yet available. [Release status](docs/RELEASE_STATUS.md) identifies the missing inputs and the differences between the manuscript and the supplied implementation.

## Data

The processed data are included; no separate download is required for the numbered Python workflow.

| File | Contents | Rows |
| --- | --- | ---: |
| [data_pred_all_data_v2.csv](data/processed/data_pred_all_data_v2.csv) | Analysis cohort, block IDs, reference RH/AGBD, simulated RH, and stored predictions | 9,790 |
| [GEDI_SIM_metrics.csv](data/processed/GEDI_SIM_metrics.csv) | Broader simulator/reference metrics table used for provenance checks | 17,707 |
| [GEDI_2019_with_Predicted_AGBD.csv](data/regional/GEDI_2019_with_Predicted_AGBD.csv) | Supplied regional GEDI and comparison-product columns | 7,832 |

See the [data dictionary](data/data_dictionary.csv), [input checksums](data/manifest.json), and [data notes](data/README.md). Original values and column names are preserved. The AGBD reference is model-derived and field-calibrated, rather than a direct field measurement at every footprint.

Raw UAV point clouds, raster products, field measurements, and the complete reference-biomass generation workflow are not included. Their required roles are listed in [data/raw/README.md](data/raw/README.md).

## Environment

The Python workflow was run with Python 3.12 and the exact package versions in [requirements.txt](requirements.txt). R is required only for the archived upstream source scripts, which are not part of the automated run.

With conda, on Windows, macOS, or Linux:

```bash
conda env create -f environment.yml
conda activate gedi-gc-wcm-vlim
```

Alternatively, install the dependencies in an existing Python 3.12 environment:

```bash
python -m pip install -r requirements.txt
```

## Run the workflow

Open a terminal in this repository folder. Run all five numbered steps and the verification:

```bash
python run_all.py
```

The same steps can be run individually:

```bash
python 01_check_data.py
python 02_calibrate_rh.py
python 03_validate_agbd_lobo.py
python 04_export_figures_and_tables.py
python 05_regional_comparison.py
python verify_reproduction.py
```

| Step | Purpose | Main output |
| --- | --- | --- |
| `01_check_data.py` | Verify file hashes, sample keys, block membership, and input coverage; recompute metrics from stored biomass predictions | `results/data_audit/` |
| `02_calibrate_rh.py` | Reproduce the original four-coefficient, reference-assisted RH calibration | `results/rh_calibration/` |
| `03_validate_agbd_lobo.py` | Refit the layer/VLIM model while holding out each complete UAV block, using the supplied original and corrected RH inputs | `results/agbd_conditional_LOBO/` |
| `04_export_figures_and_tables.py` | Export block scores, compare them with manuscript Table 3, and draw residual/scatter figures | `results/figures/` and Table 3 comparison CSVs |
| `05_regional_comparison.py` | Summarize and plot the regional CSV columns as stored | `results/regional_diagnostics/` |
| `verify_reproduction.py` | Independently recalculate scores and check prediction/sample alignment | `results/verification.json` |

Existing outputs are included. Rerunning overwrites files in the corresponding `results/` directories. Figures are exported as PNG, SVG, and PDF.

The packaged run passed the technical checks and exactly reproduced the supplied layer-model pooled RMSE values (11.651933 for original RH and 11.225270 Mg/ha for stored corrected RH). These are conditional LOBO scores. The block results do **not** reproduce manuscript Table 3 at its reported precision; the differences are exported for review. See [verification results](results/verification.json) and the [run report](docs/VERIFICATION_REPORT.md).

## Interpretation of the reproduced calculations

The supplied corrected `Pred_RH*` values are reproduced by the original four-coefficient regression, which uses UAV reference RH in calculating the correction. Step 2 reproduces this calibration; it does not establish independent correction accuracy and differs from the manuscript's written two-term GC-WCM equation.

Step 3 excludes the test block when fitting AGBD parameters, but uses the existing `Pred_RH*` predictors. It therefore provides **conditional AGBD-model LOBO validation**, not independent validation of the complete RH-correction-to-AGBD pipeline. Step 3 uses the supplied corrected columns directly; it does not consume predictions from a separate RH validation experiment. The original and revised validation designs must remain distinguishable.

Step 5 does not apply a new correction, infer an official GEDI L4A product from a similarly named column, or apply undocumented changes of units or scale. Its output is a regional data diagnostic, not a verified reproduction of manuscript Figure 14.

The supplied sources do not establish the manuscript's claimed 10-fold out-of-fold results. Stored prediction metrics and manuscript Table 2 are retained separately for reconciliation. See [calibration and validation](docs/CALIBRATION_AND_VALIDATION.md).

## Separate RH validation experiment

The supplied revision script also implements an inverse-calibration RH experiment that predicts held-out RH without using test-block reference RH:

```bash
python experiments/rh_inverse_lobo.py
```

Or include it in the automated run:

```bash
python run_all.py --include-inverse-rh
```

This is a different prediction rule from the stored reference-assisted calibration. Its results are kept in `results/rh_inverse_LOBO/` and are not substituted into the conditional AGBD experiment. Its performance may be worse than raw RH; the script preserves those results.

## Repository contents

| Location | Contents |
| --- | --- |
| `src/` | Shared utilities and the supplied revision-model implementations |
| `config/block_mapping.csv` | Mapping of flight-block identifiers to manuscript block numbers, matched using Table 1 terrain ranges |
| `data/processed/`, `data/regional/` | Unmodified processed input tables |
| `data/reported/` | Manuscript Tables 2/3 and supplied numerical revision results |
| `original_code/` | Selected original R/Python source files, preserved unchanged for provenance |
| `docs/source_code_inventory.csv` | Inclusion/exclusion decisions for the supplied source archives |
| `docs/MANUSCRIPT_OUTPUT_MAP.md` | Current manuscript figure/table mapping and reproduction coverage |
| `docs/UPSTREAM_WORKFLOW.md` | Raw-data processing stages, available source, and missing dependencies |
| `docs/GITHUB_UPLOAD_GUIDE.md` | Upload instructions and public-release preparation |

Exploratory multitask/deep-learning scripts and checkpoints were excluded from the main workflow because they were not identified as final manuscript implementations. The inventory records these decisions for author review. Archived source files retain their original local paths and interactive sections; use the numbered scripts for the automated processed-data run.

## Citation and licensing

Draft citation metadata are supplied in [CITATION.cff](CITATION.cff). Add the public repository URL, final article details, and archived-release DOI when available.

No distribution license has been selected in this preparation draft. Confirm the code license and data redistribution terms before making the final public release; see [LICENSE_STATUS.md](LICENSE_STATUS.md).
