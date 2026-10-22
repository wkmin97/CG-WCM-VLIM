# Data included in this draft

Three processed input tables are included unchanged. The numerical contents were recovered from the supplied manuscript/revision materials; data are not simulated or fabricated for this package.

| Directory | Role |
| --- | --- |
| `processed/` | Local UAV/simulator metrics and the fixed 9,790-footprint analysis cohort |
| `regional/` | Supplied regional GEDI/product extraction and stored predictions |
| `reported/` | Values transcribed programmatically from the latest manuscript tables and the supplied revision result tables; these are comparison targets, not newly validated outputs |

`manifest.json` fingerprints the nine input/comparison CSV files. `data_dictionary.csv` describes all 500 column occurrences across the three main tables. Unresolved abbreviations or units are marked as unconfirmed rather than guessed.

The analysis cohort contains 12 blocks. Its key is `site` plus `ID`; IDs must not be treated as globally unique without the block ID. Reference RH and AGBD must remain separate from simulated/predicted metrics.

The broader metrics table is not a complete raw parent of the analysis cohort: 49 cohort keys are absent. Step 1 exports those keys. No samples are silently removed or substituted.

Stored `Pred_RH*` values use UAV reference RH in their calibration. Stored `Pred_*AGBD*` columns are not evidence of independent or 10-fold validation without the generating split/prediction records.

Regional metadata require review: exact footprint geography/year, product versions, acquisition dates, units, flag filtering, and the similarly named `CFATD_AGBD`/`CFATD_AG_1` columns. Step 5 preserves the CSV values and performs no additional quality filtering or rescaling.
