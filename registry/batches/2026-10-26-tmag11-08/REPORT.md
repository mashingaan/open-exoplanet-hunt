# Bounded four-sector periodic-transit screen 08

Two previously unsearched targets completed the unchanged upstream BLS and
ten-gate screen. Neither primary signal nor any of the eight fixed aliases
survived. No private candidate hold was created.

| TIC | Sectors | Primary period, days | BLS SNR | Passed gates | Failed gates |
| --- | --- | --- | --- | --- | --- |
| 6749695 | 43, 44, 45, 71 | 12.2917548 | 6.89 | 5/10 | SNR, leave-one-epoch-out, odd/even, per-transit coherence, rotational-modulation guard |
| 9624343 | 25, 26, 52, 53 | 1.8533981 | 7.54 | 8/10 | SNR, leave-one-epoch-out |

Search [40332586047](https://github.com/mashingaan/tess-discovery/actions/runs/40332586047)
and exact replay [40333010767](https://github.com/mashingaan/tess-discovery/actions/runs/40333010767)
both succeeded at `e1dba40227f44754a26c2eb5f85eb9a9819ca04e`.
The replay compared both primary results and all fixed aliases and passed
deterministic verification and the test suite. Registry import checked the
full receipt manifest, file population, live terminal run states, exact source
commit and unchanged baseline before adding the two unique rows.

The sealed plan digest is
`81f347350fab3dfe20d9d01599dd8e4eb3c4b4f50da033eb2c025b7d21d13d5a`.
Metadata ranking and the product list were frozen before photometry. The
initial period range is 1.1 to 15 days, with alias multipliers 2, 0.5, 1.5
and 2/3. Full inputs and results remain in the compact receipt archive.

Failed rotation or odd/even gates do not establish a unique astrophysical
cause. These are rejected signals within the tested observations, products,
preprocessing and criteria, not exclusion of every possible planet. Sensitivity
calibration, independent reductions, pixel localization and complete literature
novelty review were not performed for these rejected signals. No novel planet
claim is made. No FITS products are committed.

## Public reproduction and release scope

The original pipeline commit, run receipts and compact archive are released here, so reproduction does not require private repository access.
Run `python -B registry/batches/2026-10-26-tmag11-08/replay.py` using Python 3.12, astropy 7.1.1, numpy 2.3.4, scipy 1.16.3, matplotlib 3.10.7 and requests 2.32.5, with OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1.
The replay checks the archive and scientific script hashes, downloads public MAST products, requires exact source-file hashes, and reproduces full primary and alias results against the preserved official catalog snapshot. Original neighbor metadata queries are retained rather than independently rerun.

The frozen [selection plan](plan.json) records a count-checked clean-star cohort with Tmag 11.0 to 11.2, Teff 3300 to 4300 K, radius 0.20 to 0.70 solar radii and distance 5 to 55 pc. Ranking used only metadata within the bounded twelve-target window, including independent-reduction availability, multiyear coverage and an uncalibrated CPU scheduling heuristic. Availability of another reduction is not independent signal recovery.

Upstream main was refreshed to `19e5e9e8ca72c8e584fd08a92276d285ca474888` before preparation and neither target was registered. Fresh official and external exclusions from metadata run 40338043659 contain neither target. Compact checks, timestamps and source hashes are retained in release-catalog-checks.json, release-known-checks.json and release-official-sources.json.

Credit remains with TESS, SPOC, QLP, MAST, NASA Exoplanet Archive, catalog authors, Thomas Marchand and Open Exoplanet Hunt Contributors. Original code and metadata follow the repository MIT license. Third-party data retain provider rights and acknowledgements. Codex assisted selection, execution, falsification review and preparation. No human scientific review, confirmation, external submission or invented affiliation is claimed.
