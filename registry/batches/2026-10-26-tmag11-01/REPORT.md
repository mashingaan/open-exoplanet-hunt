# Negative transit search: two previously unscanned M dwarfs

Both targets were rejected by the unchanged search and falsification pipeline.
The primary signals fail gates 3, 4, 9, and 10. None of the eight prescribed
rational-alias evaluations passes all ten gates. These are negative search
results within the recorded data and search procedure, not exclusions of all
possible planets around either star.

| TIC | Tmag | Teff (K) | Distance (pc) | Selected P (d) | BLS SNR | LOEO min / mean Z | Jackknife SNR | OOT amplitude / dip |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 115474342 | 11.1123 | 3391 | 24.1127 | 1.3476502 | 6.51 | 0.28 / 1.62 | 6.18 | 0.362 |
| 432694531 | 11.1107 | 3576 | 37.4485 | 2.5374132 | 6.98 | 1.66 / 2.68 | 5.16 | 0.445 |

## Fixed selection and data

The first two eligible stars in page 1 of the recorded MAST response were
selected with `max_targets=2`, `min_sectors=3`, `11.0 <= Tmag <= 11.2`,
`3300 <= Teff <= 4300 K`, `0.20 <= radius <= 0.70 Rsun`, and
`5 <= distance <= 55 pc`, with TIC object type `STAR`. The original 1,250-row
registry had no stars in this magnitude interval. Both IDs were absent from
that registry, the repository exclusion index, and the refreshed official
catalogs plus the repository's indexed external objects.

The official exclusions were refreshed before searching, at
`2026-10-26T07:14:17Z`, using `sync_catalogs.py` without changing the repository
catalog files. The refreshed index contained 15,027 unique TIC IDs.
`search/catalog-source-receipts.json` records public URLs and response hashes.
`search/exclusion-check.json` records the two negative membership checks.
The external exclusion scope is the repository's indexed objects, not an
exhaustive claim about every published survey.

Both targets have three discovered sectors, 17, 57, and 84, and all three were
loaded. No six-sector truncation or skipped sector was needed. The MAST
selection response and product-discovery responses are retained in
`search/mast-query-*.json`. The evaluation receipts include stellar metadata,
exact FITS URIs, SHA-256 hashes, selected flux columns, cadences, time coverage,
all ten gate results, and individual transit measurements. No FITS is in git.

## Falsification evidence

The primary odd/even depth differences are 2.27 and 1.58 sigma. The phase-0.5
anomalies are -0.39 and -2.29 sigma. The measured centroid shifts are 1.31 and
0.56 sigma, with finite centroids in every loaded sector. These passing checks
do not overcome the failed SNR, temporal prediction, jackknife, and rotational
checks, and do not establish a specific astrophysical false-positive class.

An additional fixed audit tested `2P`, `P/2`, `3P/2`, and `2P/3` for each target
on the same six FITS products, verifying their hashes first. Each alias uses
the existing targeted-period refinement of plus or minus 0.3 percent and the
same ten gates. No alias replaces the original registry ephemeris.

| TIC | 2P SNR | P/2 SNR | 3P/2 SNR | 2P/3 SNR | All ten gates passed |
| --- | --- | --- | --- | --- | --- |
| 115474342 | 6.31 | 5.49 | 5.23 | 5.21 | None |
| 432694531 | 6.25 | 5.63 | 5.22 | 5.72 | None |

The existing TIC-neighbor audit queried within 120 arcseconds. Under its
aperture-coupling and centroid-torque model, it classified all 41 and 24
cataloged neighbors as ruled out by flux or centroid constraints. This
analytical model is not pixel difference imaging or high-resolution imaging,
and cannot exclude uncataloged or unresolved companions. Raw query responses
and complete neighbor results are retained under `audit/`.

## Reproduction and validation

The scientific scripts are unchanged from upstream commit
`19e5e9e8ca72c8e584fd08a92276d285ca474888`. The exact dependency versions and
script hashes are in `search/provenance.json`.

- [Hosted search and all 10 adversarial controls](https://github.com/mashingaan/open-exoplanet-hunt/actions/runs/40284467234)
- [Hosted hash-checked alias and neighbor audit](https://github.com/mashingaan/open-exoplanet-hunt/actions/runs/40284747762)
- [Search harness at its executed commit](https://github.com/mashingaan/open-exoplanet-hunt/blob/b78437bccab9d068c0fedbfa748944f3bbdb0286/.github/run_bounded_batch.py)
- [Supplemental harness](https://github.com/mashingaan/open-exoplanet-hunt/blob/8245a0b/.github/run_supplemental_audit.py)

From the repository root, install the recorded dependencies in an isolated
environment and replay the primary searches from the saved receipts with:

```sh
python -B registry/batches/2026-10-26-tmag11-01/replay.py
python -B scripts/verify_all.py
python -B -m unittest discover -s tests -v
```

The replay downloads public MAST products into the ignored cache and refuses
changed FITS or scientific scripts. It uses the saved stellar parameters and
the original blind-search arguments. The supplemental harness records the
fixed alias evaluation procedure and the original neighbor queries remain
available even if the live catalog changes.

The inherited two-stage search selects one result from five broad per-sector
peaks. Its nominal bounds are 1.1 to 15 days, but broad seed periods stop at
12.5 days. LOEO refines periods near the full-data solution and uses the
full-data duration, so it is not a fully independent blind discovery test.
No completeness or planet-occurrence-rate inference follows from this batch.
No thresholds, detrending settings, sectors, or aliases were tuned to rescue
a signal.

Contribution prepared with OpenAI Codex under the GitHub account mashingaan.
No independent astronomical confirmation or follow-up is claimed.
