# Negative transit search: two high-SNR signals rejected by falsification

TIC 379167988 and TIC 243242290 were previously absent from the shared scan
registry and were rejected by the unchanged search pipeline. High BLS SNR
alone is insufficient. No unconfirmed candidate is promoted by this batch.

| TIC | Selected P (d) | BLS SNR | Gates passed | Failed gates |
| --- | --- | --- | --- | --- |
| 379167988 | 3.5207438 | 30.95 | 8/10 | Secondary anomaly, centroid shift |
| 243242290 | 11.7323661 | 23.20 | 6/10 | LOEO, odd/even, secondary anomaly, jackknife coherence |

## Selection, exclusions, and coverage

This is the next two-target batch from the fixed page-1 query with
`11.0 <= Tmag <= 11.2`, `3300 <= Teff <= 4300 K`,
`0.20 <= radius <= 0.70 Rsun`, `5 <= distance <= 55 pc`, and object type
`STAR`, requiring at least three sectors. The query was cached from the
preceding batch. Its raw rows, returned star order, arguments, and cache
hash are retained in `search/selection-page-1.json`.

The two preceding targets, TIC 115474342 and TIC 432694531 in
[PR #4](https://github.com/Th0rgal/open-exoplanet-hunt/pull/4), were excluded
before this run, alongside all 1,250 upstream registry targets. This PR adds
only the two new IDs and has no dependency on merging PR #4. Its leaderboard
count reflects the upstream registry plus these two rows.

Before searching, the official confirmed-planet, TOI, and CTOI exclusions
were refreshed using the unchanged `sync_catalogs.py`, preserving the
repository's catalog files. Neither target occurs in the refreshed
15,027-TIC index or the repository's indexed external objects. Source URLs,
response hashes, counts, timestamps, and membership results are under
`search/`. This exclusion scope does not cover every object in every
external survey.

TIC 379167988 has Tmag 11.0497, Teff 3850 K, radius 0.501545 Rsun, and
distance 47.2429 pc. Nine sectors were discovered: 18, 24, 25, 52, 58, 78,
79, 85, and 86. The upstream deterministic six-sector screen selected
18, 24, 52, 58, 85, and 86. Since it failed, upstream did not expand the
search to 25, 78, or 79. The registry entry is a negative result for this
recorded six-sector screen, not a claim that all nine sectors were searched.

TIC 243242290 has Tmag 11.1134, Teff 3621 K, radius 0.472980 Rsun, and
distance 42.2240 pc. Sectors 17, 42, 57, and 84 were discovered. Sector 42
failed the upstream loader with `Sector 42 has fewer than 50 valid cadences`.
The pipeline logged the failure and evaluated sectors 17, 57, and 84. The
failed sector was not evaluated or counted as a negative result.

The linked hosted log records the sector failure. Product-discovery
responses are saved in `search/mast-query-*.json`. Each successful evaluation
receipt records exact FITS URIs and SHA-256 hashes, stellar metadata, flux
columns, time coverage, finite centroid counts, all ten gate results, and
individual transit measurements. No FITS, downloaded data cache, or large
binary product is committed.

## Adversarial evidence

For TIC 379167988, the primary phase-0.5 dip is 491.7 +/- 44.6 ppm
(11.01 sigma), and the centroid shift is 0.0062 arcseconds (5.13 sigma).
The repeating signal has strong LOEO scores, but these two failed gates
prevent promotion. This evidence motivates eclipsing-binary, background,
and instrumental explanations without identifying a confirmed cause.

For TIC 243242290, odd transits have depth 167.5 ppm while even transits
have depth 3741.9 ppm, differing by 21.14 sigma. LOEO minimum and mean Z
are 0.60 and 1.28, and the jackknife SNR falls to 3.49 after removing the
strongest event. Its phase-0.5 anomaly is -2.78 sigma. These failures are
consistent with an unreliable ephemeris or outlier/alias explanation,
and do not support interpreting the best BLS peak as a planet candidate.

Rotational guards were applied to both targets and pass for the selected
periods. Every used sector has finite centroids. Instrumental and
single-sector explanations are probed through the inherited quality
selection, LOEO, centroid, and per-transit gates. These checks do not
exhaust every possible systematic effect.

A fixed supplemental audit evaluated `2P`, `P/2`, `3P/2`, and `2P/3` on
exactly the loaded FITS, checking their hashes before calculation. Each
alias uses the upstream plus or minus 0.3 percent targeted refinement and
the same ten gates. No alias was substituted into the registry result.

| TIC | 2P SNR / gates | P/2 SNR / gates | 3P/2 SNR / gates | 2P/3 SNR / gates |
| --- | --- | --- | --- | --- |
| 379167988 | 22.80 / 8 | 28.92 / 8 | 17.85 / 9 | 14.00 / 6 |
| 243242290 | 36.93 / 6 | 16.73 / 7 | 29.93 / 6 | 21.24 / 7 |

All eight alternatives fail at least one gate. In particular, the 3P/2
alias of TIC 379167988 still fails odd/even consistency, and every tested
alias of TIC 243242290 fails LOEO, odd/even consistency, and per-transit
coherence. The complete results remain in `audit/*-aliases.json`.

The existing analytical neighbor audit returns 200 and 11 cataloged
neighbors within 120 arcseconds, respectively. Under its aperture and
centroid-torque model, each is classified as excluded by flux or centroid
constraints. The full box queries completed in one page with 254 and 14
rows before angular filtering. The measured centroid failure on
TIC 379167988 remains a failure despite the neighbor model's classifications.
This model does not replace pixel difference imaging and cannot exclude
uncataloged or unresolved companions. Raw query responses and all
classifications are saved under `audit/`.

## Reproduction and limitations

Scientific scripts and thresholds are unchanged from upstream commit
`19e5e9e8ca72c8e584fd08a92276d285ca474888`. Exact dependency versions,
script hashes, run commits, and original search arguments are recorded.

- [Hosted search, verifier, and all ten adversarial controls](https://github.com/mashingaan/open-exoplanet-hunt/actions/runs/40285388770)
- [Hosted hash-checked aliases and neighbor queries](https://github.com/mashingaan/open-exoplanet-hunt/actions/runs/40285650550)
- [Executed search harness](https://github.com/mashingaan/open-exoplanet-hunt/blob/f27945b9020677dc0a2d5f5a63c6d226129522a8/.github/run_bounded_batch.py)
- [Executed supplemental harness](https://github.com/mashingaan/open-exoplanet-hunt/blob/f27945b9020677dc0a2d5f5a63c6d226129522a8/.github/run_supplemental_audit.py)

With the recorded dependencies installed in an isolated environment, from
the repository root run:

```sh
python -B registry/batches/2026-10-26-tmag11-02/replay.py
python -B scripts/verify_all.py
python -B -m unittest discover -s tests -v
```

The replay checks scientific script and FITS hashes, uses the saved stellar
metadata and original blind-search arguments, and compares the complete
primary evaluation against the receipt. It replays the actual searched
sectors rather than retrying the failed or unscreened sectors.

The inherited search selects one solution from five broad per-sector peaks.
Its nominal search bounds are 1.1 to 15 days, with broad seeds stopping at
12.5 days. LOEO refines near the full-data period using the full-data
duration and therefore is not a fully independent blind discovery test.
No thresholds, sector choices, detrending parameters, or aliases were
adjusted to make a signal pass. The results do not exclude all planets or
provide a detection-completeness estimate for either star.

Prepared with OpenAI Codex under the GitHub account mashingaan. No
independent astronomical confirmation or follow-up is claimed.
