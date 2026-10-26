# Two rejected multiyear M-dwarf signals

No novel planet claim is made. Both targets were searched with no primary or
fixed alias surviving the declared initial screen. This negative result is
limited to the recorded products, preprocessing, parameter range and gates.

| TIC | Sectors | Primary period, days | BLS SNR | Failed gates |
| --- | --- | --- | --- | --- |
| 10835890 | 3, 30, 97 | 2.0536355 | 7.13 | SNR, leave-one-epoch-out, per-transit coherence |
| 11179465 | 7, 34, 61 | 6.4731055 | 5.94 | SNR, leave-one-epoch-out, per-transit coherence, rotation guard |

Selection used a complete count-checked clean-star query with Tmag 11.0 to
11.2, Teff 3300 to 4300 K, radius 0.20 to 0.70 solar radii and distance 5 to
55 pc. The fixed twelve-target metadata window was ranked before photometry
using independent-reduction availability, multiyear coverage and a scheduling
heuristic. The [plan](plan.json) digest is
`6edecfed83feb9927066f9b3f26cfb36455fdd8059fad82b16c6a15282494fc2`.
Known planets, TOIs, CTOIs, indexed external signals, completed scans and active
deferrals were excluded. No complete literature novelty claim is made.

Unchanged upstream scientific scripts searched 1.1 to 15 days and evaluated
the fixed 2, 1/2, 3/2 and 2/3 aliases on the same final sector sets. All eight
aliases failed. TIC 10835890 used SPOC 120-second products in all three sectors,
retaining 12911, 15960 and 29082 samples. TIC 11179465 used TESS-SPOC 1800-second,
QLP 600-second and QLP 200-second products, retaining 1079, 2973 and 9059 samples.
All planned sectors were loaded without skips. The receipt archive retains
coordinates, exact product URIs and hashes, stellar assumptions, passed/failed
gates, neighbor audits, catalog snapshots, dependencies and selection evidence.

TIC 11179465 had a high catalog contamination estimate before photometry.
It remains unresolved. The failed rotation guard does not establish a unique
physical explanation. No pixel localization, injection sensitivity calibration,
independent-reduction recovery or exclusion of all planets is claimed.
A stitched mixture of reductions is not independent confirmation.

The original pipeline commit was `3cad450980c3c328154b79d609b10820a9e378a8`.
Private hosted search 40327617744 and exact replay 40328148902 succeeded.
Their compact original receipts are released here so reproducing the result
does not require access to the private orchestration repository.
The archive SHA-256 is recorded in [verification.json](verification.json).
It contains no FITS products, caches or bulk telescope data.

Public reproduction command:
`python -B registry/batches/2026-10-26-tmag11-06/replay.py`.
Use Python 3.12, astropy 7.1.1, numpy 2.3.4, scipy 1.16.3, matplotlib 3.10.7,
requests 2.32.5, OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1.
The command downloads public MAST files, requires exact source-file hashes and
reproduces complete primary and alias evaluations against the released catalog
snapshot. The preserved neighbor metadata query is not independently rerun.

Before release preparation, upstream main was refreshed to
`19e5e9e8ca72c8e584fd08a92276d285ca474888`. Neither TIC was in its registry.
Fresh catalog receipts from metadata run 40332102913 were also checked for
both targets without collisions. [External checks](release-catalog-checks.json),
[official source receipts](release-official-sources.json) and
[known-index check](release-known-checks.json) preserve their hashes and scope.

Credit belongs to TESS, SPOC, QLP, MAST, NASA Exoplanet Archive, indexed catalog
authors, Thomas Marchand and Open Exoplanet Hunt Contributors. Original code and
metadata follow the repository MIT license. Third-party data retain provider
rights and acknowledgements. Codex assisted selection, execution, falsification
review and preparation. No human scientific review, confirmation, external
submission or invented affiliation is claimed.
