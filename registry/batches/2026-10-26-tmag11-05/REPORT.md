# Two rejected Tmag 11 M-dwarf signals

No novel planet claim is made. Both targets were searched with no signal
surviving the declared initial screen and its four fixed aliases. These null
results apply only to the recorded products, preprocessing and parameter range.

| TIC | Sectors | Primary period, days | BLS SNR | Failed gates |
| --- | --- | --- | --- | --- |
| 6679510 | 7, 34, 88 | 11.3257116 | 5.63 | SNR, leave-one-epoch-out, secondary null, per-transit coherence |
| 7130810 | 7, 34, 88 | 3.5297176 | 6.03 | SNR, leave-one-epoch-out, per-transit coherence |

Selection used the complete count-checked clean-star query with Tmag 11.0 to
11.2, Teff 3300 to 4300 K, radius 0.20 to 0.70 solar radii and distance 5 to
55 pc. The fixed twelve-target metadata window was ranked before photometry,
using independent-reduction availability, multiyear coverage and a scheduling
heuristic. The two-target [plan](plan.json) digest is
`a403a13b90f83fea8af490d68e75ae8c9512707c009ac840d19013a072b5d77a`.
Known confirmed planets, TOIs, CTOIs, indexed external signals, completed scans
and active deferrals were excluded. This is not a complete literature novelty
search and no novelty assertion is necessary for these rejected signals.

The unmodified upstream scientific scripts searched 1.1 to 15 days. The primary
and fixed 2, 1/2, 3/2 and 2/3 aliases were evaluated on the same loaded sectors.
All eight aliases failed. Both stars used all three planned sectors without
skips. TIC 6679510 used SPOC 120 s, TESS-SPOC 600 s and QLP 200 s products,
retaining 16337, 3350 and 9679 points. TIC 7130810 used TESS-SPOC 1800 s,
TESS-SPOC 600 s and QLP 200 s products, retaining 1085, 3346 and 9680 points.
The receipt archive contains exact product URIs, hashes, stellar coordinates,
full gate results, neighbor audits, exclusions, dependencies and selection evidence.

The negative secondary statistic for TIC 6679510 does not identify an eclipsing
binary. Neighbor audits leave three and seven potential contaminants unresolved
under the approximate PSF model. No pixel localization, independent reduction
recovery, injection sensitivity calibration or physical planet exclusion is
claimed. A stitched mix of available reductions is not independent confirmation.

The original search commit was `524837b0d2522eecb917cdc4ea99ce26784bad93`.
Original hosted search 40309387639 and exact replay 40309858069 succeeded.
Those runs belong to a private orchestration repository. Their original
receipts are released here so reproduction does not require repository access.
The archive hash is recorded in [verification.json](verification.json).
It contains compact metadata, not FITS products, caches or telescope datasets.

Public reproduction command:
`python -B registry/batches/2026-10-26-tmag11-05/replay.py`.
Use Python 3.12, astropy 7.1.1, numpy 2.3.4, scipy 1.16.3, matplotlib 3.10.7,
requests 2.32.5, OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1.
It downloads original public MAST products, checks each SHA-256, and requires
exact primary and alias evaluations against the released catalog snapshot.
Neighbor metadata queries are preserved but not independently replayed by this command.
Public hosted [reproduction 40328329643](https://github.com/mashingaan/open-exoplanet-hunt/actions/runs/40328329643)
passed at contribution commit `ef308e3`. Its checkout commit was verified in
the log because the workflow itself runs from the fork's runner branch.
Both primary evaluations and all eight fixed aliases reproduced exactly.
Repository verification, all ten adversarial tests and the clean-tree check passed.

Before preparing this release, upstream main was refreshed to
`19e5e9e8ca72c8e584fd08a92276d285ca474888` and neither target was present in its
registry. Fresh official and enabled external catalog exclusions were checked
again on 2026-10-26. External response hashes and retrieval times are retained
in [release-catalog-checks.json](release-catalog-checks.json).

Credit belongs to TESS, SPOC, QLP, MAST, NASA Exoplanet Archive, the indexed
catalog authors, and Thomas Marchand and Open Exoplanet Hunt Contributors.
Original code and metadata follow the repository MIT license. Archive data and
third-party catalog rights remain with their original providers. Codex assisted
selection, execution, falsification review and preparation. No human scientific
review, confirmation, external submission or invented affiliation is claimed.
