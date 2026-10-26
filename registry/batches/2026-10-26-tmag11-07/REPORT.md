# Two rejected multiyear M-dwarf signals

No novel planet claim is made. Both targets were searched with no primary or fixed alias surviving the declared initial screen. This negative result is limited to the recorded products, preprocessing, parameter range and gates.

| TIC | Sectors | Primary period, days | BLS SNR | Failed gates |
| --- | --- | --- | --- | --- |
| 12131189 | 11, 38, 91 | 2.9984997 | 6.81 | SNR, leave-one-epoch-out, per-transit coherence |
| 12881049 | 29, 69, 96 | 11.9350897 | 7.42 | Multi-sector transit baseline, SNR, leave-one-epoch-out, secondary eclipse null, per-transit coherence |

Selection used a complete count-checked clean-star query with Tmag 11.0 to 11.2, Teff 3300 to 4300 K, radius 0.20 to 0.70 solar radii and distance 5 to 55 pc. The fixed twelve-target metadata window was ranked before photometry using independent-reduction availability, multiyear coverage and a scheduling heuristic. The [plan](plan.json) digest is `eaf05cf477c519a24946124789e7ff807fffc9f32c13ba80b7e2ea13d1026171`. Known planets, TOIs, CTOIs, indexed external signals, previous scans and active deferrals were excluded. No complete literature novelty claim is made.

Unchanged upstream scientific scripts searched 1.1 to 15 days and evaluated fixed 2, 1/2, 3/2 and 2/3 aliases on the same final sector sets. All eight aliases failed. A failed secondary-eclipse or coherence gate does not establish a unique physical cause. The receipt archive retains coordinates, exact product URIs and hashes, cadence, stellar assumptions, all gate results, neighbor audits, catalog snapshots, dependencies and selection evidence. No independent-reduction recovery, pixel localization, injection sensitivity calibration or exclusion of all planets is claimed.

The original pipeline commit was `e63b80feff04ea3d19289cd389f30793f29ec3c1`. Private hosted search 40329536194 and exact replay 40330252198 succeeded. Their original compact receipts are released here so reproduction does not require access to the private repository. The archive hash is recorded in [verification.json](verification.json). No FITS products, caches or bulk telescope data are committed.

Public reproduction command: `python -B registry/batches/2026-10-26-tmag11-07/replay.py`.
Use Python 3.12, astropy 7.1.1, numpy 2.3.4, scipy 1.16.3, matplotlib 3.10.7, requests 2.32.5, OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. The command downloads public MAST products, requires exact file hashes and reproduces full primary and alias results against the retained official catalog snapshot. Preserved neighbor metadata queries are not independently rerun.

Upstream main was refreshed to `19e5e9e8ca72c8e584fd08a92276d285ca474888` before preparation. Neither TIC was registered. Fresh catalog receipts from metadata run 40338043659 were checked without collisions. [External checks](release-catalog-checks.json), [official source receipts](release-official-sources.json) and [known-index check](release-known-checks.json) retain hashes and scope.

Credit belongs to TESS, SPOC, QLP, MAST, NASA Exoplanet Archive, catalog authors, Thomas Marchand and Open Exoplanet Hunt Contributors. Original code and metadata follow the repository MIT license. Third-party data retain provider rights and acknowledgements. Codex assisted selection, execution, falsification review and preparation. No human scientific review, confirmation, external submission or invented affiliation is claimed.
