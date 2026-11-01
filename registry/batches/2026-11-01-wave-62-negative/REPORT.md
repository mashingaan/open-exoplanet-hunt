# Negative TESS screen: 62 completed targets, 2 unresolved

## Results

The frozen wave planned 64 previously unscanned targets in eight partitions. It completed 62 primary screens and 248 preregistered alias evaluations. All 62 primary results and all 248 aliases were rejected by the unchanged falsification gates. No candidate survived. These results apply only to the recorded light curves and search configuration and do not show that the targets have no planets.

Partition 0 stopped with six completed targets. TIC 188567493 could not pass the unchanged valid-cadence requirement because Sector 51 had fewer than 50 valid cadences. The worker stopped before reaching TIC 148734869. Neither target is included in `registry/scanned_targets.jsonl` or the 62-result manifest. Both remain unresolved.

## Frozen search and falsification scope

The metadata-only plan used a fixed 96-probe window and selected 64 unique targets for eight groups of eight. The screen required three to twelve discovered sectors and searched periods from 1.1 to 15 days. The fixed alias multipliers were 2, 0.5, 1.5, and 2/3. The exact search-time scientific scripts and hashes are included under `scientific-source/scripts/`. No target, sector, epoch, threshold, or alias was chosen after inspecting results. Failed workers resumed their same frozen partition.

All primary rows have verdict `REJECTED_BY_FALSIFICATION`. Primary failed-gate counts are: gate 3, 48; gate 4, 60; gate 5, 9; gate 6, 15; gate 7, 9; gate 8, 7; gate 9, 43; gate 10, 22. Counts overlap because targets can fail multiple gates. All 248 aliases fail at least one gate. Alias failed-gate counts are: gate 2, 22; gate 3, 214; gate 4, 247; gate 5, 56; gate 6, 40; gate 7, 30; gate 8, 19; gate 9, 229; gate 10, 112.

The unchanged engine tests multi-sector baselines, combined BLS significance, leave-one-epoch-out repeatability, odd/even consistency, secondary eclipses, physical plausibility, centroid stability, per-transit coherence, and rotational modulation. Passing centroid or secondary-eclipse measurements are not treated as proof against blends or eclipsing binaries.

## Duplicate checks

The search-time refresh recorded no indexed catalog collision among the 64 planned TICs. Before submission, upstream main and every open scan PR were checked again. No TIC overlap was found. TIC matching does not replace a literature search or coordinate-based novelty review.

## Reproduction and limits

Run a frozen partition with `python -B registry/batches/2026-11-01-wave-62-negative/replay.py --partition 0`, replacing `0` with an index from 0 through 7. The replay loads the bundled, SHA-256-pinned `vet_candidate.py`, downloads exact public MAST products, checks their hashes, and reproduces every primary and fixed-alias result. Dependencies are pinned in `replay-requirements.txt`. GitHub Actions runs repository verification and adversarial tests once, then replays all eight partitions on hosted runners.

Neighbor-query receipts are retained but the independent MAST neighbor query is not repeated by replay.py. Centroid and secondary-eclipse checks use recorded photometric gates and are not a pixel-level background-blend or eclipsing-binary audit. This batch includes no separate momentum-dump timing audit or injection-recovery sensitivity calibration. No independent human scientific review or confirmation is claimed.

Credit belongs to TESS, MAST, SPOC, QLP, catalog providers, and Open Exoplanet Hunt contributors. Codex assisted with execution and preparation. No discovery or confirmation is claimed.
