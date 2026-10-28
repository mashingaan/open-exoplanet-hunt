# Replayed negative prefix of an interrupted TESS wave

Four targets completed the fixed periodic TESS initial screen before their
coordinated wave stopped on infrastructure errors. All four primary signals
and their 16 fixed alias evaluations were rejected by the recorded tests.
These results apply only to the recorded data, preprocessing, search bounds,
and criteria. They do not exclude planets and make no novel planet claim.

The original search used unchanged upstream scientific scripts. Three private
hosted replays independently reloaded public MAST products, checked their
hashes, and reproduced the completed evaluations and aliases. The neighboring
source audits are retained but were not independently repeated. Original
private source artifacts contain unfinished targets and are not released here.

The original engine can assign passing centroid and phase 0.5 secondary tests
when measurements are insufficient. These recorded fields must not be taken
as proof those scenarios were excluded. Each primary failed other documented
gates. LOEO uses a period selected on the full data and is a repeatability
check rather than a fully blind search. Injection sensitivity is not calibrated.

`target-catalog.json.gz` contains only entries for these four targets from the
search-time index. Its hash and the full-index hash are recorded in the
manifest and provenance. It reproduces the target-specific catalog lookups, not
the complete source catalog. A later catalog check found no indexed collision
for these targets in the listed sources. That check is not exhaustive
literature or coordinate-based novelty review.

To reproduce, use the dependency versions in `provenance.json` and run
`python -B registry/batches/2026-10-28-wave-prefix/replay.py` from a checkout
whose scientific script hashes match. Public FITS are downloaded from MAST.

Credit belongs to TESS, the relevant photometry pipelines, MAST, catalog
providers, and Open Exoplanet Hunt contributors. Codex assisted analysis and
preparation. No human scientific review or independent confirmation is claimed.
