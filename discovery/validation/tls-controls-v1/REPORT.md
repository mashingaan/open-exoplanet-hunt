# TLS synthetic controls, exact replay failure retained

The runs in this development history were executed in a private orchestration
repository, so their Actions links may require access. Complete compact receipts
are included in this public contribution without astronomical target searches
or candidate data. Its own fresh public validation is recorded separately
in `public-verification.json`. A development run's commit ID
is provenance for that execution, not a claim that this public repository
contains that private commit.

Both hosted executions used commit
`322682ee6e97fe3f6ef1ed0491a3020c916a7c08` and the unchanged specification
SHA-256 `e4c5b7a90f35d493bde2f81c71ad3f3f8d5128a6589b8e73a42a496b73215d4f`.
Run [40312711422](https://github.com/mashingaan/tess-discovery/actions/runs/40312711422)
passed the six controls. Its independent repeat,
[40313351384](https://github.com/mashingaan/tess-discovery/actions/runs/40313351384),
passed the same controls but failed the required exact scientific receipt
comparison. This original repeat remains a failure. Subsequent bounded
environment validation is recorded below. Real-target deployment remains
unauthorized.

| Fixed control | Raw TLS period, days | Raw BLS period, days | Required outcome observed in both runs |
| --- | ---: | ---: | --- |
| Positive | 3.182862978 | 3.183749618 | Both recover injected 3.1828 days, primary passes all ten gates |
| White noise | 2.429843695 | 6.845188285 | Primary and all four aliases rejected |
| Rotation | 7.994220190 | 3.187286662 | Rotation guard rejects primary, no alias survives |
| Secondary eclipse | 3.182862978 | 1.591535651 | Secondary-null gate rejects primary, no alias survives |
| Missing third-sector signal | 3.182862978 | 3.183749618 | Leave-one-epoch-out gate rejects primary, no alias survives |
| Centroid blend | 3.182862978 | 3.183749618 | Centroid gate rejects primary, no alias survives |

The positive reference and fused backends have identical period grids and
spectra within the originally specified tolerances in both executions. This
within-run tolerance is distinct from the exact between-run replay requirement.
The rotation disagreement and the BLS half-period for the secondary control
remain recorded. High detector statistics alone do not justify acceptance.

## Reproduction discrepancy

Seventeen JSON fields differ between the runs. They consist only of TLS
spectral hashes and SDE or raw SDE values in the last floating-point digits.
Input hashes, all BLS results, recovered periods, all downstream gate receipts,
alias results and control checks are exactly identical. Dependencies, Python
and platform strings, upstream script hashes and warnings also agree.

This is evidence of numerical nonidentity, not proof of its cause. Neither
original receipt captured the actual CPU feature dispatch or raw chi-square
arrays. CPU or reduction-order differences are hypotheses, not established
facts. No tolerance is introduced to relabel the failed repeat as passing.

The three retained TLS warnings include two period-grid stellar-default
fallback warnings and a warning that two of five transit epochs lack data.
The original warning collection does not identify their individual control
or call site. They must not be attributed to a specific case without further
evidence. In the pinned TLS source, the period-grid fallback preserves the
requested period bounds. The adapter additionally checks those bounds.

## Retained evidence and scope

[verification.json](verification.json) records live terminal hosted states,
source hashes and every field difference. The 54686-byte `receipts.tar.gz`
preserves complete scientific receipts, dependencies, commands, warnings and
original manifests from both runs, together with the replay failure receipt.
All 22 source and 23 replay manifested files were downloaded and hash-checked
before archival. Original spectral arrays and diagnostic plots are retained
in the hosted artifacts for their 30-day retention period. Their hashes remain
in the archived original manifests, even though those files are omitted from
this compact archive. The six-row diagnostic plot was inspected.

These are synthetic controls with matched fixed noise, not astronomical
candidates, completeness estimates or false-alarm calibration. They add no
target to the astronomical coverage registry. Real MAST recovery, independent
reductions, long baselines and sparse or heterogeneous sampling remain untested.

## Prespecified environment investigation

The next pair of runs retains the original scientific specification, all
thresholds, seeds, aliases and the exact replay criterion. Only numerical
execution and diagnostics change: generic Numba CPU target, NumPy baseline
dispatch, a separate JIT cache namespace, recorded CPU/runtime configuration
and saved raw TLS chi-square and undetrended spectral diagnostics. This tests
whether a constrained execution environment restores exact replay. It does
not retrospectively repair or erase either earlier run.

## Generic execution pair passed without changing scientific criteria

Runs [40314050728](https://github.com/mashingaan/tess-discovery/actions/runs/40314050728)
and [40314272648](https://github.com/mashingaan/tess-discovery/actions/runs/40314272648)
both completed successfully at commit
`deff0b5ae7a94503e456952ec89930dbad69ce3c`. Each passed the verifier, all 36
deterministic tests and all six fixed scientific controls. The second execution
used `--compare` against the first. Both retained the original specification
digest, input hashes and acceptance criteria.

The complete `science.json` files match byte for byte. All six compressed
spectral containers also match byte for byte, including the newly saved raw
TLS chi-square and raw power arrays. All 23 manifested files in each artifact
were downloaded and verified. Recorded CPU models differ: AMD EPYC 9V45 in
the first run and AMD EPYC 9V74 in the repeat. Both used generic Numba and
NumPy's SSE/SSE2/SSE3 baseline with no enabled dispatched features.

[generic-verification.json](generic-verification.json) and the 58522-byte
`generic-receipts.tar.gz` preserve terminal states, hashes, complete scientific
receipts, runtime configuration, dependencies, warnings and original manifests.
They supplement the earlier failure archive rather than replacing it.

This pair establishes exact replay for the tested synthetic inputs and recorded
environment on these two hosts. It does not establish the original failure's
mechanism, universal cross-platform determinism or astronomical performance.
The same three warnings remain. Real-product recovery and independent
photometric reductions are still required before operational TLS searches.

## Public pair and reusable result record

Public runs [40315351988](https://github.com/mashingaan/open-exoplanet-hunt/actions/runs/40315351988)
and [40315612271](https://github.com/mashingaan/open-exoplanet-hunt/actions/runs/40315612271)
passed at `99ae48ad78a2ddf3d4566664e393214fd4d907e4`, retained by fork tag
`tls-controls-v1-experiment`. Each passed the verifier, all 12 upstream-plus-TLS
tests and all six scientific controls. The second run's exact comparison passed.
All 23 manifested files per run were hash-checked. Complete scientific receipts
and all six spectral containers match byte for byte, and the public scientific
receipt also matches the validated generic development receipt exactly.

The stable result ID is `OEH-CONTROL-TLS-V1-20261026`. Its
[public-verification.json](public-verification.json) records status, reason
codes, checksums, terminal public run states, limitations and source equivalence
to the generic development implementation. The 58540-byte
`public-receipts.tar.gz` permanently retains the public scientific, runtime,
warning, environment and original manifest receipts. Numerical arrays and plots
remain hosted for 30 days and can be regenerated from the tagged source.

Only synthetic controls were evaluated. No new planet claim or astronomical
search coverage is made, and no candidate dossier or private observation was
released. Known objects and new observational products require their own future
validation. Original contribution code and metadata use the existing MIT
license, with external package authorship, license and citation requirements
preserved in the method note.
