# Preregistered TLS component evaluation

This contribution contains only the generic detector, synthetic controls and
their reproducibility evidence. No astronomical scan or private candidate is
included. Existing `scripts/` and the astronomical registry remain unchanged.
Development history below includes private Actions links. The complete compact
receipts are released here for inspection. A fresh public hosted pair has
also completed at this contribution's own frozen experiment commit, with
[machine-readable public verification](../validation/tls-controls-v1/public-verification.json).
The ordinary upstream suite has ten tests, with two input-contract tests added
here. The 36-test counts below describe the separate development workspace.
The cases were fixed before the original development execution. Public runs
reproduce those already known control outcomes, rather than constituting a
new blinded experiment.

Original contribution code and metadata use the repository's MIT license.
External packages retain their own licenses and citation requirements. No
third-party package source or telescope data is redistributed here.

## Public reproduction evidence

Public runs [40315351988](https://github.com/mashingaan/open-exoplanet-hunt/actions/runs/40315351988)
and [40315612271](https://github.com/mashingaan/open-exoplanet-hunt/actions/runs/40315612271)
both passed the repository verifier, 12 deterministic tests and all six fixed
controls at `99ae48ad78a2ddf3d4566664e393214fd4d907e4`. The second run's
`--compare` passed. Complete science receipts and all six spectral containers
are byte-identical. Both public hosts reported AMD EPYC 9V45. The separate
development pair tested two different EPYC models. Neither establishes
universal architecture independence.

The fork tag `tls-controls-v1-experiment` retains the executed source. A repeat
can be dispatched in `mashingaan/open-exoplanet-hunt` with that tag as the ref
and `comparison_run_id=40315351988` while the original artifact remains
available. After artifact expiry, regenerate a first run from the tag and
compare a second run to its ID. Compact receipts remain in git, with spectral
hashes and complete inputs needed to regenerate the numerical artifacts.

This evaluation adds an independent transit period-search component before
using it on astronomical targets. It does not create candidates or modify the
existing ten-gate thresholds. The fixed experiment is
[tls-v1.json](../controls/tls-v1.json), retained in git before execution.
If a declared control fails, preserve its outputs and investigate rather than
editing seeds, thresholds, periods or aliases to manufacture a pass.

## Source and method decisions

[TLS documentation](https://transitleastsquares.readthedocs.io/en/latest/Python%20interface.html)
describes normalized flux and day units, multiple backends, and sampling caveats.
The source is pinned to `7666720794a3cc064535055917500c96dbd90b52`, version
2.0.1. Its `fused` backend retains the exact statistic rather than the default
piecewise-linear approximation. A positive-control spectrum is compared against
the `numba` reference backend at the preregistered tolerances. This one comparison
does not establish equivalence for every dataset or platform.

The underlying method citation is Hippke and Heller 2019,
[DOI 10.1051/0004-6361/201834672](https://doi.org/10.1051/0004-6361/201834672).
The paper predates TLS 2.0, so it is not evidence that its newer backends have
received independent peer review. The pinned source and executed controls are
recorded separately from the publication.

The independent injection generator is
[batman 2.5.3](https://github.com/lkreidberg/batman/tree/6fc570e209fa2e3dadd2ffa4c535f528a7bef495).
Its source documentation requests Kreidberg 2015,
[DOI 10.1086/683602](https://doi.org/10.1086/683602), bibcode
`2015PASP..127.1161K`. It uses quadratic limb darkening and finite-exposure
integration with seven subsamples. The TLS search uses its own model, not
batman at runtime. The secondary-eclipse control superposes scaled transit
shapes and is not a complete binary-star model.

## Fixed controls and interpretation

Six controls share three fourteen-day observing windows at thirty-minute
cadence and matched fixed-seed PCG64 noise. The positive transit, white noise,
rotation, secondary eclipse, missing third-sector signal and centroid blend
cases are all declared before results exist. Astropy BLS and TLS search the
same normalized flux independently over 1.1 to 8 days. The positive case must
recover its injected period in both detectors. Each TLS proposal receives the
unchanged upstream ten-gate evaluation and all four predeclared aliases.
Specified negative controls must fail their declared relevant veto, as well as
leaving no ten-gate survivor at any tested alias.

The downstream evaluator performs a narrow BLS refinement around a supplied TLS
period. Thus this is independent detection followed by shared falsification,
not a fully independent TLS likelihood or parameter-estimation pipeline. Both
the original TLS result and the refined result are retained. TLS SDE, BLS SNR,
periodogram half-width and TLS white-noise FAP remain distinct quantities.
None is treated as an astrophysical false-positive probability or confirmation.

The adapter refuses invalid samples instead of permitting silent TLS cleaning,
requires positive sorted times and normalized flux, checks the executed source
commit and bounds threads. Its current readiness guard is 50000 samples and a
120-day baseline. Longer, sparse and mixed-cadence datasets still need dedicated
validation, as do actual MAST recovery and independent photometric reductions.
Passing these six controls would not provide a completeness or false-alarm
calibration from one noise realization.

The first complete pair of executions passed every scientific control but
failed exact replay in TLS spectral hashes and last-digit SDE values.
The [retained report](../validation/tls-controls-v1/REPORT.md) includes both
results and all seventeen differing fields. The subsequent environment pair
passed all six controls and exact replay at the same immutable specification.
Its complete science receipts and all six spectral containers are byte-identical
on two different AMD EPYC models. Real-target deployment remains unauthorized.
This environment keeps the exact comparison and selects `NUMBA_CPU_NAME=generic` and
`NPY_ENABLE_CPU_FEATURES=SSE3`, giving the pinned x86 NumPy wheel its SSE/SSE2/SSE3
baseline with no additional dispatched targets. The pinned
[NumPy environment parser](https://github.com/numpy/numpy/blob/v2.3.4/numpy/_core/src/common/npy_cpu_features.c)
preserves baseline features and disables unselected dispatch targets.
[Numba documents](https://numba.readthedocs.io/en/stable/reference/envvars.html)
the generic CPU target. These results establish reproducibility for this bounded
pair, not the mechanism of the original discrepancy. Runtime and CPU receipts
accompany subsequent runs. The original failed repeat remains in the record.

## Bounded hosted execution

The manual `Independent TLS detector controls` workflow has a fifteen-minute
limit and one numerical thread. Dependencies are pinned in
[tls-requirements.txt](../tls-requirements.txt). TLS is installed from the exact
Git commit and its direct-url provenance is checked before execution. Binary
distributions are required for batman, Numba and llvmlite. Pip and immutable JIT
caches are hosted. The ordinary discovery workflows do not install TLS.

The workflow runs the repository verifier and deterministic tests before the
control suite. It retains per-control detector and falsification receipts,
warnings, numerical periodograms, provenance and diagnostic plots, including
failed controls. These numerical artifacts are kept in Actions, with no FITS
or bulk data committed. The astronomical scan registry must remain byte-identical.

For an independent repeat, dispatch the same frozen branch with
`comparison_run_id=FIRST_RUN_ID`. It generates the original inputs again and
requires equality with the complete prior `science.json`. Any disagreement
remains a recorded failure pending investigation. Component registration alone
does not establish readiness for real-target analysis.

## Bootstrap compatibility failure retained

Run [40312269010](https://github.com/mashingaan/tess-discovery/actions/runs/40312269010)
at `bffa0b373523e0627c4aa81562f43b861ce3521a` installed the pinned numerical
dependencies and passed the repository verifier and all 36 tests. The control
process failed before generating data because batman imports `distutils`,
which Python 3.12 no longer provides. This was a runtime packaging failure,
with no detector results or changes to the preregistered experiment.

The environment now explicitly pins setuptools 80.9.0. Its maintained vendored
distutils is described in the
[setuptools porting guide](https://setuptools.pypa.io/en/latest/deprecated/distutils-legacy.html).
The control specification, generators, seeds, thresholds and expected vetoes
remain unchanged. The failed run is retained rather than relabeled a pass.

The two method citations were checked against original package documentation
and Crossref DOI records. A compact
[metadata receipt](../reference-data/tls-method-citations.json) retains their
verified titles, authors, dates and response hashes. This verifies citation
metadata, not the current software's numerical behavior.
