"""Bounded independent transit search, with no automatic scientific promotion."""
from hashlib import sha256
from importlib.metadata import distribution, version
import json

import numpy as np


TLS_COMMIT = "7666720794a3cc064535055917500c96dbd90b52"


def array_hash(*arrays):
    digest = sha256()
    for array in arrays:
        data = np.asarray(array, dtype="<f8")
        digest.update(json.dumps(data.shape).encode() + b"\n")
        digest.update(data.tobytes(order="C"))
    return digest.hexdigest()


def installed_provenance():
    dist = distribution("transitleastsquares")
    direct = json.loads(dist.read_text("direct_url.json") or "{}")
    if version("transitleastsquares") != "2.0.1" or direct.get("vcs_info", {}).get("commit_id") != TLS_COMMIT:
        raise ValueError("TLS must be installed from the exact reviewed source commit")
    return {"version": version("transitleastsquares"), "direct_url": direct,
            "role": "independent detector, never planetary confirmation"}


def search(t, flux, error, parameters):
    arrays = [np.asarray(values, dtype=float) for values in (t, flux, error)]
    if any(values.ndim != 1 or not np.all(np.isfinite(values)) for values in arrays):
        raise ValueError("TLS input arrays must be finite and one-dimensional")
    if len({len(values) for values in arrays}) != 1 or len(t) < 50:
        raise ValueError("TLS input arrays must be aligned with at least fifty cadences")
    if np.any(np.diff(arrays[0]) <= 0) or np.any(arrays[0] <= 0):
        raise ValueError("TLS times must be positive and strictly increasing, with the origin recorded separately")
    if np.any(arrays[1] <= 0) or np.any(arrays[2] <= 0):
        raise ValueError("TLS flux and uncertainty must be positive")
    if len(t) > 50000 or arrays[0][-1] - arrays[0][0] > 120:
        raise ValueError("This readiness adapter is not validated for more than 50000 samples or 120 days")
    if not 0.99 <= np.median(arrays[1]) <= 1.01:
        raise ValueError("TLS flux must already be normalized, never silently rescaled")
    if type(parameters["use_threads"]) is not int or parameters["use_threads"] != 1 or not 0.4 <= parameters["period_min"] < parameters["period_max"] <= 15:
        raise ValueError("Search exceeds the bounded adapter policy")
    installed_provenance()
    from transitleastsquares import transitleastsquares

    model = transitleastsquares(*arrays, verbose=False)
    # TLS can clean inputs internally. Never permit silent sample removal.
    if not np.array_equal(model.t, arrays[0]) or not np.array_equal(model.y, arrays[1]):
        raise ValueError("TLS silently changed the input sample population")
    result = model.power(**parameters)
    if np.min(result.periods) < parameters["period_min"] - 1e-10 or np.max(result.periods) > parameters["period_max"] + 1e-10:
        raise ValueError("TLS period grid expanded beyond the preregistered search bounds")
    scalars = {}
    for key in ("period", "period_uncertainty", "T0", "duration", "depth", "SDE", "SDE_raw", "snr", "FAP", "odd_even_mismatch", "distinct_transit_count"):
        value = float(result[key])
        scalars[key] = value if np.isfinite(value) else None
    return {"scalars": scalars, "input_sha256": array_hash(*arrays),
            "spectrum_sha256": array_hash(result.periods, result.power),
            "period_grid_points": len(result.periods), "parameters": parameters,
            "status": "DETECTOR_OUTPUT_REQUIRES_FALSIFICATION",
            "parameter_notes": {"period_uncertainty": "TLS periodogram half-width, not a posterior",
                                "uncertainties": "TLS normalizes dy internally as relative weights",
                                "FAP": "TLS white-noise lookup, not astrophysical false-positive probability",
                                "depth": "TLS returned flux level at the transit bottom"}}, result
