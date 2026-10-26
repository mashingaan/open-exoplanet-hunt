"""Execute preregistered synthetic TLS/BLS controls without scanning real targets."""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import warnings

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import batman
import numba
import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as const
from astropy.timeseries import BoxLeastSquares

from discovery.adapters.tls import array_hash, installed_provenance, search
import vet_candidate as vet

OUT = ROOT / "tls-output"
SPEC = ROOT / "discovery/controls/tls-v1.json"


def write(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def generate(spec, case):
    config = spec["generator"]
    rng = np.random.Generator(np.random.PCG64(config["flux_seed"]))
    centroids = np.random.Generator(np.random.PCG64(config["centroid_seed"]))
    params = batman.TransitParams()
    params.t0 = config["epoch_days"]
    params.per = config["period_days"]
    params.rp = config["radius_ratio"]
    seconds = params.per * 86400
    semimajor_m = (const.G.value * const.M_sun.value * config["stellar_mass_solar"] * (seconds / (2 * np.pi)) ** 2) ** (1 / 3)
    params.a = semimajor_m / (const.R_sun.value * config["stellar_radius_solar"])
    params.inc = config["inclination_deg"]
    params.ecc = config["eccentricity"]
    params.w = config["omega_deg"]
    params.u = config["limb_darkening"]
    params.limb_dark = config["limb_darkening_law"]
    cadence = config["exposure_seconds"] / 86400
    sectors = []
    for index, (start, end) in enumerate(config["windows_days"]):
        t = np.arange(start, end, cadence)
        model = batman.TransitModel(params, t, supersample_factor=config["supersample_factor"],
                                    exp_time=cadence, nthreads=1)
        transit = model.light_curve(params) - 1.0
        signal = case["transit_scale"] * transit
        if index in case.get("remove_signal_sector_indices", []):
            signal = np.zeros_like(t)
        if case.get("secondary_scale"):
            secondary = deepcopy(params)
            secondary.t0 += secondary.per / 2
            second_model = batman.TransitModel(secondary, t, supersample_factor=config["supersample_factor"],
                                               exp_time=cadence, nthreads=1)
            signal += case["secondary_scale"] * (second_model.light_curve(secondary) - 1.0)
        if case.get("sine_amplitude"):
            signal += case["sine_amplitude"] * np.sin(2 * np.pi * (t - params.t0) / params.per)
        flux = 1 + signal + rng.normal(0, config["flux_sigma"], len(t))
        dx = centroids.normal(0, config["centroid_sigma_pixels"], len(t))
        dy = centroids.normal(0, config["centroid_sigma_pixels"], len(t))
        dx[transit < -1e-8] += case.get("centroid_shift_pixels", 0)
        sectors.append({
            "sector": index + 1, "pipeline": "SYNTHETIC_BATMAN_CONTROL", "t": t,
            "f": flux, "f_sap": flux.copy(), "e": np.full(len(t), config["flux_sigma"]),
            "dx": dx, "dy": dy, "exptime_s": config["exposure_seconds"],
            "uri": f"synthetic://tls-v1/{case['id']}/{index}", "flux_col": "SYNTHETIC",
            "n_points": len(t), "mad_ppm": config["flux_sigma"] * 1e6,
            "sha256": array_hash(t, flux, dx, dy),
        })
    # This inherited synthetic fixture identity is never written to a scan registry.
    meta = {"tic_id": 999999999, "rad_rsun": config["stellar_radius_solar"],
            "mass_msun": config["stellar_mass_solar"], "ra_deg": 80.0, "dec_deg": -1.0,
            "tmag": 9.7, "teff_k": 4000.0, "dist_pc": 35.0, "contratio": 0.01}
    return meta, sectors


def clean_vetting(meta, sectors, period):
    result = vet.evaluate_signal_from_arrays(meta, sectors, target_period=period,
                                            period_min=1.1, period_max=8.0, allow_known=False)
    return {key: value for key, value in result.items() if key != "_arrays"}


def run(spec):
    figures, axes = plt.subplots(len(spec["cases"]), 3, figsize=(15, 3 * len(spec["cases"])))
    records = []
    for index, case in enumerate(spec["cases"]):
        print(f"Control {case['id']}, fixed inputs and parameters", flush=True)
        meta, sectors = generate(spec, case)
        t, flux, error = [np.concatenate([sector[key] for sector in sectors]) for key in ("t", "f", "e")]
        tls, raw = search(t, flux, error, spec["tls"]["power"])
        cfg = spec["bls"]
        periods = 1 / np.linspace(1 / cfg["period_max"], 1 / cfg["period_min"], cfg["frequency_grid_points"])[::-1]
        durations = np.linspace(cfg["duration_min_days"], cfg["duration_max_days"], cfg["duration_grid_points"])
        bls_raw = BoxLeastSquares(t, flux, dy=error).power(periods, durations, objective=cfg["objective"])
        best = int(np.argmax(bls_raw.power))
        bls = {"period": float(bls_raw.period[best]), "snr": float(bls_raw.power[best]),
               "epoch": float(bls_raw.transit_time[best]), "duration": float(bls_raw.duration[best]),
               "spectrum_sha256": array_hash(bls_raw.period, bls_raw.power)}
        write(f"detectors-{case['id']}.json", {"tls": tls, "bls": bls})
        np.savez_compressed(OUT / f"control-{case['id']}-spectra.npz",
                            tls_period=raw.periods, tls_power=raw.power, tls_chi2=raw.chi2,
                            tls_power_raw=raw.power_raw, bls_period=bls_raw.period, bls_power=bls_raw.power)
        primary = clean_vetting(meta, sectors, tls["scalars"]["period"])
        aliases = [{"multiplier": multiplier,
                    "result": clean_vetting(meta, sectors, tls["scalars"]["period"] * multiplier)}
                   for multiplier in spec["fixed_alias_multipliers"]]
        checks = {}
        if case["id"] == "positive":
            truth = spec["generator"]["period_days"]
            tolerance = spec["positive_recovery_relative_period_tolerance"]
            checks["tls_period_recovery"] = abs(tls["scalars"]["period"] / truth - 1) <= tolerance
            checks["bls_period_recovery"] = abs(bls["period"] / truth - 1) <= tolerance
            checks["primary_passes_all_ten"] = primary["all_gates_passed"]
        else:
            checks["no_primary_or_alias_survives"] = not any(
                result["all_gates_passed"] for result in [primary] + [alias["result"] for alias in aliases])
            if case.get("primary_veto_any"):
                checks["specified_primary_veto"] = any(not primary["gates"][name]["pass"] for name in case["primary_veto_any"])
        reference = None
        if case["id"] == spec["tls"]["reference_control"]:
            reference_cfg = {**spec["tls"]["power"], "backend": spec["tls"]["reference_backend"]}
            reference, reference_raw = search(t, flux, error, reference_cfg)
            checks["reference_period_grid_identical"] = bool(np.array_equal(raw.periods, reference_raw.periods))
            checks["reference_spectrum_matches"] = bool(np.allclose(
                raw.power, reference_raw.power,
                atol=spec["tls"]["backend_spectrum_absolute_tolerance"],
                rtol=spec["tls"]["backend_spectrum_relative_tolerance"], equal_nan=False))
        record = {"case": case, "input_sha256": array_hash(t, flux, error),
                  "sector_input_hashes": [sector["sha256"] for sector in sectors],
                  "tls": tls, "bls": bls, "reference_tls": reference,
                  "primary_upstream_vetting": primary, "fixed_alias_vetting": aliases,
                  "checks": checks, "control_passed": all(checks.values())}
        records.append(record)
        # Commit each result before proceeding so a failed or interrupted suite
        # retains every completed diagnostic rather than discarding failures.
        write(f"control-{case['id']}.json", record)
        phase = ((t - raw.T0 + raw.period / 2) % raw.period) - raw.period / 2
        axes[index, 0].scatter(phase, (flux - 1) * 1e6, s=2, alpha=0.35)
        axes[index, 0].set_title(case["id"] + ": TLS folded data")
        axes[index, 1].plot(raw.periods, raw.power)
        axes[index, 1].set_title(f"TLS P={raw.period:.6f} d, SDE={raw.SDE:.2f}")
        axes[index, 2].plot(bls_raw.period, bls_raw.power)
        axes[index, 2].set_title(f"BLS P={bls['period']:.6f} d, SNR={bls['snr']:.2f}")
        print(case["id"], checks, flush=True)
    for row in axes:
        row[0].set_xlabel("Days from TLS mid-transit")
        row[0].set_ylabel("Relative flux, ppm")
        for axis in row[1:]:
            axis.set_xlabel("Trial period, days")
    figures.suptitle("Synthetic controls only, not astronomical candidates", y=1.002)
    figures.tight_layout()
    figures.savefig(OUT / "diagnostics.png", dpi=130, bbox_inches="tight")
    plt.close(figures)
    return {"experiment_id": spec["experiment_id"], "spec_sha256": sha256(SPEC.read_bytes()).hexdigest(),
            "all_controls_passed": all(record["control_passed"] for record in records), "controls": records}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--compare", type=Path, help="Require exact equality with a prior complete science receipt")
    args = parser.parse_args()
    if OUT.exists() and any(OUT.iterdir()):
        raise ValueError("Controls require a fresh output directory")
    OUT.mkdir(exist_ok=True)
    spec = json.loads(SPEC.read_text())
    if version("batman-package") != spec["generator"]["version"]:
        raise ValueError("Independent injection generator version differs from the preregistration")
    registry_before = (ROOT / "registry/scanned_targets.jsonl").read_bytes()
    provenance = {"mode": "SYNTHETIC_VALIDATION_NO_ASTRONOMICAL_SEARCH",
                  "started_utc": datetime.now(timezone.utc).isoformat(), "run_id": os.environ.get("GITHUB_RUN_ID"),
                  "pipeline_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
                  "command_argv": ["python", "-B", *sys.argv], "python": sys.version, "platform": platform.platform(),
                  "spec": spec, "spec_sha256": sha256(SPEC.read_bytes()).hexdigest(),
                  "dependencies": subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True).splitlines(),
                  "tls_source": installed_provenance(), "batman_distribution_version": version("batman-package"),
                  "batman_module_version": getattr(batman, "__version__", None),
                  "scientific_scripts": {path.name: sha256(path.read_bytes()).hexdigest() for path in (ROOT / "scripts").glob("*.py")}}
    write("provenance.json", provenance)
    # Keep host-specific diagnostics outside the exact scientific comparison.
    # NumPy's private feature table is inspected only for the pinned version.
    from numpy._core import _multiarray_umath
    cpuinfo = Path("/proc/cpuinfo")
    cpu_lines = cpuinfo.read_text().splitlines() if cpuinfo.exists() else []
    write("runtime.json", {
        "environment": {key: os.environ.get(key) for key in (
            "NUMBA_CPU_NAME", "NUMBA_CPU_FEATURES", "NPY_ENABLE_CPU_FEATURES",
            "NPY_DISABLE_CPU_FEATURES", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
            "NUMBA_NUM_THREADS", "OPENBLAS_CORETYPE")},
        "numba_cpu_name": numba.config.CPU_NAME,
        "numba_cpu_features": numba.config.CPU_FEATURES,
        "numpy_cpu_baseline": list(_multiarray_umath.__cpu_baseline__),
        "numpy_cpu_dispatch": list(_multiarray_umath.__cpu_dispatch__),
        "numpy_cpu_features": dict(_multiarray_umath.__cpu_features__),
        "numpy_configuration": np.show_config(mode="dicts"),
        "host_cpu_description": sorted({line.strip() for line in cpu_lines
                                        if line.startswith(("vendor_id", "model name", "flags"))}),
    })
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            science = run(spec)
            write("science.json", science)
            if args.compare and science != json.loads(args.compare.read_text()):
                raise ValueError("Independent control replay differs, preserve both receipts for investigation")
            if not science["all_controls_passed"]:
                raise ValueError("A preregistered control failed, deployment remains blocked")
        except Exception as error:
            write("failure.json", {"type": type(error).__name__, "message": str(error),
                                   "deployment_authorized": False})
            raise
        finally:
            write("warnings.json", [{"category": item.category.__name__, "message": str(item.message)} for item in caught])
            if (ROOT / "registry/scanned_targets.jsonl").read_bytes() != registry_before:
                raise ValueError("Synthetic controls must never mutate astronomical scan coverage")
            files = {path.name: sha256(path.read_bytes()).hexdigest() for path in OUT.iterdir() if path.is_file()}
            write("manifest.json", files)


if __name__ == "__main__":
    main()
