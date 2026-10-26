"""Fork-only execution harness retaining receipts for the unchanged hunter."""
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import hunt_sectors as hunt
import vet_candidate as vet
import sync_catalogs as catalogs

OUT = ROOT / "batch-output"
OUT.mkdir(exist_ok=True)


def save(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


config = {
    "max_targets": int(os.environ["BATCH_MAX_TARGETS"]),
    "min_sectors": 3,
    "tmag_min": float(os.environ["BATCH_TMAG_MIN"]),
    "tmag_max": float(os.environ["BATCH_TMAG_MAX"]),
    "teff_min": 3300.0,
    "teff_max": 4300.0,
    "dist_max": 55.0,
    "page": int(os.environ["BATCH_PAGE"]),
}
assert 1 <= config["max_targets"] <= 2
baseline_ids = hunt.load_known_and_scanned()[1]
save("provenance.json", {
    "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
    "run_url": f"https://github.com/{os.environ['GITHUB_REPOSITORY']}/actions/runs/{os.environ['GITHUB_RUN_ID']}",
    "python": platform.python_version(),
    "dependencies": subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True).splitlines(),
    "script_sha256": {p.name: digest(p) for p in (ROOT / "scripts").glob("*.py")},
    "baseline_registry_sha256": digest(hunt.SCANNED_LOG),
    "baseline_known_tics_sha256": digest(vet.KNOWN_TICS_GZ),
    "config": config,
    "period_min_d": 1.1,
    "period_max_d": 15.0,
    "selection_policy": "Upstream run_hunt order and fixed six-sector initial screen, expand only as upstream specifies",
})

# Refresh official exclusions before any target is searched, preserving the
# repository catalog and recording the downloaded source bytes for audit.
fresh = ROOT / ".cache/fresh-catalog"
fresh.mkdir(parents=True, exist_ok=True)
catalogs.CATALOG_DIR = fresh
catalogs.SUMMARY_PATH = fresh / "catalog_summary.json"
catalogs.KNOWN_TICS_GZ = fresh / "known_tics.json.gz"
fetch_original = catalogs.fetch_bytes
source_receipts = []


def fetch_receipt(url, timeout=60):
    raw = fetch_original(url, timeout)
    source_receipts.append({"url": url, "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)})
    save("catalog-source-receipts.json", source_receipts)
    return raw


catalogs.fetch_bytes = fetch_receipt
summary = catalogs.build_catalog()
save("fresh-catalog-summary.json", summary)
vet.KNOWN_TICS_GZ = catalogs.KNOWN_TICS_GZ
hunt.KNOWN_TICS_GZ = catalogs.KNOWN_TICS_GZ

invoke_original = vet.mast_invoke
query_count = 0


def invoke_receipt(request_obj, timeout=60, retries=3):
    global query_count
    result = invoke_original(request_obj, timeout, retries)
    query_count += 1
    save(f"mast-query-{query_count:04d}.json", {"request": request_obj, "response": result})
    return result


vet.mast_invoke = invoke_receipt
hunt.mast_invoke = invoke_receipt
evaluate_original = hunt.evaluate_signal_from_arrays
evaluation_count = 0


def evaluate_receipt(meta, sectors, **kwargs):
    global evaluation_count
    evaluation_count += 1
    result = evaluate_original(meta, sectors, **kwargs)
    receipts = [{k: s[k] for k in (
        "sector", "pipeline", "exptime_s", "uri", "flux_col", "sha256", "n_points", "mad_ppm"
    )} for s in sectors]
    coverage = [{
        "sector": s["sector"],
        "n_finite_centroids": int((vet.np.isfinite(s["dx"]) & vet.np.isfinite(s["dy"])).sum()),
        "start_btjd": float(s["t"].min()),
        "end_btjd": float(s["t"].max()),
    } for s in sectors]
    save(f"TIC-{meta['tic_id']}-evaluation-{evaluation_count}.json", {
        "stellar_parameters": meta,
        "sectors_used": receipts,
        "coverage": coverage,
        "search_arguments": kwargs,
        "result": {k: v for k, v in result.items() if k != "_arrays"},
    })
    return result


hunt.evaluate_signal_from_arrays = evaluate_receipt
results = hunt.run_hunt(**config)
assert len(results) == config["max_targets"], "Incomplete batch, inspect query receipts and logs"
assert all(r["tic_id"] not in baseline_ids for r in results)
save("new-registry-rows.json", results)
save("exclusion-check.json", {
    str(r["tic_id"]): vet.check_known_catalog(r["tic_id"]) for r in results
})
print(json.dumps(results, indent=2))
