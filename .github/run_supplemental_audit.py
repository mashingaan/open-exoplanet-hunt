"""Fixed rational aliases and neighbor checks on exactly the recorded FITS."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import vet_candidate as vet
import audit_neighbors_and_aliases as audit

OUT = ROOT / "audit-output"
OUT.mkdir(exist_ok=True)
SOURCE = ROOT / "prior-batch/batch-output"


def save(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


provenance = json.loads((SOURCE / "provenance.json").read_text())
for name, expected in provenance["script_sha256"].items():
    assert hashlib.sha256((ROOT / "scripts" / name).read_bytes()).hexdigest() == expected
save("provenance.json", {
    "source": provenance,
    "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
    "run_url": f"https://github.com/{os.environ['GITHUB_REPOSITORY']}/actions/runs/{os.environ['GITHUB_RUN_ID']}",
    "fixed_alias_multipliers": [2.0, 0.5, 1.5, 2.0 / 3.0],
    "policy": "Audit original searched sectors and verify every FITS hash before evaluation. Alias results do not replace the search result.",
})
query_original = audit.mast_invoke
queries = []


def query_receipt(request_obj, **kwargs):
    response = query_original(request_obj, **kwargs)
    queries.append({"request": request_obj, "response": response})
    save("neighbor-query-receipts.json", queries)
    return response


audit.mast_invoke = query_receipt
rows = json.loads((SOURCE / "new-registry-rows.json").read_text())
for row in rows:
    tic = row["tic_id"]
    matches = sorted(SOURCE.glob(f"TIC-{tic}-evaluation-*.json"))
    receipt = json.loads(matches[-1].read_text())
    sectors = []
    for prod in receipt["sectors_used"]:
        loaded = vet.download_and_load_sector(tic, prod)
        assert loaded["sha256"] == prod["sha256"], f"FITS hash changed for TIC {tic} sector {prod['sector']}"
        sectors.append(loaded)
    cert = {
        "tic_id": tic,
        "stellar_parameters": receipt["stellar_parameters"],
        "primary_signal": receipt["result"],
    }
    neighbors = audit.audit_target_neighbors(cert)
    save(f"TIC-{tic}-neighbors.json", neighbors)
    aliases = []
    for mult in (2.0, 0.5, 1.5, 2.0 / 3.0):
        p = row["best_period_d"] * mult
        result = vet.evaluate_signal_from_arrays(
            receipt["stellar_parameters"], sectors,
            target_period=p, period_min=max(0.4, p * 0.98),
            period_max=p * 1.02, allow_known=False,
        )
        aliases.append({
            "multiplier": mult,
            "result": {k: v for k, v in result.items() if k != "_arrays"},
        })
        save(f"TIC-{tic}-aliases.json", aliases)
        print(tic, mult, result["ephemeris"], result["all_gates_passed"], flush=True)
