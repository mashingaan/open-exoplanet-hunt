"""Replay this fixed negative batch from public MAST FITS with hash checks."""
from hashlib import sha256
import json
from pathlib import Path
import sys

BATCH = Path(__file__).resolve().parent
ROOT = BATCH.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import vet_candidate as vet

provenance = json.loads((BATCH / "search/provenance.json").read_text())
for name, expected in provenance["script_sha256"].items():
    actual = sha256((ROOT / "scripts" / name).read_bytes()).hexdigest()
    if actual != expected:
        raise RuntimeError(f"Scientific script hash mismatch: {name}. Use the recorded upstream commit and LF line endings.")

for path in sorted((BATCH / "search").glob("TIC-*-evaluation-*.json")):
    receipt = json.loads(path.read_text())
    tic = receipt["stellar_parameters"]["tic_id"]
    sectors = []
    for product in receipt["sectors_used"]:
        data = vet.download_and_load_sector(tic, product)
        if data["sha256"] != product["sha256"]:
            raise RuntimeError(f"FITS hash mismatch: TIC {tic}, sector {product['sector']}")
        sectors.append(data)
    result = vet.evaluate_signal_from_arrays(
        receipt["stellar_parameters"], sectors, **receipt["search_arguments"]
    )
    actual = {key: value for key, value in result.items() if key != "_arrays"}
    if actual != receipt["result"]:
        raise RuntimeError(f"Evaluation differs for TIC {tic}. Compare dependencies and the recorded catalog snapshot.")
    print(f"PASS TIC {tic}: exact recorded primary evaluation reproduced")
