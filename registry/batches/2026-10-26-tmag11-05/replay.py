"""Replay recorded private evaluations without selecting or promoting targets."""
import tarfile
import tempfile
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
import vet_candidate as vet


def main():
    batch = Path(sys.argv[1]).resolve()
    manifest = json.loads((batch / "manifest.json").read_text())
    for name, expected in manifest.items():
        path = (batch / name).resolve()
        if not path.is_relative_to(batch) or sha256(path.read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"Receipt integrity failure: {name}")
    provenance = json.loads((batch / "provenance.json").read_text())
    for name, expected in provenance["scientific_scripts"].items():
        if sha256((ROOT / "scripts" / name).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"Scientific script changed: {name}")
    snapshot = batch / "official-known-tics.json.gz"
    if snapshot.exists():
        vet.KNOWN_TICS_GZ = snapshot
    else:
        # The initial private run predates retention of the full official index.
        # Its gate decisions remain recorded. Require exact equality below and
        # disclose that catalog replay uses the inherited repository snapshot.
        print("Catalog replay limitation: using inherited repository index, not the live search-time index")
    receipts = sorted(batch.glob("TIC-*-evaluation-*.json"))
    if not receipts:
        raise RuntimeError("No completed scientific evaluations to reproduce")
    for path in receipts:
        receipt = json.loads(path.read_text())
        tic = receipt["tic_id"]
        sectors = []
        for product in receipt["sectors_used"]:
            sector = vet.download_and_load_sector(tic, product)
            if sector["sha256"] != product["sha256"]:
                raise RuntimeError(f"Public FITS changed: TIC {tic}, sector {product['sector']}")
            sectors.append(sector)
        result = vet.evaluate_signal_from_arrays(receipt["stellar_parameters"], sectors, **receipt["search_arguments"])
        clean = {key: value for key, value in result.items() if key != "_arrays"}
        if clean != receipt["result"]:
            raise RuntimeError(f"Scientific result differs: {path.name}")
        aliases_path = batch / f"TIC-{tic}-aliases.json"
        # Supplemental aliases use the final sector set, not a superseded screen.
        later = [item for item in receipts if item.name.startswith(f"TIC-{tic}-") and item.name > path.name]
        if aliases_path.exists() and not later:
            for alias in json.loads(aliases_path.read_text()):
                period = receipt["result"]["ephemeris"]["period_d"] * alias["multiplier"]
                result = vet.evaluate_signal_from_arrays(
                    receipt["stellar_parameters"], sectors, target_period=period,
                    period_min=max(0.4, period * 0.98), period_max=period * 1.02, allow_known=False,
                )
                clean = {key: value for key, value in result.items() if key != "_arrays"}
                if clean != alias["result"]:
                    raise RuntimeError(f"Alias result differs: TIC {tic}, multiplier {alias['multiplier']}")
        print(f"Exact scientific evaluation reproduced: {path.name}")


if __name__ == "__main__":
    archive_path = Path(__file__).with_name("receipts.tar.gz")
    expected = "01ac07a210269a39362c9802cf24b44cafc460fba735a7e16e3612d2e958653d"
    if sha256(archive_path.read_bytes()).hexdigest() != expected:
        raise RuntimeError("Published receipt archive hash mismatch")
    with tempfile.TemporaryDirectory() as directory:
        with tarfile.open(archive_path) as archive:
            archive.extractall(directory, filter="data")
        sys.argv = [sys.argv[0], directory]
        main()
