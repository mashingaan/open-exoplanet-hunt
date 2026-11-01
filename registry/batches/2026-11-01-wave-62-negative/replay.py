"""Reproduce one frozen negative TESS wave partition from public MAST data."""
import argparse
from hashlib import sha256
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile

def load_engine(batch, provenance):
    source_root=batch/"scientific-source"/"scripts"
    bundled=provenance["scientific_engine"]["bundled_script_sha256"]
    for name,expected in bundled.items():
        if sha256((source_root/name).read_bytes()).hexdigest()!=expected:
            raise RuntimeError(f"Pinned scientific engine source hash mismatch: {name}")
    root=Path(tempfile.gettempdir())/"tess-wave-62-engine"
    scripts=root/"scripts"
    scripts.mkdir(parents=True,exist_ok=True)
    path=scripts/"vet_candidate.py"
    path.write_bytes((source_root/"vet_candidate.py").read_bytes())
    sys.path.insert(0,str(scripts))
    spec=importlib.util.spec_from_file_location("vet_candidate",path)
    if spec is None or spec.loader is None: raise RuntimeError("Could not load the pinned scientific engine")
    vet=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=vet
    spec.loader.exec_module(vet)
    vet.KNOWN_TICS_GZ=batch/"target-catalog.json.gz"
    cache=os.environ.get("TESS_REPLAY_CACHE")
    vet.CACHE_DIR=Path(cache) if cache else root/".cache"/"mast"
    vet.CACHE_DIR.mkdir(parents=True,exist_ok=True)
    return vet

def replay_target(vet,receipts,batch):
    for path in receipts:
        receipt=json.loads(path.read_text())
        tic=receipt["tic_id"]
        sectors=[]
        for product in receipt["sectors_used"]:
            sector=vet.download_and_load_sector(tic,product)
            if sector["sha256"]!=product["sha256"]: raise RuntimeError(f"Public FITS changed for TIC {tic}, sector {product['sector']}")
            sectors.append(sector)
        result=vet.evaluate_signal_from_arrays(receipt["stellar_parameters"],sectors,**receipt["search_arguments"])
        clean={key:value for key,value in result.items() if key!="_arrays"}
        if clean!=receipt["result"]: raise RuntimeError(f"Scientific result differs: {path.name}")
        aliases_path=batch/f"TIC-{tic}-aliases.json"
        if aliases_path.exists() and path==receipts[-1]:
            for alias in json.loads(aliases_path.read_text()):
                period=receipt["result"]["ephemeris"]["period_d"]*alias["multiplier"]
                result=vet.evaluate_signal_from_arrays(receipt["stellar_parameters"],sectors,target_period=period,period_min=max(0.4,period*0.98),period_max=period*1.02,allow_known=False)
                clean={key:value for key,value in result.items() if key!="_arrays"}
                if clean!=alias["result"]: raise RuntimeError(f"Alias result differs for TIC {tic}, multiplier {alias['multiplier']}")
    return [path.name for path in receipts]

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--partition",type=int,choices=range(8),required=True)
    args=parser.parse_args()
    batch=Path(__file__).resolve().parent
    manifest=json.loads((batch/"manifest.json").read_text())
    for name,expected in manifest.items():
        path=(batch/name).resolve()
        if not path.is_relative_to(batch) or sha256(path.read_bytes()).hexdigest()!=expected: raise RuntimeError(f"Receipt integrity failure: {name}")
    actual={p.relative_to(batch).as_posix() for p in batch.rglob("*") if p.is_file() and p.name!="manifest.json"}
    if actual!=set(manifest): raise RuntimeError("Batch has missing or unmanifested files")
    provenance=json.loads((batch/"provenance.json").read_text())
    if sha256((batch/"target-catalog.json.gz").read_bytes()).hexdigest()!=provenance["scientific_engine"]["target_catalog_sha256"]: raise RuntimeError("Target catalog snapshot hash mismatch")
    vet=load_engine(batch,provenance)
    allowed=set(provenance["completed_targets_by_partition"][str(args.partition)])
    receipts=sorted(p for p in batch.glob("TIC-*-evaluation-*.json") if int(json.loads(p.read_text())["tic_id"]) in allowed)
    groups={}
    for path in receipts:
        tic=int(json.loads(path.read_text())["tic_id"])
        groups.setdefault(tic,[]).append(path)
    if not groups: raise RuntimeError("No completed scientific evaluations selected for replay")
    for paths in groups.values():
        for name in replay_target(vet,paths,batch): print(f"Exact public-MAST replay passed: {name}",flush=True)

if __name__=="__main__": main()
