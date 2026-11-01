#!/usr/bin/env python3
"""Download and index confirmed exoplanets, NASA TESS TOIs, ExoFOP CTOIs,
and curated external AI/ML/citizen-science discoveries into a unified fast lookup index.
"""
import sys
if sys.flags.optimize:
    raise ValueError("Assertions must remain enabled (do not run with -O)")

import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parent.parent
CATALOG_DIR = ROOT / "catalog"
EXTERNAL_PATH = CATALOG_DIR / "external_discoveries.json"
SUMMARY_PATH = CATALOG_DIR / "catalog_summary.json"
KNOWN_TICS_GZ = CATALOG_DIR / "known_tics.json.gz"

PS_URL = (
    "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?"
    "query=select+pl_name,hostname,tic_id,disc_year,discoverymethod,disc_facility,"
    "pl_orbper,pl_rade,ra,dec,sy_dist+from+pscomppars&format=json"
)
TOI_URL = (
    "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?"
    "query=select+toi,tid,tfopwg_disp,pl_orbper,pl_tranmid,pl_trandep,pl_trandurh,"
    "pl_rade,st_tmag,st_teff,st_rad,st_dist+from+toi&format=json"
)
CTOI_URL = "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv"


def fetch_bytes(url: str, timeout: int = 60) -> bytes:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "open-exoplanet-hunt/1.0 (https://github.com/Th0rgal/open-exoplanet-hunt)"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def normalize_tic(val) -> int | None:
    if val is None:
        return None
    s = str(val).strip()
    if s.upper().startswith("TIC"):
        s = s[3:].strip()
    s = s.split(".")[0].strip()
    if not s or not s.isdigit():
        return None
    return int(s)


def build_catalog() -> dict:
    CATALOG_DIR.mkdir(parents=True, exist_ok=True)
    print("Fetching NASA Exoplanet Archive pscomppars...")
    ps_raw = fetch_bytes(PS_URL)
    ps_rows = json.loads(ps_raw.decode("utf-8"))

    print("Fetching NASA Exoplanet Archive TOI table...")
    toi_raw = fetch_bytes(TOI_URL)
    toi_rows = json.loads(toi_raw.decode("utf-8"))

    print("Fetching ExoFOP-TESS CTOI table...")
    ctoi_raw = fetch_bytes(CTOI_URL)
    ctoi_text = ctoi_raw.decode("utf-8", errors="replace")
    ctoi_reader = csv.DictReader(io.StringIO(ctoi_text))
    ctoi_rows = list(ctoi_reader)

    known_tics: dict[str, dict] = {}

    def ensure_tic(tic_int: int) -> dict:
        key = str(tic_int)
        if key not in known_tics:
            known_tics[key] = {
                "tic_id": tic_int,
                "confirmed_planets": [],
                "tois": [],
                "ctois": [],
                "external_candidates": [],
            }
        return known_tics[key]

    confirmed_with_tic = 0
    for row in ps_rows:
        tic_int = normalize_tic(row.get("tic_id"))
        if tic_int is not None:
            confirmed_with_tic += 1
            entry = ensure_tic(tic_int)
            entry["confirmed_planets"].append({
                "pl_name": row.get("pl_name"),
                "period_d": row.get("pl_orbper"),
                "rade": row.get("pl_rade"),
                "disc_year": row.get("disc_year"),
                "method": row.get("discoverymethod"),
                "facility": row.get("disc_facility"),
            })

    toi_dispositions: dict[str, int] = {}
    for row in toi_rows:
        tic_int = normalize_tic(row.get("tid"))
        disp = str(row.get("tfopwg_disp") or "UNKNOWN")
        toi_dispositions[disp] = toi_dispositions.get(disp, 0) + 1
        if tic_int is not None:
            entry = ensure_tic(tic_int)
            entry["tois"].append({
                "toi": row.get("toi"),
                "disp": disp,
                "period_d": row.get("pl_orbper"),
                "t0_btjd": row.get("pl_tranmid"),
                "depth_ppm": row.get("pl_trandep"),
                "dur_h": row.get("pl_trandurh"),
                "rade": row.get("pl_rade"),
            })

    for row in ctoi_rows:
        tic_int = normalize_tic(row.get("TIC ID") or row.get("tic_id") or row.get("TIC"))
        if tic_int is not None:
            entry = ensure_tic(tic_int)
            period_val = row.get("Period (days)") or row.get("Period")
            try:
                period_f = float(period_val) if period_val else None
            except ValueError:
                period_f = None
            depth_val = row.get("Depth ppm") or row.get("Depth (ppm)")
            try:
                depth_f = float(depth_val) if depth_val else None
            except ValueError:
                depth_f = None
            entry["ctois"].append({
                "ctoi": row.get("CTOI"),
                "disp": row.get("TFOPWG Disposition") or row.get("User Disposition"),
                "period_d": period_f,
                "depth_ppm": depth_f,
            })

    external_data = {"projects": [], "indexed_candidates": []}
    if EXTERNAL_PATH.exists():
        external_data = json.loads(EXTERNAL_PATH.read_text(encoding="utf-8"))
        for cand in external_data.get("indexed_candidates", []):
            tic_int = normalize_tic(cand.get("tic_id"))
            if tic_int is not None:
                entry = ensure_tic(tic_int)
                entry["external_candidates"].append({
                    "name": cand.get("name"),
                    "period_d": cand.get("period_d"),
                    "depth_ppm": cand.get("depth_ppm"),
                    "rade": cand.get("rade"),
                    "source": cand.get("source"),
                    "status": cand.get("status"),
                    "url": cand.get("url"),
                })

    payload_bytes = json.dumps(known_tics, sort_keys=True, separators=(",", ":")).encode("utf-8")
    with gzip.open(KNOWN_TICS_GZ, "wb", compresslevel=9) as gz:
        gz.write(payload_bytes)

    summary = {
        "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "counts": {
            "nasa_confirmed_exoplanets_total": len(ps_rows),
            "nasa_confirmed_with_tic": confirmed_with_tic,
            "nasa_tois_total": len(toi_rows),
            "exofop_ctois_total": len(ctoi_rows),
            "external_indexed_candidates": len(external_data.get("indexed_candidates", [])),
            "external_indexed_projects": len(external_data.get("projects", [])),
            "unique_known_tics_indexed": len(known_tics),
        },
        "toi_dispositions": dict(sorted(toi_dispositions.items())),
        "sources": {
            "nasa_pscomppars_sha256": sha256(ps_raw).hexdigest(),
            "nasa_toi_sha256": sha256(toi_raw).hexdigest(),
            "exofop_ctoi_sha256": sha256(ctoi_raw).hexdigest(),
            "known_tics_gz_sha256": sha256(KNOWN_TICS_GZ.read_bytes()).hexdigest(),
        },
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    build_catalog()
