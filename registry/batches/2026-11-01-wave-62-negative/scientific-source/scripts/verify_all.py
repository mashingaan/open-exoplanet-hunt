#!/usr/bin/env python3
"""Master verification and leaderboard compiler for open-exoplanet-hunt.

Verifies:
1. Assertions are enabled (rejects -O / PYTHONOPTIMIZE).
2. Catalog index integrity (catalog_summary.json, external_discoveries.json, known_tics.json.gz).
3. Every candidate dossier under candidates/TIC-*/certificate.json (all 10 gates, SHA-256 receipts, REPORT.md, diagnostics.png).
4. Scanned target registry integrity (registry/scanned_targets.jsonl).
5. Regenerates registry/leaderboard.json.
"""
import sys
if sys.flags.optimize:
    raise ValueError("Assertions must remain enabled (do not run with -O)")

import gzip
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CATALOG_SUMMARY = ROOT / "catalog/catalog_summary.json"
EXTERNAL_DISC = ROOT / "catalog/external_discoveries.json"
KNOWN_TICS_GZ = ROOT / "catalog/known_tics.json.gz"
CANDIDATES_DIR = ROOT / "candidates"
SCANNED_LOG = ROOT / "registry/scanned_targets.jsonl"
LEADERBOARD_JSON = ROOT / "registry/leaderboard.json"

REQUIRED_GATES = (
    "gate1_catalog_dedup",
    "gate2_multi_sector_baseline",
    "gate3_combined_bls_snr",
    "gate4_leave_one_epoch_out",
    "gate5_odd_even_consistency",
    "gate6_secondary_eclipse_null",
    "gate7_physical_plausibility",
    "gate8_centroid_stability",
    "gate9_per_transit_coherence",
    "gate10_rotational_modulation_guard",
)


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def verify_repository() -> dict:
    require(CATALOG_SUMMARY.exists(), "Missing catalog/catalog_summary.json")
    require(EXTERNAL_DISC.exists(), "Missing catalog/external_discoveries.json")
    require(KNOWN_TICS_GZ.exists(), "Missing catalog/known_tics.json.gz")

    summary = json.loads(CATALOG_SUMMARY.read_text(encoding="utf-8"))
    ext = json.loads(EXTERNAL_DISC.read_text(encoding="utf-8"))
    actual_gz_sha = sha256(KNOWN_TICS_GZ.read_bytes()).hexdigest()
    require(
        summary["sources"]["known_tics_gz_sha256"] == actual_gz_sha,
        "known_tics.json.gz SHA-256 mismatch with catalog_summary.json",
    )

    with gzip.open(KNOWN_TICS_GZ, "rt", encoding="utf-8") as gz:
        known_db = json.load(gz)
    require(
        len(known_db) == summary["counts"]["unique_known_tics_indexed"],
        "Unique TIC count mismatch in known_tics.json.gz",
    )

    verified_candidates = []
    if CANDIDATES_DIR.exists():
        for cdir in sorted(CANDIDATES_DIR.glob("TIC-*")):
            if not cdir.is_dir():
                continue
            cert_path = cdir / "certificate.json"
            report_path = cdir / "REPORT.md"
            plot_path = cdir / "diagnostics.png"
            require(cert_path.exists(), f"Missing {cert_path}")
            require(report_path.exists(), f"Missing {report_path}")
            require(plot_path.exists(), f"Missing {plot_path}")

            cert = json.loads(cert_path.read_text(encoding="utf-8"))
            tic_id = int(cert["tic_id"])
            require(cdir.name == f"TIC-{tic_id}", f"Directory name {cdir.name} != TIC-{tic_id}")
            require(len(cert["sectors_used"]) >= 2, f"TIC {tic_id}: < 2 sectors used")

            sig = cert["primary_signal"]
            gates = sig["gates"]
            for gname in REQUIRED_GATES:
                require(gname in gates, f"TIC {tic_id}: missing gate {gname}")

            all_pass = all(bool(gates[g]["pass"]) for g in REQUIRED_GATES)
            require(bool(sig["all_gates_passed"]) == all_pass, f"TIC {tic_id}: all_gates_passed mismatch")

            status = cert["status"]
            if status in ("unreported_candidate", "candidate_pass_all_10", "reference_benchmark_ddt100"):
                require(all_pass, f"TIC {tic_id} has status {status} but did not pass all 10 gates")

            eph = sig["ephemeris"]
            verified_candidates.append({
                "tic_id": tic_id,
                "status": status,
                "all_gates_passed": all_pass,
                "gates_passed": sum(1 for g in REQUIRED_GATES if gates[g]["pass"]),
                "period_d": eph["period_d"],
                "t0_btjd": eph["t0_btjd"],
                "weighted_dip_ppm": eph["weighted_dip_ppm"],
                "inferred_rp_rearth": eph["inferred_rp_rearth"],
                "combined_bls_snr": eph["combined_bls_snr"],
                "loeo_min_z_sigma": gates["gate4_leave_one_epoch_out"]["min_holdout_z_sigma"],
                "loeo_mean_z_sigma": gates["gate4_leave_one_epoch_out"]["mean_holdout_z_sigma"],
                "sectors": [s["sector"] for s in cert["sectors_used"]],
                "certificate_sha256": sha256(cert_path.read_bytes()).hexdigest(),
            })

    scanned_rows = []
    seen_scanned: set[int] = set()
    if SCANNED_LOG.exists():
        for line_no, line in enumerate(SCANNED_LOG.read_text(encoding="utf-8").splitlines(), start=1):
            if not line.strip():
                continue
            row = json.loads(line)
            tid = int(row["tic_id"])
            require(tid not in seen_scanned, f"Duplicate TIC {tid} in scanned_targets.jsonl at line {line_no}")
            seen_scanned.add(tid)
            scanned_rows.append(row)

    LEADERBOARD_JSON.parent.mkdir(parents=True, exist_ok=True)
    leaderboard = {
        "schema_version": 1,
        "catalog_counts": summary["counts"],
        "external_projects_indexed": len(ext.get("projects", [])),
        "total_targets_scanned_in_repo": len(scanned_rows),
        "verified_dossiers_count": len(verified_candidates),
        "verified_dossiers": sorted(
            verified_candidates,
            key=lambda x: (not x["all_gates_passed"], -x["combined_bls_snr"]),
        ),
    }
    LEADERBOARD_JSON.write_text(json.dumps(leaderboard, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"PASS: verified {summary['counts']['unique_known_tics_indexed']} indexed known TICs, "
        f"{len(verified_candidates)} candidate dossiers, and {len(scanned_rows)} scanned targets."
    )
    return leaderboard


if __name__ == "__main__":
    verify_repository()
