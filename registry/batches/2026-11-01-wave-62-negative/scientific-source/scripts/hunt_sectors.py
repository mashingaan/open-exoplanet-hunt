#!/usr/bin/env python3
"""Autonomous multi-sector TESS K/M-dwarf exoplanet hunter.

Queries the MAST TESS Input Catalog for nearby K and M dwarfs without any known
NASA exoplanet, TOI, or ExoFOP CTOI, stitches multi-year SPOC / TESS-SPOC / QLP
light curves across all available sectors, and runs the 10-gate adversarial
falsification pipeline.
"""
import sys
if sys.flags.optimize:
    raise ValueError("Assertions must remain enabled (do not run with -O)")

import argparse
import gzip
import json
from pathlib import Path
import time

from vet_candidate import (
    KNOWN_TICS_GZ,
    ROOT,
    discover_sector_products,
    download_and_load_sector,
    evaluate_signal_from_arrays,
     generate_markdown_report,
    mast_invoke,
    render_diagnostic_plot,
)
from hashlib import sha256

REGISTRY_DIR = ROOT / "registry"
SCANNED_LOG = REGISTRY_DIR / "scanned_targets.jsonl"


def load_known_and_scanned() -> tuple[set[int], set[int]]:
    known: set[int] = set()
    if KNOWN_TICS_GZ.exists():
        with gzip.open(KNOWN_TICS_GZ, "rt", encoding="utf-8") as gz:
            db = json.load(gz)
            known = {int(k) for k in db}
    scanned: set[int] = set()
    if SCANNED_LOG.exists():
        for line in SCANNED_LOG.read_text(encoding="utf-8").splitlines():
            if line.strip():
                row = json.loads(line)
                scanned.add(int(row["tic_id"]))
    return known, scanned


def query_candidate_km_dwarfs(
    tmag_min: float = 9.0,
    tmag_max: float = 11.2,
    teff_min: float = 3200.0,
    teff_max: float = 4300.0,
    rad_max: float = 0.70,
    dist_max: float = 60.0,
    page_size: int = 200,
    page: int = 1,
) -> list[dict]:
    cache_dir = ROOT / ".cache/mast"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"tic_page_{page}_{page_size}.json"
    if cache_file.exists():
        rows = json.loads(cache_file.read_text(encoding="utf-8"))
    else:
        res = mast_invoke({
            "service": "Mast.Catalogs.Filtered.Tic.Rows",
            "format": "json",
            "params": {
                "columns": "ID,ra,dec,Tmag,Teff,rad,mass,d,plx,lum,logg,contratio",
                "filters": [
                    {"paramName": "Tmag", "values": [{"min": tmag_min, "max": tmag_max}]},
                    {"paramName": "Teff", "values": [{"min": teff_min, "max": teff_max}]},
                    {"paramName": "rad", "values": [{"min": 0.20, "max": rad_max}]},
                    {"paramName": "d", "values": [{"min": 5.0, "max": dist_max}]},
                    {"paramName": "objType", "values": ["STAR"]},
                ],
            },
            "pagesize": page_size,
            "page": page,
        }, timeout=60, retries=3)
        rows = res.get("data", [])
        cache_file.write_text(json.dumps(rows), encoding="utf-8")
    out = []
    for r in rows:
        if r.get("ID") is None or r.get("rad") is None or r.get("Teff") is None:
            continue
        out.append({
            "tic_id": int(r["ID"]),
            "ra_deg": float(r["ra"]) if r.get("ra") is not None else None,
            "dec_deg": float(r["dec"]) if r.get("dec") is not None else None,
            "tmag": float(r["Tmag"]) if r.get("Tmag") is not None else None,
            "teff_k": float(r["Teff"]),
            "rad_rsun": float(r["rad"]),
            "mass_msun": float(r["mass"]) if r.get("mass") is not None else float(r["rad"]),
            "dist_pc": float(r["d"]) if r.get("d") is not None else None,
            "plx_mas": float(r["plx"]) if r.get("plx") is not None else None,
            "contratio": float(r["contratio"]) if r.get("contratio") is not None else None,
        })
    return out


def append_scan_log(entry: dict) -> None:
    REGISTRY_DIR.mkdir(parents=True, exist_ok=True)
    with SCANNED_LOG.open("a", encoding="utf-8") as f:
        f.write(json.dumps(entry, sort_keys=True) + "\n")


def run_hunt(
    max_targets: int = 12,
    min_sectors: int = 3,
    tmag_min: float = 9.2,
    tmag_max: float = 11.0,
    teff_min: float = 3300.0,
    teff_max: float = 4300.0,
    dist_max: float = 55.0,
    page: int = 1,
) -> list[dict]:
    known, scanned = load_known_and_scanned()
    results = []
    processed = 0
    cur_page = page
    while processed < max_targets and cur_page < page + 10:
        stars = query_candidate_km_dwarfs(
            tmag_min=tmag_min,
            tmag_max=tmag_max,
            teff_min=teff_min,
            teff_max=teff_max,
            dist_max=dist_max,
            page=cur_page,
        )
        if not stars:
            break
        print(f"Queried {len(stars)} nearby K/M dwarfs from MAST TIC (page={cur_page}).")
        found_on_page = 0
        for star in stars:
            if processed >= max_targets:
                break
            tic_id = star["tic_id"]
            if tic_id in known or tic_id in scanned:
                continue

            try:
                prods = discover_sector_products(tic_id)
            except Exception as exc:
                print(f"  [skip] TIC {tic_id}: MAST product query error: {exc}")
                continue

            if len(prods) < min_sectors:
                continue
            found_on_page += 1

            all_prods = list(prods)
            # Limit initial fast screen to at most 6 well-separated sectors; expand to all sectors if promising
            if len(prods) > 6:
                idxs = [0, 1, len(prods) // 2 - 1, len(prods) // 2, len(prods) - 2, len(prods) - 1]
                prods = [prods[i] for i in sorted(set(idxs))]

            print(f"[{processed + 1}/{max_targets}] Scanning TIC {tic_id} (Tmag={star['tmag']:.2f}, Teff={star['teff_k']:.0f}K, R={star['rad_rsun']:.2f}Rsun, d={star['dist_pc']:.1f}pc, sectors={[p['sector'] for p in prods]})...")
            try:
                sectors_data = []
                for p in prods:
                    try:
                        sectors_data.append(download_and_load_sector(tic_id, p))
                    except Exception as sec_exc:
                        print(f"  [skip-sector] TIC {tic_id} s{p['sector']}: {sec_exc}")
                if len(sectors_data) < min_sectors:
                    scanned.add(tic_id)
                    continue
                eval_res = evaluate_signal_from_arrays(
                    star,
                    sectors_data,
                    target_period=None,
                    period_min=1.1,
                    period_max=15.0,
                    allow_known=False,
                )
            except Exception as exc:
                print(f"  [error] TIC {tic_id}: {exc}")
                continue

            def _compute_verdict(res_dict: dict) -> tuple[str, int, list[str]]:
                e_d = res_dict["ephemeris"]
                g_d = res_dict["gates"]
                n_p = sum(1 for g in g_d.values() if g["pass"])
                f_g = [k for k, g in g_d.items() if not g["pass"]]
                v_str = "CANDIDATE_PASS_ALL_10" if res_dict["all_gates_passed"] else (
                    "WATCHLIST_8_PLUS_GATES"
                    if (
                        n_p >= 9
                        and g_d["gate4_leave_one_epoch_out"]["min_holdout_z_sigma"] >= 1.5
                        and g_d["gate4_leave_one_epoch_out"]["mean_holdout_z_sigma"] >= 3.5
                        and g_d["gate5_odd_even_consistency"]["diff_sigma"] < 4.0
                        and abs(g_d["gate6_secondary_eclipse_null"]["secondary_snr"]) < 4.0
                        and g_d["gate7_physical_plausibility"]["pass"]
                        and (
                            g_d["gate9_per_transit_coherence"]["pass"]
                            or (
                                g_d["gate9_per_transit_coherence"]["positive_fraction"] >= 0.80
                                and g_d["gate9_per_transit_coherence"]["jackknife_snr"] >= 6.0
                            )
                        )
                        and g_d["gate10_rotational_modulation_guard"]["pass"]
                        and e_d["combined_bls_snr"] >= 7.5
                    )
                    else "REJECTED_BY_FALSIFICATION"
                )
                return v_str, n_p, f_g

            verdict, n_passed, failed_gates = _compute_verdict(eval_res)
            if verdict in ("CANDIDATE_PASS_ALL_10", "WATCHLIST_8_PLUS_GATES") and len(all_prods) > len(prods):
                loaded_sec_nums = {s["sector"] for s in sectors_data}
                extra_prods = [p for p in all_prods if p["sector"] not in loaded_sec_nums]
                for p in extra_prods:
                    try:
                        sectors_data.append(download_and_load_sector(tic_id, p))
                    except Exception as sec_exc:
                        print(f"  [skip-sector] TIC {tic_id} s{p['sector']}: {sec_exc}")
                sectors_data.sort(key=lambda s: s["sector"])
                print(f"  [full-sector-audit] Re-evaluating TIC {tic_id} across all {len(sectors_data)} sectors {[s['sector'] for s in sectors_data]}...")
                eval_res = evaluate_signal_from_arrays(
                    star,
                    sectors_data,
                    target_period=None,
                    period_min=1.1,
                    period_max=15.0,
                    allow_known=False,
                )
                verdict, n_passed, failed_gates = _compute_verdict(eval_res)

            processed += 1
            scanned.add(tic_id)
            eph = eval_res["ephemeris"]
            gates = eval_res["gates"]

            log_entry = {
                "scanned_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "tic_id": tic_id,
                "tmag": star["tmag"],
                "teff_k": star["teff_k"],
                "rad_rsun": star["rad_rsun"],
                "dist_pc": star["dist_pc"],
                "sectors": [s["sector"] for s in sectors_data],
                "best_period_d": eph["period_d"],
                "t0_btjd": eph["t0_btjd"],
                "weighted_dip_ppm": eph["weighted_dip_ppm"],
                "inferred_rp_rearth": eph["inferred_rp_rearth"],
                "combined_bls_snr": eph["combined_bls_snr"],
                "loeo_min_z": gates["gate4_leave_one_epoch_out"]["min_holdout_z_sigma"],
                "loeo_mean_z": gates["gate4_leave_one_epoch_out"]["mean_holdout_z_sigma"],
                "gates_passed": n_passed,
                "failed_gates": failed_gates,
                "verdict": verdict,
            }
            append_scan_log(log_entry)
            results.append(log_entry)
            print(
                f"  -> {verdict} ({n_passed}/10 gates) | P={eph['period_d']:.6f} d | "
                f"dip={eph['weighted_dip_ppm']:.0f} ppm (Rp={eph['inferred_rp_rearth']:.2f} R⊕) | "
                f"SNR={eph['combined_bls_snr']:.2f} | LOEO mean Z={gates['gate4_leave_one_epoch_out']['mean_holdout_z_sigma']:.2f}σ | "
                f"failed={failed_gates}"
            )

            # Save full dossier if it passes all 10 gates or is a high-priority watchlist target
            if verdict in ("CANDIDATE_PASS_ALL_10", "WATCHLIST_8_PLUS_GATES"):
                out_dir = ROOT / f"candidates/TIC-{tic_id}"
                out_dir.mkdir(parents=True, exist_ok=True)
                render_diagnostic_plot(star, eval_res, out_dir / "diagnostics.png")
                primary_clean = {k: v for k, v in eval_res.items() if k != "_arrays"}
                cert = {
                    "schema_version": 1,
                    "tic_id": tic_id,
                    "status": verdict.lower(),
                    "generated_utc": log_entry["scanned_utc"],
                    "stellar_parameters": star,
                    "sectors_used": [
                        {
                            "sector": s["sector"],
                            "pipeline": s["pipeline"],
                            "exptime_s": s["exptime_s"],
                            "uri": s["uri"],
                            "flux_col": s["flux_col"],
                            "sha256": s["sha256"],
                            "n_points": s["n_points"],
                            "mad_ppm": round(s["mad_ppm"], 1),
                        }
                        for s in sectors_data
                    ],
                    "primary_signal": primary_clean,
                    "secondary_candidates": [],
                    "script_sha256": sha256((ROOT / "scripts/vet_candidate.py").read_bytes()).hexdigest(),
                }
                (out_dir / "certificate.json").write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n", encoding="utf-8")
                (out_dir / "REPORT.md").write_text(generate_markdown_report(cert), encoding="utf-8")

        cur_page += 1

    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Scan unreported multi-sector TESS K/M dwarfs.")
    parser.add_argument("--max-targets", type=int, default=10, help="Number of multi-sector stars to scan in this batch")
    parser.add_argument("--min-sectors", type=int, default=3, help="Minimum TESS sectors required")
    parser.add_argument("--tmag-min", type=float, default=9.2)
    parser.add_argument("--tmag-max", type=float, default=11.0)
    parser.add_argument("--teff-min", type=float, default=3300.0)
    parser.add_argument("--teff-max", type=float, default=4300.0)
    parser.add_argument("--dist-max", type=float, default=55.0)
    parser.add_argument("--page", type=int, default=1)
    args = parser.parse_args()

    run_hunt(
        max_targets=args.max_targets,
        min_sectors=args.min_sectors,
        tmag_min=args.tmag_min,
        tmag_max=args.tmag_max,
        teff_min=args.teff_min,
        teff_max=args.teff_max,
        dist_max=args.dist_max,
        page=args.page,
    )


if __name__ == "__main__":
    main()
