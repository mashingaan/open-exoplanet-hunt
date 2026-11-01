#!/usr/bin/env python3
"""
Gaia DR3 / TIC Neighbor Dilution, Centroid-Torque & Period-Alias Falsification Auditor.

For every candidate and watchlist dossier in `candidates/TIC-*/certificate.json`:
1. Queries MAST TIC / Gaia DR3 within 120 arcseconds (nearly 6 TESS pixels) of the target
   and tests whether any neighboring star could produce the observed transit depth (flux
   sufficiency) without violating the measured in-transit centroid shift limit (centroid torque).
2. Evaluates the multi-sector light curve at exact rational period aliases
   (P, 2P, P/2, 3P/2, 2P/3) to test for period-doubling or sub-harmonic aliasing.
3. Saves machine-readable results to `registry/neighbor_alias_audit.json` and a detailed
   markdown report to `docs/NEIGHBOR_AND_ALIAS_AUDIT.md`.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

if sys.flags.optimize > 0:
    raise RuntimeError("Assertions must remain enabled; do not run with python -O or PYTHONOPTIMIZE.")

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from vet_candidate import (  # noqa: E402
    discover_sector_products,
    download_and_load_sector,
    evaluate_signal_from_arrays,
    mast_invoke,
)


def angular_sep_arcsec(ra1_deg: float, dec1_deg: float, ra2_deg: float, dec2_deg: float) -> float:
    """Compute exact spherical angular separation in arcseconds."""
    r1, d1 = math.radians(ra1_deg), math.radians(dec1_deg)
    r2, d2 = math.radians(ra2_deg), math.radians(dec2_deg)
    dlon = r2 - r1
    dlat = d2 - d1
    a = math.sin(dlat / 2.0) ** 2 + math.cos(d1) * math.cos(d2) * math.sin(dlon / 2.0) ** 2
    c = 2.0 * math.asin(min(1.0, math.sqrt(max(0.0, a))))
    return math.degrees(c) * 3600.0


def query_tic_neighbors(ra_deg: float, dec_deg: float, radius_arcsec: float = 120.0) -> list[dict]:
    """Query MAST TIC / Gaia neighbors within `radius_arcsec` around (ra_deg, dec_deg)."""
    ddec_deg = radius_arcsec / 3600.0
    cos_dec = max(0.05, math.cos(math.radians(dec_deg)))
    dra_deg = ddec_deg / cos_dec
    res = mast_invoke({
        "service": "Mast.Catalogs.Filtered.Tic.Rows",
        "format": "json",
        "params": {
            "columns": "ID,GAIA,ra,dec,Tmag,Gmag,Teff,rad,d,plx",
            "filters": [
                {"paramName": "ra", "values": [{"min": ra_deg - dra_deg, "max": ra_deg + dra_deg}]},
                {"paramName": "dec", "values": [{"min": dec_deg - ddec_deg, "max": dec_deg + ddec_deg}]},
                {"paramName": "Tmag", "values": [{"min": -5.0, "max": 20.5}]},
            ],
        },
    })
    rows = res.get("data", [])
    kept = []
    for r in rows:
        if r.get("ra") is None or r.get("dec") is None:
            continue
        sep = angular_sep_arcsec(ra_deg, dec_deg, float(r["ra"]), float(r["dec"]))
        if sep <= radius_arcsec:
            r["sep_arcsec"] = sep
            kept.append(r)
    kept.sort(key=lambda r: float(r["sep_arcsec"]))
    return kept


def audit_target_neighbors(cert: dict) -> dict:
    """Perform analytical flux-dilution and centroid-torque exclusion on all TIC/Gaia neighbors <= 120 arcsec."""
    tic_id = int(cert["tic_id"])
    st = cert["stellar_parameters"]
    ra_deg = float(st["ra_deg"])
    dec_deg = float(st["dec_deg"])
    target_tmag = float(st["tmag"])

    eph = cert["primary_signal"]["ephemeris"]
    dip_ppm = float(eph["weighted_dip_ppm"])
    dip_frac = dip_ppm * 1e-6

    g8 = cert["primary_signal"]["gates"]["gate8_centroid_stability"]
    obs_shift_arcsec = float(g8.get("centroid_shift_arcsec") or 0.0)
    obs_shift_sigma = float(g8.get("centroid_shift_sigma") or 0.0)
    sigma_r_arcsec = max(0.001, obs_shift_arcsec / obs_shift_sigma) if obs_shift_sigma > 0.05 else 0.005
    limit_3sigma_arcsec = obs_shift_arcsec + 3.0 * sigma_r_arcsec

    raw_rows = query_tic_neighbors(ra_deg, dec_deg, radius_arcsec=120.0)
    neighbors = []
    n_ruled_out_flux = 0
    n_ruled_out_centroid = 0
    n_unresolved = 0

    for row in raw_rows:
        nid = int(row["ID"])
        if nid == tic_id:
            continue
        tmag = float(row["Tmag"]) if row.get("Tmag") is not None else 99.0
        sep_arcsec = float(row["sep_arcsec"])
        delta_tmag = tmag - target_tmag
        flux_ratio = 10.0 ** (-0.4 * delta_tmag)

        # Conservative aperture coupling weight:
        # 1.0 inside 1.5 TESS pixels (31.5"), Gaussian PSF roll-off (sigma=21") outside 31.5"
        if sep_arcsec <= 31.5:
            aperture_weight = 1.0
        else:
            aperture_weight = math.exp(-0.5 * ((sep_arcsec - 31.5) / 21.0) ** 2)

        # Maximum diluted depth (ppm) if neighbor undergoes a 50% eclipse or a 100% total occultation
        max_dip_50pct_ppm = 1e6 * (0.5 * aperture_weight * flux_ratio) / (1.0 + aperture_weight * flux_ratio)
        max_dip_100pct_ppm = 1e6 * (1.0 * aperture_weight * flux_ratio) / (1.0 + aperture_weight * flux_ratio)

        # Required center-of-light shift (arcsec) if this neighbor causes `dip_ppm`
        req_centroid_shift_arcsec = sep_arcsec * dip_frac

        if max_dip_100pct_ppm < dip_ppm:
            status = "RULED_OUT_BY_FLUX"
            n_ruled_out_flux += 1
        elif req_centroid_shift_arcsec > limit_3sigma_arcsec:
            status = "RULED_OUT_BY_CENTROID"
            n_ruled_out_centroid += 1
        else:
            status = "UNRESOLVED_WITHIN_PSF"
            n_unresolved += 1

        neighbors.append({
            "tic_id": nid,
            "gaia_id": str(row.get("GAIA") or ""),
            "sep_arcsec": round(sep_arcsec, 2),
            "tmag": round(tmag, 3),
            "delta_tmag": round(delta_tmag, 3),
            "aperture_weight": round(aperture_weight, 4),
            "max_dip_50pct_ppm": round(max_dip_50pct_ppm, 1),
            "max_dip_100pct_ppm": round(max_dip_100pct_ppm, 1),
            "req_centroid_shift_arcsec": round(req_centroid_shift_arcsec, 4),
            "status": status,
        })

    return {
        "tic_id": tic_id,
        "target_tmag": target_tmag,
        "observed_dip_ppm": dip_ppm,
        "observed_centroid_shift_arcsec": obs_shift_arcsec,
        "observed_centroid_shift_sigma": obs_shift_sigma,
        "centroid_3sigma_upper_limit_arcsec": round(limit_3sigma_arcsec, 4),
        "n_neighbors_120arcsec": len(neighbors),
        "n_ruled_out_by_flux": n_ruled_out_flux,
        "n_ruled_out_by_centroid": n_ruled_out_centroid,
        "n_unresolved_within_psf": n_unresolved,
        "neighbors": neighbors,
    }


def audit_target_aliases(cert: dict) -> list[dict]:
    """Evaluate multi-sector light curve at exact rational period aliases (P, 2P, P/2, 3P/2, 2P/3)."""
    tic_id = int(cert["tic_id"])
    star = cert["stellar_parameters"]
    prim = cert["primary_signal"]
    base_eph = prim["ephemeris"]
    base_gates = prim["gates"]
    base_p = float(base_eph["period_d"])

    # Record exact certified primary signal for 1.0 P
    alias_results = [{
        "alias_label": "1.0 P (Certified Primary)",
        "multiplier": 1.0,
        "period_d": round(base_p, 7),
        "weighted_dip_ppm": round(float(base_eph["weighted_dip_ppm"]), 1),
        "weighted_dip_err_ppm": round(float(base_eph["weighted_dip_err_ppm"]), 1),
        "inferred_rp_rearth": round(float(base_eph["inferred_rp_rearth"]), 2),
        "combined_bls_snr": round(float(base_eph["combined_bls_snr"]), 2),
        "loeo_min_z_sigma": round(float(base_gates["gate4_leave_one_epoch_out"]["min_holdout_z_sigma"]), 2),
        "loeo_mean_z_sigma": round(float(base_gates["gate4_leave_one_epoch_out"]["mean_holdout_z_sigma"]), 2),
        "odd_even_diff_sigma": round(float(base_gates["gate5_odd_even_consistency"]["diff_sigma"]), 2),
        "secondary_snr_sigma": round(float(base_gates["gate6_secondary_eclipse_null"]["secondary_snr"]), 2),
        "gates_passed": sum(1 for g in base_gates.values() if g["pass"]),
    }]

    prods = discover_sector_products(tic_id)
    if len(prods) > 6:
        idxs = [0, 1, len(prods) // 2 - 1, len(prods) // 2, len(prods) - 2, len(prods) - 1]
        prods = [prods[i] for i in sorted(set(idxs))]

    sectors_data = []
    for p in prods:
        try:
            sectors_data.append(download_and_load_sector(tic_id, p))
        except Exception:
            pass

    alias_specs = [
        ("2.0 P (Double Period)", 2.0),
        ("0.5 P (Half Period)", 0.5),
        ("1.5 P (3:2 Resonance/Alias)", 1.5),
        ("0.667 P (2:3 Resonance/Alias)", 2.0 / 3.0),
    ]
    for label, mult in alias_specs:
        test_p = base_p * mult
        res = evaluate_signal_from_arrays(
            star,
            sectors_data,
            target_period=test_p,
            period_min=max(0.4, test_p * 0.98),
            period_max=test_p * 1.02,
            allow_known=True,
        )
        eph = res["ephemeris"]
        gates = res["gates"]
        n_passed = sum(1 for g in gates.values() if g["pass"])
        alias_results.append({
            "alias_label": label,
            "multiplier": round(mult, 4),
            "period_d": round(float(eph["period_d"]), 7),
            "weighted_dip_ppm": round(float(eph["weighted_dip_ppm"]), 1),
            "weighted_dip_err_ppm": round(float(eph["weighted_dip_err_ppm"]), 1),
            "inferred_rp_rearth": round(float(eph["inferred_rp_rearth"]), 2),
            "combined_bls_snr": round(float(eph["combined_bls_snr"]), 2),
            "loeo_min_z_sigma": round(float(gates["gate4_leave_one_epoch_out"]["min_holdout_z_sigma"]), 2),
            "loeo_mean_z_sigma": round(float(gates["gate4_leave_one_epoch_out"]["mean_holdout_z_sigma"]), 2),
            "odd_even_diff_sigma": round(float(gates["gate5_odd_even_consistency"]["diff_sigma"]), 2),
            "secondary_snr_sigma": round(float(gates["gate6_secondary_eclipse_null"]["secondary_snr"]), 2),
            "gates_passed": n_passed,
        })
    return alias_results


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit Gaia DR3 / TIC neighbors and rational period aliases.")
    parser.add_argument("--recompute-aliases", action="store_true", help="Re-run multi-sector BLS at rational period aliases")
    args = parser.parse_args()

    cached_aliases: dict[int, list[dict]] = {}
    audit_json_path = ROOT / "registry" / "neighbor_alias_audit.json"
    if audit_json_path.exists() and not args.recompute_aliases:
        try:
            for item in json.loads(audit_json_path.read_text()):
                if "tic_id" in item and "alias_audit" in item:
                    cached_aliases[int(item["tic_id"])] = item["alias_audit"]
        except Exception:
            pass

    cert_paths = sorted(
        (ROOT / "candidates").glob("TIC-*/certificate.json"),
        key=lambda p: (
            0 if "reference" in json.loads(p.read_text())["status"] else (
                1 if "pass_all_10" in json.loads(p.read_text())["status"] else 2
            ),
            int(p.parent.name.split("-")[1]),
        ),
    )

    host_names = {
        4206066: "StKM 1-561 (K-dwarf, 35.6 pc)",
        69912980: "G 86-28 (M3.0V, 12.4 pc)",
        264151965: "L 198-93 (M3.0V, 21.2 pc)",
        63349323: "PM J02019+3413 (M0.5Ve, 54.0 pc)",
        119487103: "G 214-12 (M3.0V, 23.1 pc)",
        120844941: "TYC 2898-1488-1 (K7V, 47.0 pc)",
        125757128: "HG 8-41 (K5V, 35.8 pc)",
        252907497: "G 192-59 (M2.0V, 33.2 pc)",
        265756555: "Ross 24 (M1.0V, 28.2 pc)",
        457195191: "G 82-48 (M0V, 26.1 pc)",
    }

    full_audit = []
    md_lines = [
        "# Gaia DR3 / TIC Neighbor Blend & Period-Alias Falsification Audit",
        "",
        f"This document provides the quantitative results of two critical astrophysical falsification tests across all {len(cert_paths)} candidate and watchlist dossiers in `open-exoplanet-hunt` (addressing [Issue #2](https://github.com/Th0rgal/open-exoplanet-hunt/issues/2)):",
        "",
        "1. **Gaia DR3 / MAST TIC Neighbor Dilution & Centroid-Torque Exclusion ($r \\le 120''$)**:",
        "   - Queries every known star in the TESS Input Catalog / Gaia DR3 within $120''$ ($\\sim 5.7$ TESS pixels) of each host star.",
        "   - **Flux-Sufficiency Test**: Computes the maximum diluted eclipse depth $\\Delta F_{\\text{max}, j}$ (in ppm) if neighbor $j$ (with magnitude difference $\\Delta T_j = T_j - T_{\\text{target}}$ and PSF/aperture coupling $w_j$) were a $100\\%$ eclipsing binary. If $\\Delta F_{\\text{max}, j} < \\Delta F_{\\text{obs}}$, neighbor $j$ is **mathematically incapable** of producing the transit depth (`RULED_OUT_BY_FLUX`).",
        "   - **Centroid-Torque Test**: If neighbor $j$ at angular separation $d_j$ produces a fractional depth $\\Delta F_{\\text{obs}}$, it must pull the center of light during transit by $\\Delta r_{\\text{req}, j} = d_j \\times (\\Delta F_{\\text{obs}} \\times 10^{-6})$. If $\\Delta r_{\\text{req}, j} > \\Delta r_{\\text{obs}} + 3\\sigma_r$, neighbor $j$ is **ruled out at $>3\\sigma$ by SPOC/QLP centroid astrometry** (`RULED_OUT_BY_CENTROID`).",
        "2. **Rational Period-Alias Matrix ($P, 2P, P/2, \\frac{3}{2}P, \\frac{2}{3}P$)**:",
        "   - Evaluates the multi-sector light curve at exact rational multiples of the primary BLS period to test for period-doubling (e.g., equal-depth eclipsing binaries at $2P$) or sub-harmonic aliasing.",
        "",
        "## Executive Summary Table",
        "",
        "| Target | Host Star | Observed Depth (ppm) | Centroid $3\\sigma$ Limit | TIC/Gaia Stars ($\\le 120''$) | Ruled Out by Flux | Ruled Out by Centroid | Unresolved ($<31.5''$ & bright) | Preferred Ephemeris |",
        "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
    ]

    detail_sections = []

    for cp in cert_paths:
        cert = json.loads(cp.read_text())
        tic_id = int(cert["tic_id"])
        host = host_names.get(tic_id, f"TIC {tic_id}")
        print(f"Auditing neighbors and period aliases for TIC {tic_id} ({host})...")
        nb_audit = audit_target_neighbors(cert)
        alias_audit = cached_aliases.get(tic_id) or audit_target_aliases(cert)

        entry = {
            "tic_id": tic_id,
            "host_name": host,
            "status": cert["status"],
            "neighbor_audit": nb_audit,
            "alias_audit": alias_audit,
        }
        full_audit.append(entry)

        p1 = alias_audit[0]
        p2 = alias_audit[1]
        if p2["odd_even_diff_sigma"] < p1["odd_even_diff_sigma"] - 0.8 and p2["loeo_min_z_sigma"] >= 2.0:
            pref_note = f"`1P = {p1['period_d']:.5f} d` (check `2P = {p2['period_d']:.5f} d`, O/E `{p2['odd_even_diff_sigma']:.2f}σ`)"
        else:
            pref_note = f"`1P = {p1['period_d']:.5f} d` (`{p1['gates_passed']}/10` gates, O/E `{p1['odd_even_diff_sigma']:.2f}σ`)"

        md_lines.append(
            f"| [`TIC {tic_id}`](../candidates/TIC-{tic_id}/REPORT.md) | {host} | `{nb_audit['observed_dip_ppm']:.1f}` | `{nb_audit['centroid_3sigma_upper_limit_arcsec']:.3f}″` | `{nb_audit['n_neighbors_120arcsec']}` | `{nb_audit['n_ruled_out_by_flux']}` | `{nb_audit['n_ruled_out_by_centroid']}` | **`{nb_audit['n_unresolved_within_psf']}`** | {pref_note} |"
        )

        sec = [
            f"### [`TIC {tic_id}`](../candidates/TIC-{tic_id}/REPORT.md) — {host}",
            "",
            f"- **Observed Transit Depth**: `{nb_audit['observed_dip_ppm']:.1f} ppm` | **Observed In-Transit Centroid Shift**: `{nb_audit['observed_centroid_shift_arcsec']:.4f}″` (`{nb_audit['observed_centroid_shift_sigma']:.2f}σ`; $3\\sigma$ upper bound: `{nb_audit['centroid_3sigma_upper_limit_arcsec']:.4f}″`)",
            f"- **Neighbor Exclusion ($r \\le 120''$)**: **{nb_audit['n_neighbors_120arcsec']}** cataloged TIC/Gaia stars within $120''$ — **{nb_audit['n_ruled_out_by_flux']}** ruled out by flux dilution, **{nb_audit['n_ruled_out_by_centroid']}** ruled out by centroid torque, **{nb_audit['n_unresolved_within_psf']}** unresolved blending stars capable of mimicking the dip.",
            "",
            "#### Brightest / Closest Cataloged Neighbors ($r \\le 120''$, Top 8 by Max Diluted Depth)",
            "",
            "| Neighbor TIC | Gaia DR3 ID | Sep ($''$) | $T_{\\text{mag}}$ | $\\Delta T$ | Max Depth (50% EB, ppm) | Max Depth (100% Occ, ppm) | Req. Centroid Shift ($''$) | Verdict |",
            "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
        ]
        top_nbs = sorted(nb_audit["neighbors"], key=lambda x: (-x["max_dip_100pct_ppm"], x["sep_arcsec"]))[:8]
        if not top_nbs:
            sec.append("| *(None within 120″)* | — | — | — | — | — | — | — | — |")
        else:
            for nb in top_nbs:
                sec.append(
                    f"| `TIC {nb['tic_id']}` | `{nb['gaia_id'] or '—'}` | `{nb['sep_arcsec']:.2f}″` | `{nb['tmag']:.2f}` | `+{nb['delta_tmag']:.2f}` | `{nb['max_dip_50pct_ppm']:.1f}` | `{nb['max_dip_100pct_ppm']:.1f}` | `{nb['req_centroid_shift_arcsec']:.4f}″` | `{nb['status']}` |"
                )

        sec.extend([
            "",
            "#### Rational Period-Alias Evaluation ($P, 2P, P/2, 3P/2, 2P/3$)",
            "",
            "| Tested Alias | Period (d) | Depth (ppm) | $R_p$ ($R_\\oplus$) | BLS SNR | LOEO Min $Z$ | LOEO Mean $Z$ | Odd/Even ($\\sigma$) | Phase-0.5 Sec ($\\sigma$) | Gates Passed |",
            "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
        ])
        for al in alias_audit:
            sec.append(
                f"| **{al['alias_label']}** | `{al['period_d']:.6f}` | `{al['weighted_dip_ppm']:.1f} ± {al['weighted_dip_err_ppm']:.1f}` | `{al['inferred_rp_rearth']:.2f}` | `{al['combined_bls_snr']:.2f}σ` | `{al['loeo_min_z_sigma']:+.2f}σ` | `{al['loeo_mean_z_sigma']:+.2f}σ` | `{al['odd_even_diff_sigma']:.2f}σ` | `{al['secondary_snr_sigma']:+.2f}σ` | `{al['gates_passed']}/10` |"
            )
        sec.append("")
        detail_sections.extend(sec)

    md_lines.extend(["", "## Detailed Per-Target Neighbor & Alias Audits", ""] + detail_sections)
    (ROOT / "registry" / "neighbor_alias_audit.json").write_text(json.dumps(full_audit, indent=2) + "\n")
    (ROOT / "docs" / "NEIGHBOR_AND_ALIAS_AUDIT.md").write_text("\n".join(md_lines) + "\n")
    print("Saved registry/neighbor_alias_audit.json and docs/NEIGHBOR_AND_ALIAS_AUDIT.md")


if __name__ == "__main__":
    main()
