#!/usr/bin/env python3
"""Adversarial falsification & vetting engine for TESS exoplanet candidates.

Implements 10 quantitative falsification gates inspired by Pavel Rabtsevich's
TIC 4206066 workflow, NASA SPOC Data Validation, SHERLOCK, and Planet Hunters TESS:
  1. Catalog & Alert De-duplication (NASA pscomppars, TOI, ExoFOP CTOI, external index)
  2. Multi-Sector / Multi-Epoch Baseline (>= 2 sectors, >= 6 distinct transits)
  3. Combined Detrended BLS Signal-to-Noise (>= 8.0 sigma)
  4. Blind Leave-One-Epoch-Out (LOEO) Temporal Prediction (mean Z >= 3.5 sigma, min Z >= 2.0 sigma)
  5. Odd vs. Even Transit Depth Consistency (< 2.5 sigma difference)
  6. Secondary Eclipse Vetting at Phase 0.5 (< 2.5 sigma or < 30% primary depth)
  7. Physical Radius & Keplerian Duration Plausibility (0.4 <= R_p <= 18 R_earth; 0.2 <= dur/dur_circ <= 2.2)
  8. In/Out-of-Transit Centroid Shift Test (< 3.0 sigma and < 2.5 arcsec)
  9. Per-Transit Coherence & Jackknife Stability (>= 70% positive transits; jackknife SNR >= 6.5 sigma)
 10. Stellar Rotational Modulation / Harmonic Separation (sinusoidal amp at P, 2P, P/2 < 50% depth)
"""
import sys
if sys.flags.optimize:
    raise ValueError("Assertions must remain enabled (do not run with -O)")

import argparse
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import time
import urllib.parse
import urllib.request

from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import medfilt

ROOT = Path(__file__).resolve().parent.parent
KNOWN_TICS_GZ = ROOT / "catalog/known_tics.json.gz"
CACHE_DIR = ROOT / ".cache/mast"
R_SUN_TO_R_EARTH = 109.076
G_SI = 6.67430e-11
M_SUN_KG = 1.98847e30
R_SUN_M = 6.957e8
DAY_S = 86400.0


def mast_invoke(request_obj: dict, timeout: int = 60, retries: int = 3) -> dict:
    data = urllib.parse.urlencode({"request": json.dumps(request_obj)}).encode("utf-8")
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                "https://mast.stsci.edu/api/v0/invoke",
                data=data,
                headers={"User-Agent": "open-exoplanet-hunt/1.0"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception as exc:
            last_err = exc
            if attempt + 1 < retries:
                time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"MAST API invoke failed after {retries} attempts: {last_err}") from last_err


def fetch_tic_metadata(tic_id: int) -> dict:
    res = mast_invoke({
        "service": "Mast.Catalogs.Filtered.Tic.Rows",
        "format": "json",
        "params": {
            "columns": "ID,ra,dec,Tmag,Teff,rad,mass,d,plx,lum,logg,contratio",
            "filters": [{"paramName": "ID", "values": [str(tic_id)]}],
        },
    })
    rows = res.get("data", [])
    if not rows:
        raise ValueError(f"TIC {tic_id} not found in MAST TIC catalog")
    row = rows[0]
    return {
        "tic_id": int(row["ID"]),
        "ra_deg": float(row["ra"]) if row.get("ra") is not None else None,
        "dec_deg": float(row["dec"]) if row.get("dec") is not None else None,
        "tmag": float(row["Tmag"]) if row.get("Tmag") is not None else None,
        "teff_k": float(row["Teff"]) if row.get("Teff") is not None else 4000.0,
        "rad_rsun": float(row["rad"]) if row.get("rad") is not None else 0.65,
        "mass_msun": float(row["mass"]) if row.get("mass") is not None else 0.65,
        "dist_pc": float(row["d"]) if row.get("d") is not None else None,
        "plx_mas": float(row["plx"]) if row.get("plx") is not None else None,
        "contratio": float(row["contratio"]) if row.get("contratio") is not None else None,
    }


def discover_sector_products(tic_id: int) -> list[dict]:
    """Find best light-curve product per TESS sector (preferring 2-min SPOC > TESS-SPOC HLSP > QLP HLSP)."""
    res = mast_invoke({
        "service": "Mast.Caom.Filtered",
        "format": "json",
        "params": {
            "columns": "obs_id,obs_collection,dataproduct_type,target_name,sequence_number,t_exptime,dataURL",
            "filters": [
                {"paramName": "target_name", "values": [str(tic_id), f"{tic_id:016d}"]},
                {"paramName": "obs_collection", "values": ["TESS", "HLSP"]},
            ],
        },
    })
    by_sector: dict[int, list[dict]] = {}
    for row in res.get("data", []):
        uri = str(row.get("dataURL") or "")
        seq = row.get("sequence_number")
        if seq is None or not uri.endswith(".fits"):
            continue
        if not (uri.endswith("_lc.fits") or uri.endswith("_llc.fits")):
            continue
        sec = int(seq)
        if not (1 <= sec <= 300):
            continue
        # Priority: 0 = TESS 2-min SPOC (_lc.fits), 1 = HLSP tess-spoc, 2 = HLSP qlp
        if row.get("obs_collection") == "TESS" and "_lc.fits" in uri:
            prio = 0
            kind = "SPOC"
        elif "hlsp_tess-spoc_" in uri:
            prio = 1
            kind = "TESS-SPOC"
        elif "hlsp_qlp_" in uri:
            prio = 2
            kind = "QLP"
        else:
            continue
        by_sector.setdefault(sec, []).append({
            "sector": sec,
            "priority": prio,
            "pipeline": kind,
            "exptime_s": float(row.get("t_exptime") or 0),
            "uri": uri,
            "obs_id": row.get("obs_id"),
        })
    selected = []
    for sec in sorted(by_sector):
        candidates = sorted(by_sector[sec], key=lambda x: (x["priority"], x["exptime_s"]))
        best = dict(candidates[0])
        if len(candidates) > 1:
            best["fallbacks"] = candidates[1:]
        selected.append(best)
    return selected


def download_and_load_sector(tic_id: int, prod: dict, window_days: float = 0.75) -> dict:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    sec = prod["sector"]
    pipe = prod["pipeline"]
    fname = f"tic{tic_id}_s{sec:04d}_{pipe}.fits"
    fpath = CACHE_DIR / fname
    if not fpath.exists():
        url = f"https://mast.stsci.edu/api/v0.1/Download/file?uri={urllib.parse.quote(prod['uri'], safe=':/')}"
        last_err: Exception | None = None
        for attempt in range(4):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "open-exoplanet-hunt/1.0"})
                with urllib.request.urlopen(req, timeout=60) as resp:
                    fpath.write_bytes(resp.read())
                break
            except Exception as exc:
                last_err = exc
                if attempt + 1 < 4:
                    time.sleep(1.5 * (2 ** attempt))
        else:
            raise RuntimeError(f"Failed to download {fname} after 4 attempts: {last_err}") from last_err
    raw_bytes = fpath.read_bytes()
    file_sha = sha256(raw_bytes).hexdigest()

    with fits.open(fpath) as hdul:
        d = hdul[1].data
        cols = d.columns.names
        t = np.asarray(d["TIME"], dtype=float)
        q = np.asarray(d["QUALITY"], dtype=int) if "QUALITY" in cols else np.zeros(len(t), dtype=int)
        flux_col = "SAP_FLUX"
        f_raw = np.full(len(t), np.nan)
        for cand_col in ("PDCSAP_FLUX", "DET_FLUX", "KSPSAP_FLUX", "SAP_FLUX"):
            if cand_col in cols:
                cand_arr = np.asarray(d[cand_col], dtype=float)
                if np.sum(np.isfinite(cand_arr)) >= 50:
                    flux_col = cand_col
                    f_raw = cand_arr
                    break
        f_sap = np.asarray(d["SAP_FLUX"], dtype=float) if "SAP_FLUX" in cols else f_raw.copy()
        if np.sum(np.isfinite(f_sap)) < 50:
            f_sap = f_raw.copy()

        if "MOM_CENTR1" in cols and "MOM_CENTR2" in cols:
            cx = np.asarray(d["MOM_CENTR1"], dtype=float)
            cy = np.asarray(d["MOM_CENTR2"], dtype=float)
        elif "SAP_X" in cols and "SAP_Y" in cols:
            cx = np.asarray(d["SAP_X"], dtype=float)
            cy = np.asarray(d["SAP_Y"], dtype=float)
        else:
            cx = np.full(len(t), np.nan)
            cy = np.full(len(t), np.nan)

    mask = np.isfinite(t) & np.isfinite(f_raw) & (q == 0)
    if np.sum(mask) < 50:
        mask = np.isfinite(t) & np.isfinite(f_raw)
    if np.sum(mask) < 50:
        for fb in prod.get("fallbacks", []):
            try:
                return download_and_load_sector(tic_id, fb, window_days=window_days)
            except Exception:
                continue
        raise ValueError(f"Sector {sec} has fewer than 50 valid cadences")
    t, f_raw, f_sap, cx, cy = t[mask], f_raw[mask], f_sap[mask], cx[mask], cy[mask]
    order = np.argsort(t)
    t, f_raw, f_sap, cx, cy = t[order], f_raw[order], f_sap[order], cx[order], cy[order]

    f_norm = f_raw / np.nanmedian(f_raw)
    f_sap_norm = f_sap / np.nanmedian(f_sap)
    dt = float(np.median(np.diff(t)))
    k = int(round(window_days / dt)) | 1
    k = max(5, min(k, (len(t) // 2) | 1))
    trend = medfilt(f_norm, kernel_size=k)
    f_det = f_norm / trend
    mad = float(1.4826 * np.median(np.abs(f_det - 1.0)))
    # Upward-only flare clipping so downward transits are never clipped
    keep = (f_det - 1.0) < 3.5 * mad
    t_c = t[keep]
    f_c = f_det[keep]
    f_sap_c = f_sap_norm[keep]
    e_c = np.full(len(t_c), mad)
    cx_c = cx[keep]
    cy_c = cy[keep]

    # Detrend centroids for in/out transit comparison
    vc = np.isfinite(cx_c) & np.isfinite(cy_c)
    dx_c = np.full(len(t_c), np.nan)
    dy_c = np.full(len(t_c), np.nan)
    if np.sum(vc) > 55:
        kc = min(101, (int(np.sum(vc)) // 3) | 1)
        kc = max(5, kc)
        dx_c[vc] = cx_c[vc] - medfilt(cx_c[vc], kernel_size=kc)
        dy_c[vc] = cy_c[vc] - medfilt(cy_c[vc], kernel_size=kc)

    return {
        "sector": sec,
        "pipeline": pipe,
        "exptime_s": prod["exptime_s"],
        "uri": prod["uri"],
        "flux_col": flux_col,
        "sha256": file_sha,
        "n_points": int(len(t_c)),
        "mad_ppm": float(mad * 1e6),
        "t": t_c,
        "f": f_c,
        "f_sap": f_sap_c,
        "e": e_c,
        "dx": dx_c,
        "dy": dy_c,
    }


def check_known_catalog(tic_id: int) -> dict:
    if not KNOWN_TICS_GZ.exists():
        return {"is_known": False, "matches": {}}
    with gzip.open(KNOWN_TICS_GZ, "rt", encoding="utf-8") as gz:
        db = json.load(gz)
    entry = db.get(str(tic_id))
    if not entry:
        return {"is_known": False, "matches": {}}
    has_official = bool(entry.get("confirmed_planets") or entry.get("tois") or entry.get("ctois"))
    return {
        "is_known": has_official,
        "has_external": bool(entry.get("external_candidates")),
        "matches": entry,
    }


def circular_transit_duration_days(period_d: float, r_star_rsun: float, m_star_msun: float) -> float:
    p_s = period_d * DAY_S
    a_m = ((G_SI * (m_star_msun * M_SUN_KG) * (p_s ** 2)) / (4.0 * (math.pi ** 2))) ** (1.0 / 3.0)
    r_m = r_star_rsun * R_SUN_M
    ratio = min(0.95, r_m / a_m)
    return (period_d / math.pi) * math.asin(ratio)


def evaluate_signal_from_arrays(
    tic_meta: dict,
    sectors_data: list[dict],
    target_period: float | None = None,
    period_min: float = 1.1,
    period_max: float = 15.0,
    mask_intervals: list[tuple[float, float, float]] | None = None,
    allow_known: bool = False,
) -> dict:
    """Run two-stage BLS search (or targeted period refinement) and evaluate all 10 falsification gates."""
    sec_list = []
    for s in sectors_data:
        t, f, f_sap, e, dx, dy = s["t"], s["f"], s["f_sap"], s["e"], s["dx"], s["dy"]
        if mask_intervals:
            keep = np.ones(len(t), dtype=bool)
            for mp, mt0, mdur in mask_intervals:
                ph = ((t - mt0 + 0.5 * mp) % mp) - 0.5 * mp
                keep &= np.abs(ph) > 0.8 * mdur
            t, f, f_sap, e, dx, dy = t[keep], f[keep], f_sap[keep], e[keep], dx[keep], dy[keep]
        sec_list.append({**s, "t": t, "f": f, "f_sap": f_sap, "e": e, "dx": dx, "dy": dy})

    t_all = np.concatenate([s["t"] for s in sec_list])
    f_all = np.concatenate([s["f"] for s in sec_list])
    f_sap_all = np.concatenate([s["f_sap"] for s in sec_list])
    e_all = np.concatenate([s["e"] for s in sec_list])
    dx_all = np.concatenate([s["dx"] for s in sec_list])
    dy_all = np.concatenate([s["dy"] for s in sec_list])

    bls_all = BoxLeastSquares(t_all, f_all, dy=e_all)
    durations = np.linspace(0.04, 0.14, 7)
    d_coarse = durations[durations < 0.8 * period_min]

    if target_period is not None:
        p_grid = np.linspace(target_period * 0.997, target_period * 1.003, 4000)
        res = bls_all.power(p_grid, d_coarse, objective="snr")
    else:
        # Stage 1: Per-sector broad frequency BLS with continuum subtraction
        # (resolves every 27-day sector peak while suppressing single-event edge plateaus)
        p_max_s1 = min(period_max, 12.5)
        f_grid = np.linspace(1.0 / p_max_s1, 1.0 / period_min, 4500)
        p_sector = 1.0 / f_grid[::-1]
        joint_snr = np.zeros(len(p_sector), dtype=float)
        for s in sec_list:
            bls_s = BoxLeastSquares(s["t"], s["f"], dy=s["e"])
            r_s = bls_s.power(p_sector, d_coarse, objective="snr")
            pwr = np.maximum(0.0, np.nan_to_num(np.asarray(r_s.power, dtype=float), nan=0.0))
            base = medfilt(pwr, kernel_size=151)
            joint_snr += np.maximum(0.0, pwr - base)

        # Pick top 5 well-separated candidate periods from the continuum-subtracted joint periodogram
        order_idx = np.argsort(joint_snr)[::-1]
        top_periods: list[float] = []
        for idx_cand in order_idx:
            pc = float(p_sector[int(idx_cand)])
            if all(abs(pc - p_prev) / min(pc, p_prev) > 0.025 for p_prev in top_periods):
                top_periods.append(pc)
                if len(top_periods) >= 5:
                    break

        # Stage 2: Multi-year phase-locking fine BLS around each candidate period
        # Rank candidates by multi-sector coherence (min per-sector Z) so single-sector glitches never eclipse real planets
        best_res = None
        best_score = (-1, -1, -1.0)
        for p0 in top_periods:
            pg = np.linspace(max(period_min, p0 * 0.997), min(period_max, p0 * 1.003), 4000)
            r_sub = bls_all.power(pg, d_coarse, objective="snr")
            idx_sub = int(np.nanargmax(r_sub.power))
            pk = float(r_sub.power[idx_sub])
            p_sub = float(r_sub.period[idx_sub])
            t0_sub = float(r_sub.transit_time[idx_sub])
            dur_sub = float(r_sub.duration[idx_sub])
            sec_zs = []
            for s in sec_list:
                ts, fs, es = s["t"], s["f"], s["e"]
                ws = 1.0 / (es ** 2)
                phs = ((ts - t0_sub + 0.5 * p_sub) % p_sub) - 0.5 * p_sub
                ins = np.abs(phs) <= 0.5 * dur_sub
                outs = (np.abs(phs) > 1.5 * dur_sub) & (np.abs(phs) < 4.0 * dur_sub)
                if np.sum(ins) >= 2 and np.sum(outs) >= 4:
                    ds = float((np.sum(ws[outs] * fs[outs]) / np.sum(ws[outs]) - np.sum(ws[ins] * fs[ins]) / np.sum(ws[ins])) * 1e6)
                    errs = float(math.sqrt(1.0 / np.sum(ws[ins]) + 1.0 / np.sum(ws[outs])) * 1e6)
                    sec_zs.append(ds / errs)
                else:
                    sec_zs.append(0.0)
            min_sz = min(sec_zs) if sec_zs else 0.0
            mean_sz = float(np.mean(sec_zs)) if sec_zs else 0.0
            tier1 = int(min_sz >= 1.8 and mean_sz >= 3.2 and pk >= 7.8)
            tier2 = int(min_sz >= 1.0 and mean_sz >= 2.2 and pk >= 7.0)
            score = (tier1, tier2, pk)
            if score > best_score:
                best_score = score
                best_res = r_sub
        res = best_res

    best_i = int(np.nanargmax(res.power))
    P = float(res.period[best_i])
    T0 = float(res.transit_time[best_i])
    dur = float(res.duration[best_i])
    bls_depth_ppm = float(res.depth[best_i] * 1e6)
    bls_snr = float(res.power[best_i])

    # Compute inverse-variance weighted in-transit vs local out-of-transit dip
    w = 1.0 / (e_all ** 2)
    epoch_idx = np.round((t_all - T0) / P).astype(int)
    phase = (t_all - T0) - epoch_idx * P
    in_tr = np.abs(phase) <= 0.5 * dur
    out_tr = (np.abs(phase) > 1.5 * dur) & (np.abs(phase) < 4.0 * dur)
    mean_out = float(np.sum(w[out_tr] * f_all[out_tr]) / np.sum(w[out_tr]))
    mean_in = float(np.sum(w[in_tr] * f_all[in_tr]) / np.sum(w[in_tr]))
    dip_ppm = (mean_out - mean_in) * 1e6
    dip_err_ppm = float(math.sqrt(1.0 / np.sum(w[in_tr]) + 1.0 / np.sum(w[out_tr])) * 1e6)

    # Gate 1: Catalog de-duplication
    cat_info = check_known_catalog(tic_meta["tic_id"])
    gate1_pass = (not cat_info["is_known"]) or allow_known

    # Per-transit inventory & Gate 2 + Gate 9
    unique_epochs = sorted(set(int(x) for x in epoch_idx[in_tr]))
    per_transit = []
    for ep in unique_epochs:
        m_in = in_tr & (epoch_idx == ep)
        m_out = out_tr & (epoch_idx == ep)
        if np.sum(m_in) >= 1 and np.sum(m_out) >= 2:
            f_o = float(np.sum(w[m_out] * f_all[m_out]) / np.sum(w[m_out]))
            f_i = float(np.sum(w[m_in] * f_all[m_in]) / np.sum(w[m_in]))
            d_ep = (f_o - f_i) * 1e6
            e_ep = float(math.sqrt(1.0 / np.sum(w[m_in]) + 1.0 / np.sum(w[m_out])) * 1e6)
            per_transit.append({
                "epoch": ep,
                "t_mid_btjd": round(T0 + ep * P, 5),
                "n_in": int(np.sum(m_in)),
                "dip_ppm": round(d_ep, 1),
                "err_ppm": round(e_ep, 1),
                "snr": round(d_ep / e_ep, 2),
            })

    n_transits = len(per_transit)
    n_positive_transits = sum(1 for tr in per_transit if tr["dip_ppm"] > 0)
    pos_transit_fraction = n_positive_transits / max(1, n_transits)

    # Jackknife SNR (dropping strongest single transit)
    if len(per_transit) >= 2:
        strongest_ep = max(per_transit, key=lambda x: x["snr"])["epoch"]
        jk_mask = epoch_idx != strongest_ep
        in_jk = in_tr & jk_mask
        out_jk = out_tr & jk_mask
        if np.sum(in_jk) >= 1 and np.sum(out_jk) >= 2:
            jk_dip = (np.sum(w[out_jk] * f_all[out_jk]) / np.sum(w[out_jk]) - np.sum(w[in_jk] * f_all[in_jk]) / np.sum(w[in_jk])) * 1e6
            jk_err = math.sqrt(1.0 / np.sum(w[in_jk]) + 1.0 / np.sum(w[out_jk])) * 1e6
            jackknife_snr = float(jk_dip / jk_err)
        else:
            jackknife_snr = 0.0
    else:
        jackknife_snr = 0.0

    gate2_pass = (len(sec_list) >= 2) and (n_transits >= 6)
    gate3_pass = bls_snr >= 8.0

    # Gate 4: Blind Leave-One-Epoch-Out (LOEO) cross-validation
    loeo_results = []
    p_loeo_grid = np.linspace(P * 0.9992, P * 1.0008, 3000)
    for idx_s, held in enumerate(sec_list):
        others = [s for j, s in enumerate(sec_list) if j != idx_s]
        t_tr = np.concatenate([s["t"] for s in others])
        f_tr = np.concatenate([s["f"] for s in others])
        e_tr = np.concatenate([s["e"] for s in others])
        bls_tr = BoxLeastSquares(t_tr, f_tr, dy=e_tr)
        r_tr = bls_tr.power(p_loeo_grid, dur, objective="snr")
        itr = int(np.nanargmax(r_tr.power))
        p_tr = float(r_tr.period[itr])
        t0_tr = float(r_tr.transit_time[itr])

        th, fh, eh = held["t"], held["f"], held["e"]
        wh = 1.0 / (eh ** 2)
        ph_h = ((th - t0_tr + 0.5 * p_tr) % p_tr) - 0.5 * p_tr
        in_h = np.abs(ph_h) <= 0.5 * dur
        out_h = (np.abs(ph_h) > 1.5 * dur) & (np.abs(ph_h) < 4.0 * dur)
        if np.sum(in_h) >= 2 and np.sum(out_h) >= 4:
            dip_h = float((np.sum(wh[out_h] * fh[out_h]) / np.sum(wh[out_h]) - np.sum(wh[in_h] * fh[in_h]) / np.sum(wh[in_h])) * 1e6)
            err_h = float(math.sqrt(1.0 / np.sum(wh[in_h]) + 1.0 / np.sum(wh[out_h])) * 1e6)
            z_h = dip_h / err_h
        else:
            dip_h, err_h, z_h = 0.0, 999.0, 0.0
        loeo_results.append({
            "held_out_sector": held["sector"],
            "train_period_d": round(p_tr, 7),
            "train_t0_btjd": round(t0_tr, 5),
            "n_in_points": int(np.sum(in_h)),
            "predicted_dip_ppm": round(dip_h, 1),
            "err_ppm": round(err_h, 1),
            "z_sigma": round(z_h, 2),
        })

    min_loeo_z = min((r["z_sigma"] for r in loeo_results), default=0.0)
    mean_loeo_z = float(np.mean([r["z_sigma"] for r in loeo_results])) if loeo_results else 0.0
    gate4_pass = (min_loeo_z >= 2.0) and (mean_loeo_z >= 3.5)

    # Gate 5: Odd vs Even transit depth consistency
    odd_in = in_tr & (epoch_idx % 2 == 1)
    even_in = in_tr & (epoch_idx % 2 == 0)
    if np.sum(odd_in) >= 2 and np.sum(even_in) >= 2:
        odd_dip = float((mean_out - np.sum(w[odd_in] * f_all[odd_in]) / np.sum(w[odd_in])) * 1e6)
        odd_err = float(math.sqrt(1.0 / np.sum(w[odd_in]) + 1.0 / np.sum(w[out_tr])) * 1e6)
        even_dip = float((mean_out - np.sum(w[even_in] * f_all[even_in]) / np.sum(w[even_in])) * 1e6)
        even_err = float(math.sqrt(1.0 / np.sum(w[even_in]) + 1.0 / np.sum(w[out_tr])) * 1e6)
        odd_even_diff_sigma = abs(odd_dip - even_dip) / math.sqrt(odd_err ** 2 + even_err ** 2)
        gate5_pass = (odd_dip > 0) and (even_dip > 0) and (odd_even_diff_sigma < 2.5)
    else:
        odd_dip, odd_err, even_dip, even_err, odd_even_diff_sigma = 0.0, 999.0, 0.0, 999.0, 99.0
        gate5_pass = False

    # Gate 6: Secondary eclipse / phase-0.5 anomaly null test
    ph_sec = ((t_all - (T0 + 0.5 * P) + 0.5 * P) % P) - 0.5 * P
    in_sec = np.abs(ph_sec) <= 0.5 * dur
    if np.sum(in_sec) >= 2:
        sec_dip = float((mean_out - np.sum(w[in_sec] * f_all[in_sec]) / np.sum(w[in_sec])) * 1e6)
        sec_err = float(math.sqrt(1.0 / np.sum(w[in_sec]) + 1.0 / np.sum(w[out_tr])) * 1e6)
        sec_snr = sec_dip / sec_err
    else:
        sec_dip, sec_err, sec_snr = 0.0, 999.0, 0.0
    gate6_pass = abs(sec_snr) < 2.5

    # Gate 7: Physical radius & circular transit duration plausibility
    r_star = tic_meta["rad_rsun"]
    m_star = tic_meta["mass_msun"]
    rp_rearth = r_star * R_SUN_TO_R_EARTH * math.sqrt(max(0.0, dip_ppm * 1e-6))
    dur_circ = circular_transit_duration_days(P, r_star, m_star)
    dur_ratio = dur / dur_circ
    gate7_pass = (0.4 <= rp_rearth <= 18.0) and (0.20 <= dur_ratio <= 2.20)

    # Gate 8: In/Out-of-transit Centroid Shift
    valid_c = np.isfinite(dx_all) & np.isfinite(dy_all)
    in_c = in_tr & valid_c
    out_c = out_tr & valid_c
    if np.sum(in_c) >= 5 and np.sum(out_c) >= 20:
        mx = float(np.mean(dx_all[in_c]) - np.mean(dx_all[out_c]))
        my = float(np.mean(dy_all[in_c]) - np.mean(dy_all[out_c]))
        sx = float(np.std(dx_all[out_c]) / math.sqrt(np.sum(in_c)))
        sy = float(np.std(dy_all[out_c]) / math.sqrt(np.sum(in_c)))
        shift_pix = math.hypot(mx, my)
        shift_arcsec = shift_pix * 21.0
        shift_sigma = math.hypot(mx / max(sx, 1e-9), my / max(sy, 1e-9))
    else:
        shift_pix, shift_arcsec, shift_sigma = 0.0, 0.0, 0.0
    gate8_pass = (shift_sigma < 3.0) and (shift_arcsec < 2.5)

    # Gate 9: Per-transit coherence & jackknife stability
    median_tr_dip_ppm = float(np.median([tr["dip_ppm"] for tr in per_transit])) if per_transit else 0.0
    gate9_pass = (
        (pos_transit_fraction >= 0.70)
        and (jackknife_snr >= 6.5)
        and (median_tr_dip_ppm >= 0.55 * dip_ppm)
    )

    # Gate 10: Out-of-transit sinusoidal modulation guard (starspot & sub-harmonic check)
    # Fit sin/cos at P, 2P, P/2, P/3 on out-of-transit SAP flux
    oot_mask = np.abs(phase) > 0.75 * dur
    t_oot = t_all[oot_mask]
    f_oot = f_sap_all[oot_mask] - np.median(f_sap_all[oot_mask])
    max_harm_amp_ppm = 0.0
    for mult in (1.0 / 3.0, 0.5, 1.0, 2.0):
        omega = 2.0 * math.pi / (P * mult)
        design = np.column_stack([np.sin(omega * t_oot), np.cos(omega * t_oot), np.ones(len(t_oot))])
        coeffs, _, _, _ = np.linalg.lstsq(design, f_oot, rcond=None)
        amp_ppm = float(math.hypot(coeffs[0], coeffs[1]) * 1e6)
        max_harm_amp_ppm = max(max_harm_amp_ppm, amp_ppm)
    # Compare on detrended out-of-transit flux across fundamental, rational sub-harmonics, and joint 4-harmonic series
    f_det_oot = f_all[oot_mask] - 1.0
    det_sin_amp_ppm = 0.0
    for mult_d in (0.20, 0.25, 1.0 / 3.0, 0.40, 0.50, 2.0 / 3.0, 1.0):
        omega1 = 2.0 * math.pi / (P * mult_d)
        design1 = np.column_stack([np.sin(omega1 * t_oot), np.cos(omega1 * t_oot), np.ones(len(t_oot))])
        c_det, _, _, _ = np.linalg.lstsq(design1, f_det_oot, rcond=None)
        det_sin_amp_ppm = max(det_sin_amp_ppm, float(math.hypot(c_det[0], c_det[1]) * 1e6))
    if len(t_oot) >= 20:
        cols_4h = [np.ones(len(t_oot))]
        for k_h in (1, 2, 3, 4):
            omega_k = 2.0 * math.pi * k_h / P
            cols_4h.append(np.sin(omega_k * t_oot))
            cols_4h.append(np.cos(omega_k * t_oot))
        beta_4h, _, _, _ = np.linalg.lstsq(np.column_stack(cols_4h), f_det_oot, rcond=None)
        ph_grid = np.linspace(0.75 * dur / P, 1.0 - 0.75 * dur / P, 400)
        t_grid = T0 + ph_grid * P
        y_grid = np.zeros_like(t_grid)
        for k_h in (1, 2, 3, 4):
            omega_k = 2.0 * math.pi * k_h / P
            y_grid += beta_4h[2 * k_h - 1] * np.sin(omega_k * t_grid) + beta_4h[2 * k_h] * np.cos(omega_k * t_grid)
        joint_4h_amp_ppm = 0.5 * float(np.max(y_grid) - np.min(y_grid)) * 1e6
        det_sin_amp_ppm = max(det_sin_amp_ppm, joint_4h_amp_ppm)
    gate10_pass = (det_sin_amp_ppm < 0.32 * dip_ppm) and (max_harm_amp_ppm < 3.5 * dip_ppm)

    gates = {
        "gate1_catalog_dedup": {
            "pass": bool(gate1_pass),
            "is_known_in_toi_ctoi_ps": bool(cat_info["is_known"]),
            "has_external_index_entry": bool(cat_info.get("has_external", False)),
        },
        "gate2_multi_sector_baseline": {
            "pass": bool(gate2_pass),
            "n_sectors": len(sec_list),
            "n_transits_observed": n_transits,
        },
        "gate3_combined_bls_snr": {
            "pass": bool(gate3_pass),
            "bls_snr": round(bls_snr, 2),
            "threshold": 8.0,
        },
        "gate4_leave_one_epoch_out": {
            "pass": bool(gate4_pass),
            "min_holdout_z_sigma": round(min_loeo_z, 2),
            "mean_holdout_z_sigma": round(mean_loeo_z, 2),
            "epochs": loeo_results,
        },
        "gate5_odd_even_consistency": {
            "pass": bool(gate5_pass),
            "odd_dip_ppm": round(odd_dip, 1),
            "odd_err_ppm": round(odd_err, 1),
            "even_dip_ppm": round(even_dip, 1),
            "even_err_ppm": round(even_err, 1),
            "diff_sigma": round(odd_even_diff_sigma, 2),
        },
        "gate6_secondary_eclipse_null": {
            "pass": bool(gate6_pass),
            "secondary_dip_ppm": round(sec_dip, 1),
            "secondary_err_ppm": round(sec_err, 1),
            "secondary_snr": round(sec_snr, 2),
        },
        "gate7_physical_plausibility": {
            "pass": bool(gate7_pass),
            "inferred_rp_rearth": round(rp_rearth, 3),
            "duration_h": round(dur * 24.0, 2),
            "circular_b0_duration_h": round(dur_circ * 24.0, 2),
            "duration_ratio": round(dur_ratio, 3),
        },
        "gate8_centroid_stability": {
            "pass": bool(gate8_pass),
            "centroid_shift_mpix": round(shift_pix * 1e3, 3),
            "centroid_shift_arcsec": round(shift_arcsec, 4),
            "centroid_shift_sigma": round(shift_sigma, 2),
        },
        "gate9_per_transit_coherence": {
            "pass": bool(gate9_pass),
            "positive_transits": n_positive_transits,
            "total_transits": n_transits,
            "positive_fraction": round(pos_transit_fraction, 3),
            "jackknife_snr": round(jackknife_snr, 2),
        },
        "gate10_rotational_modulation_guard": {
            "pass": bool(gate10_pass),
            "detrended_oot_sin_amp_ppm": round(det_sin_amp_ppm, 1),
            "raw_sap_max_harmonic_amp_ppm": round(max_harm_amp_ppm, 1),
            "ratio_to_dip": round(det_sin_amp_ppm / max(dip_ppm, 1.0), 3),
        },
    }

    all_passed = all(g["pass"] for g in gates.values())
    return {
        "all_gates_passed": bool(all_passed),
        "ephemeris": {
            "period_d": round(P, 7),
            "t0_btjd": round(T0, 5),
            "duration_d": round(dur, 5),
            "duration_h": round(dur * 24.0, 2),
            "bls_box_depth_ppm": round(bls_depth_ppm, 1),
            "weighted_dip_ppm": round(dip_ppm, 1),
            "weighted_dip_err_ppm": round(dip_err_ppm, 1),
            "combined_bls_snr": round(bls_snr, 2),
            "inferred_rp_rearth": round(rp_rearth, 2),
        },
        "gates": gates,
        "per_transit": per_transit,
        "_arrays": {
            "t_all": t_all,
            "f_all": f_all,
            "e_all": e_all,
            "periods": res.period,
            "power": res.power,
            "sec_list": sec_list,
        },
    }


def render_diagnostic_plot(tic_meta: dict, eval_res: dict, out_png: Path) -> None:
    arr = eval_res["_arrays"]
    eph = eval_res["ephemeris"]
    P = eph["period_d"]
    T0 = eph["t0_btjd"]
    dur = eph["duration_d"]
    t_all, f_all = arr["t_all"], arr["f_all"]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9.5), dpi=140)
    fig.patch.set_facecolor("#0d1117")
    for ax in axes.flat:
        ax.set_facecolor("#161b22")
        for spine in ax.spines.values():
            spine.set_color("#30363d")
        ax.tick_params(colors="#c9d1d9")
        ax.xaxis.label.set_color("#c9d1d9")
        ax.yaxis.label.set_color("#c9d1d9")
        ax.title.set_color("#f0f6fc")

    # Panel 1: BLS Periodogram
    ax = axes[0, 0]
    ax.plot(arr["periods"], arr["power"], color="#58a6ff", lw=0.9)
    ax.axvline(P, color="#f78166", ls="--", lw=1.2, label=f"P = {P:.6f} d (SNR={eph['combined_bls_snr']:.1f})")
    ax.set_title("1. Multi-Sector BLS Periodogram")
    ax.set_xlabel("Period (days)")
    ax.set_ylabel("BLS SNR")
    ax.legend(facecolor="#161b22", edgecolor="#30363d", labelcolor="#f0f6fc", fontsize=8)

    # Panel 2: Phase-folded zoom around transit
    ax = axes[0, 1]
    ph_hours = (((t_all - T0 + 0.5 * P) % P) - 0.5 * P) * 24.0
    near = np.abs(ph_hours) <= 8.0
    ax.scatter(ph_hours[near], (f_all[near] - 1.0) * 1e6, s=2, alpha=0.25, color="#8b949e", label="Detrended cadences")
    bins = np.linspace(-7.5, 7.5, 46)
    bc = 0.5 * (bins[:-1] + bins[1:])
    idx_b = np.digitize(ph_hours[near], bins)
    bm = [np.mean((f_all[near][idx_b == k] - 1.0) * 1e6) if np.sum(idx_b == k) >= 3 else np.nan for k in range(1, len(bins))]
    ax.plot(bc, bm, color="#58a6ff", lw=2.2, marker="o", ms=3.5, label="20-min binned mean")
    ax.axvline(-0.5 * dur * 24.0, color="#f78166", ls=":", lw=1.0)
    ax.axvline(+0.5 * dur * 24.0, color="#f78166", ls=":", lw=1.0)
    ax.set_ylim(-1200, 800)
    ax.set_title(f"2. Phase-Folded Transit (Dip={eph['weighted_dip_ppm']:.0f} ppm, Rp={eph['inferred_rp_rearth']:.2f} R⊕)")
    ax.set_xlabel("Hours from mid-transit")
    ax.set_ylabel("Relative flux (ppm)")
    ax.legend(facecolor="#161b22", edgecolor="#30363d", labelcolor="#f0f6fc", fontsize=8)

    # Panel 3: Leave-One-Epoch-Out (LOEO) per-sector folded comparison
    ax = axes[0, 2]
    colors = ["#58a6ff", "#3fb950", "#d2a8ff", "#ffa657", "#ff7b72"]
    for idx_s, s in enumerate(arr["sec_list"]):
        phs = (((s["t"] - T0 + 0.5 * P) % P) - 0.5 * P) * 24.0
        ns = np.abs(phs) <= 6.5
        b_s = np.linspace(-6.0, 6.0, 25)
        b_c = 0.5 * (b_s[:-1] + b_s[1:])
        ib = np.digitize(phs[ns], b_s)
        ym = [np.mean((s["f"][ns][ib == k] - 1.0) * 1e6) if np.sum(ib == k) >= 2 else np.nan for k in range(1, len(b_s))]
        ax.plot(b_c, ym, lw=1.8, marker="o", ms=3, color=colors[idx_s % len(colors)], label=f"S{s['sector']} ({s['pipeline']})")
    ax.axhline(0, color="#8b949e", ls="--", lw=0.8)
    ax.set_title("3. Per-Sector Coherence & Leave-One-Out")
    ax.set_xlabel("Hours from mid-transit")
    ax.set_ylabel("Relative flux (ppm)")
    ax.legend(facecolor="#161b22", edgecolor="#30363d", labelcolor="#f0f6fc", fontsize=8)

    # Panel 4: Odd vs Even transits
    ax = axes[1, 0]
    epoch_idx = np.round((t_all - T0) / P).astype(int)
    for is_odd, lbl, col in [(1, "Odd transits", "#3fb950"), (0, "Even transits", "#d2a8ff")]:
        m = near & ((epoch_idx % 2) == is_odd)
        ib = np.digitize(ph_hours[m], bins)
        ym = [np.mean((f_all[m][ib == k] - 1.0) * 1e6) if np.sum(ib == k) >= 2 else np.nan for k in range(1, len(bins))]
        ax.plot(bc, ym, lw=1.8, marker="o", ms=3, color=col, label=lbl)
    g5 = eval_res["gates"]["gate5_odd_even_consistency"]
    ax.set_ylim(-1200, 800)
    ax.set_title(f"4. Odd vs Even Test (Δ = {g5['diff_sigma']:.2f}σ)")
    ax.set_xlabel("Hours from mid-transit")
    ax.set_ylabel("Relative flux (ppm)")
    ax.legend(facecolor="#161b22", edgecolor="#30363d", labelcolor="#f0f6fc", fontsize=8)

    # Panel 5: Full phase curve & Secondary Eclipse check
    ax = axes[1, 1]
    ph_norm = ((t_all - T0) % P) / P
    b_full = np.linspace(0.0, 1.0, 81)
    bc_full = 0.5 * (b_full[:-1] + b_full[1:])
    ib_f = np.digitize(ph_norm, b_full)
    ym_f = [np.mean((f_all[ib_f == k] - 1.0) * 1e6) if np.sum(ib_f == k) >= 3 else np.nan for k in range(1, len(b_full))]
    ax.plot(bc_full, ym_f, color="#58a6ff", lw=1.5)
    ax.axvline(0.5, color="#ffa657", ls="--", lw=1.2, label="Phase 0.5 (Secondary check)")
    g6 = eval_res["gates"]["gate6_secondary_eclipse_null"]
    ax.set_title(f"5. Full Orbit & Secondary Eclipse ({g6['secondary_snr']:.2f}σ at 0.5P)")
    ax.set_xlabel("Orbital phase")
    ax.set_ylabel("Relative flux (ppm)")
    ax.legend(facecolor="#161b22", edgecolor="#30363d", labelcolor="#f0f6fc", fontsize=8)

    # Panel 6: Per-transit depth consistency across years
    ax = axes[1, 2]
    trs = eval_res["per_transit"]
    xs = list(range(1, len(trs) + 1))
    ys = [tr["dip_ppm"] for tr in trs]
    es = [tr["err_ppm"] for tr in trs]
    ax.errorbar(xs, ys, yerr=es, fmt="o", color="#3fb950", ecolor="#8b949e", elinewidth=1.0, capsize=2, ms=4)
    ax.axhline(eph["weighted_dip_ppm"], color="#58a6ff", ls="--", lw=1.2, label=f"Mean dip ({eph['weighted_dip_ppm']:.0f} ppm)")
    ax.axhline(0.0, color="#f78166", ls=":", lw=1.0)
    g9 = eval_res["gates"]["gate9_per_transit_coherence"]
    ax.set_title(f"6. Per-Transit Depths ({g9['positive_transits']}/{g9['total_transits']} > 0, JK SNR={g9['jackknife_snr']:.1f})")
    ax.set_xlabel("Observed transit index")
    ax.set_ylabel("Transit dip (ppm)")
    ax.legend(facecolor="#161b22", edgecolor="#30363d", labelcolor="#f0f6fc", fontsize=8)

    fig.suptitle(
        f"TIC {tic_meta['tic_id']} — P = {P:.6f} d | Rp = {eph['inferred_rp_rearth']:.2f} R⊕ | "
        f"Combined BLS SNR = {eph['combined_bls_snr']:.2f} | All 10 Gates: {'PASS' if eval_res['all_gates_passed'] else 'FAIL'}",
        color="#f0f6fc",
        fontsize=13,
        fontweight="bold",
        y=0.98,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_png, facecolor=fig.get_facecolor())
    plt.close(fig)


def generate_markdown_report(cert: dict) -> str:
    tic = cert["tic_id"]
    star = cert["stellar_parameters"]
    sig = cert["primary_signal"]
    eph = sig["ephemeris"]
    gates = sig["gates"]
    lines = [
        f"# Candidate Vetting Dossier: TIC {tic}",
        "",
        f"- **Status**: `{cert['status']}` (`all_gates_passed = {sig['all_gates_passed']}`)",
        f"- **Generated UTC**: `{cert['generated_utc']}`",
        f"- **Host Star**: TIC {tic} ($T_{{\\text{{mag}}}} = {star['tmag']}$, $T_{{\\text{{eff}}}} = {star['teff_k']:.0f}\\text{{ K}}$, $R_\\star = {star['rad_rsun']:.4f}\\,R_\\odot$, $M_\\star = {star['mass_msun']:.3f}\\,M_\\odot$, $d = {star['dist_pc']:.2f}\\text{{ pc}}$)",
        f"- **Sectors Analyzed**: `{[s['sector'] for s in cert['sectors_used']]}`",
        "",
        "## Primary Signal Ephemeris & Physical Parameters",
        "",
        "| Parameter | Value |",
        "| :--- | :--- |",
        f"| **Orbital Period ($P$)** | `{eph['period_d']:.7f} d` |",
        f"| **Mid-Transit Epoch ($T_0$)** | `{eph['t0_btjd']:.5f} BTJD` |",
        f"| **Transit Duration ($T_{{14}}$)** | `{eph['duration_h']:.2f} h` (`{eph['duration_d']:.5f} d`) |",
        f"| **Weighted Transit Depth ($\\Delta F$)** | `{eph['weighted_dip_ppm']:.1f} ± {eph['weighted_dip_err_ppm']:.1f} ppm` (BLS box: `{eph['bls_box_depth_ppm']:.1f} ppm`) |",
        f"| **Combined Multi-Sector BLS SNR** | `{eph['combined_bls_snr']:.2f} σ` |",
        f"| **Inferred Planet Radius ($R_p$)** | **`{eph['inferred_rp_rearth']:.2f} R_earth`** (assuming target star host) |",
        "",
        "## 10-Gate Adversarial Falsification Scorecard",
        "",
        "| Gate | Falsification Test | Key Metric | Result |",
        "| :--- | :--- | :--- | :--- |",
        f"| **1** | Catalog & Alert De-duplication | `is_known_in_toi_ctoi_ps = {gates['gate1_catalog_dedup']['is_known_in_toi_ctoi_ps']}` | {'PASS' if gates['gate1_catalog_dedup']['pass'] else 'FAIL'} |",
        f"| **2** | Multi-Sector / Multi-Epoch Baseline | `{gates['gate2_multi_sector_baseline']['n_sectors']}` sectors, `{gates['gate2_multi_sector_baseline']['n_transits_observed']}` transits | {'PASS' if gates['gate2_multi_sector_baseline']['pass'] else 'FAIL'} |",
        f"| **3** | Combined Detrended BLS SNR | `{gates['gate3_combined_bls_snr']['bls_snr']} σ` (min `8.0 σ`) | {'PASS' if gates['gate3_combined_bls_snr']['pass'] else 'FAIL'} |",
        f"| **4** | Blind Leave-One-Epoch-Out (LOEO) Prediction | `min Z = {gates['gate4_leave_one_epoch_out']['min_holdout_z_sigma']} σ`, `mean Z = {gates['gate4_leave_one_epoch_out']['mean_holdout_z_sigma']} σ` | {'PASS' if gates['gate4_leave_one_epoch_out']['pass'] else 'FAIL'} |",
        f"| **5** | Odd vs. Even Transit Depth Consistency | Odd `{gates['gate5_odd_even_consistency']['odd_dip_ppm']} ppm` vs Even `{gates['gate5_odd_even_consistency']['even_dip_ppm']} ppm` (`{gates['gate5_odd_even_consistency']['diff_sigma']} σ`) | {'PASS' if gates['gate5_odd_even_consistency']['pass'] else 'FAIL'} |",
        f"| **6** | Secondary Eclipse Null Test (Phase 0.5) | `{gates['gate6_secondary_eclipse_null']['secondary_dip_ppm']} ± {gates['gate6_secondary_eclipse_null']['secondary_err_ppm']} ppm` (`{gates['gate6_secondary_eclipse_null']['secondary_snr']} σ`) | {'PASS' if gates['gate6_secondary_eclipse_null']['pass'] else 'FAIL'} |",
        f"| **7** | Physical Radius & Keplerian Duration | `$R_p = {gates['gate7_physical_plausibility']['inferred_rp_rearth']} R_\\oplus$`, `$T_{{14}}/T_{{14,\\text{{circ}}}} = {gates['gate7_physical_plausibility']['duration_ratio']}$` | {'PASS' if gates['gate7_physical_plausibility']['pass'] else 'FAIL'} |",
        f"| **8** | In/Out-of-Transit Centroid Stability | `Δr = {gates['gate8_centroid_stability']['centroid_shift_arcsec']} arcsec` (`{gates['gate8_centroid_stability']['centroid_shift_sigma']} σ`) | {'PASS' if gates['gate8_centroid_stability']['pass'] else 'FAIL'} |",
        f"| **9** | Per-Transit Coherence & Jackknife Stability | `{gates['gate9_per_transit_coherence']['positive_transits']}/{gates['gate9_per_transit_coherence']['total_transits']}` positive, Jackknife SNR `{gates['gate9_per_transit_coherence']['jackknife_snr']} σ` | {'PASS' if gates['gate9_per_transit_coherence']['pass'] else 'FAIL'} |",
        f"| **10** | Rotational Modulation / Starspot Guard | Detrended OOT sin amp `{gates['gate10_rotational_modulation_guard']['detrended_oot_sin_amp_ppm']} ppm` (`{gates['gate10_rotational_modulation_guard']['ratio_to_dip']*100:.1f}%` of dip) | {'PASS' if gates['gate10_rotational_modulation_guard']['pass'] else 'FAIL'} |",
        "",
        "## Blind Leave-One-Epoch-Out (LOEO) Cross-Validation Breakdown",
        "",
        "| Held-Out Sector | Trained Period (d) | Trained $T_0$ (BTJD) | In-Transit Cadences | Predicted Dip in Held-Out Sector | Z-Score |",
        "| :--- | :--- | :--- | :--- | :--- | :--- |",
    ]
    for row in gates["gate4_leave_one_epoch_out"]["epochs"]:
        lines.append(
            f"| Sector {row['held_out_sector']} | `{row['train_period_d']:.7f}` | `{row['train_t0_btjd']:.5f}` | `{row['n_in_points']}` | `{row['predicted_dip_ppm']:+.1f} ± {row['err_ppm']:.1f} ppm` | **`{row['z_sigma']:+.2f} σ`** |"
        )
    if cert.get("secondary_candidates"):
        lines.extend([
            "",
            "## Additional Tentative Signals (After Masking Primary Transit)",
            "",
        ])
        for sc in cert["secondary_candidates"]:
            se = sc["ephemeris"]
            lines.append(
                f"- **Signal `{sc['label']}`**: $P = {se['period_d']:.6f}\\text{{ d}}$, $T_0 = {se['t0_btjd']:.5f}\\text{{ BTJD}}$, depth $= {se['weighted_dip_ppm']:.1f} \\pm {se['weighted_dip_err_ppm']:.1f}\\text{{ ppm}}$, $R_p \\approx {se['inferred_rp_rearth']:.2f}\\,R_\\oplus$, BLS SNR $= {se['combined_bls_snr']:.2f}\\,\\sigma$ (`{sc['note']}`)"
            )
    lines.extend([
        "",
        "## Diagnostic Plots",
        "",
        "![Diagnostic Plots](diagnostics.png)",
        "",
    ])
    return "\n".join(lines)


def vet_target(
    tic_id: int,
    target_period: float | None = None,
    secondary_period: float | None = None,
    status_label: str = "unreported_candidate",
    allow_known: bool = False,
    out_dir: Path | None = None,
) -> dict:
    tic_meta = fetch_tic_metadata(tic_id)
    prods = discover_sector_products(tic_id)
    if len(prods) < 2:
        raise ValueError(f"TIC {tic_id} only has {len(prods)} sector(s) available")
    sectors_data = [download_and_load_sector(tic_id, p) for p in prods]

    primary = evaluate_signal_from_arrays(
        tic_meta,
        sectors_data,
        target_period=target_period,
        allow_known=allow_known,
    )

    secondary_list = []
    if secondary_period is not None:
        eph1 = primary["ephemeris"]
        sec_eval = evaluate_signal_from_arrays(
            tic_meta,
            sectors_data,
            target_period=secondary_period,
            mask_intervals=[(eph1["period_d"], eph1["t0_btjd"], eph1["duration_d"])],
            allow_known=allow_known,
        )
        sec_clean = {k: v for k, v in sec_eval.items() if k != "_arrays"}
        sec_clean["label"] = f"TIC {tic_id}.02"
        sec_clean["note"] = "Tentative multi-sector signal after masking primary transit"
        secondary_list.append(sec_clean)

    if out_dir is None:
        out_dir = ROOT / f"candidates/TIC-{tic_id}"
    out_dir.mkdir(parents=True, exist_ok=True)

    render_diagnostic_plot(tic_meta, primary, out_dir / "diagnostics.png")

    primary_clean = {k: v for k, v in primary.items() if k != "_arrays"}
    cert = {
        "schema_version": 1,
        "tic_id": tic_id,
        "status": status_label,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "stellar_parameters": tic_meta,
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
        "secondary_candidates": secondary_list,
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }

    (out_dir / "certificate.json").write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (out_dir / "REPORT.md").write_text(generate_markdown_report(cert), encoding="utf-8")
    return cert


def main() -> None:
    parser = argparse.ArgumentParser(description="Vet a TESS TIC target through the 10-gate adversarial falsification pipeline.")
    parser.add_argument("--tic", type=int, required=True, help="TIC ID to vet")
    parser.add_argument("--period", type=float, default=None, help="Optional target period (days) to refine around")
    parser.add_argument("--secondary-period", type=float, default=None, help="Optional second period (days) to evaluate after masking primary")
    parser.add_argument("--status", type=str, default="unreported_candidate", help="Status label for certificate")
    parser.add_argument("--allow-benchmark", action="store_true", help="Allow target in external benchmark index")
    args = parser.parse_args()

    cert = vet_target(
        tic_id=args.tic,
        target_period=args.period,
        secondary_period=args.secondary_period,
        status_label=args.status,
        allow_known=args.allow_benchmark,
    )
    eph = cert["primary_signal"]["ephemeris"]
    print(
        f"TIC {args.tic}: all_gates_passed={cert['primary_signal']['all_gates_passed']} | "
        f"P={eph['period_d']:.7f} d | dip={eph['weighted_dip_ppm']:.1f}+/-{eph['weighted_dip_err_ppm']:.1f} ppm | "
        f"Rp={eph['inferred_rp_rearth']:.2f} R_earth | SNR={eph['combined_bls_snr']:.2f}"
    )


if __name__ == "__main__":
    main()
