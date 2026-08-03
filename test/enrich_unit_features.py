#!/usr/bin/env python3
"""Enrich per-lobe unit-assignment features with physical patch statistics.

Label-free: reads only feature/patch TSVs produced by the fit pipeline; never
benchmark truth, sequences, or class counts. The extra columns describe the
2D shape of each lobe's backward/difference patch (raw and residual) and the
local chain geometry.

Patch conventions (see test/extract_lobe_patches_bwd.jl):
  - patches are 9x9, z-scored per patch (median/std); pixel order row-major
  - t: along the chain, u: transverse; coordinates -4..4 (pixel units)
  - bwd_raw / bwd_res : backward-channel raw / residual (Gaussian model)
  - diff_raw / diff_res : forward-backward difference raw / residual
    (topography cancels in diff; model cancels in diff_res)

New columns per family F in {bwd_res, diff_res, diff_raw, bwd_raw}:
  F_en        sum p^2                (residual/raw signal energy)
  F_com_u     sum p*u / sum|p|       (transverse centroid)
  F_com_t     sum p*t / sum|p|       (along-chain centroid)
  F_sku       sum p*u^3 / sum|p|     (transverse skewness, 3rd moment)
  F_skt       sum p*t^3 / sum|p|     (along-chain skewness)
  F_asym_u    (sum_{u>0} p - sum_{u<0} p) / sum|p|   (left/right parity)
  F_asym_t    (sum_{t>0} p - sum_{t<0} p) / sum|p|   (up/down parity)
and for the two residual families only:
  F_kurt_u    sum p*u^4 / sum|p|     (transverse kurtosis)
  F_ring1..3  annular means r=1..3 / sum|p| (radial localization)

Geometry columns (from the feature table itself):
  spacing_next_nm  distance to the next lobe along the chain
  spacing_asym     (prev - next) / (prev + next)
  chain_curv_deg   angle between the two chain segments at this lobe
  elongation       sigma_parallel_nm / sigma_perp_nm

Usage:
  python3 test/enrich_unit_features.py --patches PATCHES.tsv \
      --features FEATURES.tsv --out OUT.tsv
"""

import argparse
import csv
import math

import numpy as np


def load_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def patch_moments(vals, coords):
    """Moments of a normalized 9x9 patch. Returns dict or None if not finite."""
    try:
        p = np.array(vals, dtype=float)
    except ValueError:
        return None
    if not np.all(np.isfinite(p)):
        return None
    side = int(round(math.sqrt(len(p))))
    p = p.reshape(side, side)
    grid = np.array(coords[:side])
    # row-major: index i -> row (t), j -> col (u)
    T, U = np.meshgrid(grid, grid, indexing="ij")
    w = np.abs(p)
    s = w.sum()
    if s <= 1e-12:
        return None
    m = {
        "en": float((p * p).sum()),
        "com_u": float((p * U).sum() / s),
        "com_t": float((p * T).sum() / s),
        "sku": float((p * U**3).sum() / s),
        "skt": float((p * T**3).sum() / s),
        "asym_u": float((p * np.sign(U)).sum() / s),
        "asym_t": float((p * np.sign(T)).sum() / s),
    }
    m["kurt_u"] = float((p * U**4).sum() / s)
    r = np.sqrt(T**2 + U**2)
    for ring in (1, 2, 3):
        ring_w = (r > ring - 0.5) & (r <= ring + 0.5)
        m[f"ring{ring}"] = float((p * ring_w).sum() / s)
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--patches", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    patch_rows = load_tsv(args.patches)
    feat_rows = load_tsv(args.features)
    print(f"patches: {len(patch_rows)} rows, features: {len(feat_rows)} rows")

    coords = list(range(-4, 5))  # 9x9 patch, pixel units

    fam_cols = {}
    fam_sizes = {}
    for fam in ("bwd_res", "diff_res", "diff_raw", "bwd_raw"):
        cols = [c for c in patch_rows[0] if c.startswith(f"{fam}_p")]
        fam_sizes[fam] = len(cols)
        fam_cols[fam] = cols
        if len(cols) != 81:
            raise SystemExit(f"{fam}: expected 81 patch cols, got {len(cols)}")

    # precompute per-lobe patch moments
    patch_mom = {}
    for row in patch_rows:
        key = (row["file"], int(row["lobe"]))
        moms = {}
        ok = True
        for fam in ("bwd_res", "diff_res", "diff_raw", "bwd_raw"):
            m = patch_moments([row[c] for c in fam_cols[fam]], coords)
            if m is None:
                ok = False
                break
            moms[fam] = m
        if ok:
            patch_mom[key] = moms

    print(f"patches with finite moments: {len(patch_mom)}/{len(patch_rows)}")

    # per-file sorted positions for chain geometry
    by_file = {}
    for row in feat_rows:
        by_file.setdefault(row["file"], []).append(row)
    for rows in by_file.values():
        rows.sort(key=lambda r: float(r["t_nm"]))

    out_cols = list(feat_rows[0].keys())
    new_cols = []
    for fam in ("bwd_res", "diff_res", "diff_raw", "bwd_raw"):
        for k in ("en", "com_u", "com_t", "sku", "skt", "asym_u", "asym_t"):
            new_cols.append(f"{fam}_{k}")
        if fam in ("bwd_res", "diff_res"):
            for k in ("kurt_u", "ring1", "ring2", "ring3"):
                new_cols.append(f"{fam}_{k}")
    new_cols += ["spacing_next_nm", "spacing_asym", "chain_curv_deg", "elongation"]
    if "skew_ratio" in feat_rows[0]:
        new_cols.append("split_log_skew")
    out_cols += new_cols

    missing_patch = 0
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for row in feat_rows:
            key = (row["file"], int(row["lobe"]))
            moms = patch_mom.get(key)
            if moms is None:
                missing_patch += 1
                w.writerow(row)  # keep the row; new cols stay absent -> NA
                continue
            for fam in ("bwd_res", "diff_res", "diff_raw", "bwd_raw"):
                for k, v in moms[fam].items():
                    row[f"{fam}_{k}"] = f"{v:.8g}"

            # chain geometry
            rows = by_file[row["file"]]
            i = next(idx for idx, r in enumerate(rows) if r is row)
            prev_t = float(rows[i - 1]["t_nm"]) if i > 0 else math.nan
            prev_u = float(rows[i - 1]["u_nm"]) if i > 0 else math.nan
            next_t = float(rows[i + 1]["t_nm"]) if i + 1 < len(rows) else math.nan
            next_u = float(rows[i + 1]["u_nm"]) if i + 1 < len(rows) else math.nan
            t = float(row["t_nm"]); u = float(row["u_nm"])
            row["spacing_next_nm"] = (f"{math.hypot(next_t - t, next_u - u):.8g}"
                                      if math.isfinite(next_t) else "NA")
            sp = (math.hypot(t - prev_t, u - prev_u)
                  if math.isfinite(prev_t) else math.nan)
            sn = (math.hypot(next_t - t, next_u - u)
                  if math.isfinite(next_t) else math.nan)
            if math.isfinite(sp) and math.isfinite(sn) and sp + sn > 0:
                row["spacing_asym"] = f"{(sp - sn) / (sp + sn):.8g}"
            else:
                row["spacing_asym"] = "NA"
            if (math.isfinite(prev_t) and math.isfinite(next_t)):
                v1 = (t - prev_t, u - prev_u)
                v2 = (next_t - t, next_u - u)
                n1 = math.hypot(*v1); n2 = math.hypot(*v2)
                if n1 > 0 and n2 > 0:
                    cosang = max(-1.0, min(1.0, (v1[0] * v2[0] + v1[1] * v2[1]) / (n1 * n2)))
                    row["chain_curv_deg"] = f"{math.degrees(math.acos(cosang)):.8g}"
                else:
                    row["chain_curv_deg"] = "NA"
            else:
                row["chain_curv_deg"] = "NA"
            spar = float(row["sigma_parallel_nm"]); sperp = float(row["sigma_perp_nm"])
            if math.isfinite(spar) and math.isfinite(sperp) and sperp > 0:
                row["elongation"] = f"{spar / sperp:.8g}"
            else:
                row["elongation"] = "NA"
            if "split_log_skew" in out_cols:
                skew = row.get("skew_ratio", "NA")
                try:
                    skew_f = float(skew)
                except ValueError:
                    skew_f = math.nan
                if math.isfinite(skew_f) and skew_f > 0:
                    row["split_log_skew"] = f"{math.log(skew_f):.8g}"
                else:
                    row["split_log_skew"] = "NA"
            w.writerow(row)

    print(f"wrote {args.out}: {len(out_cols)} cols, {len(feat_rows) - missing_patch} rows "
          f"({missing_patch} rows skipped for missing/NaN patches)")


if __name__ == "__main__":
    main()
