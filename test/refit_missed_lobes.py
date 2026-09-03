#!/usr/bin/env python3
"""Chantier option 1: targeted Gaussian refit of the expected missing lobe.

For each short-N chain, estimate the expected next-lobe position (last lobe
+ mean chain spacing along the chain direction), fit a 2D Gaussian in a
window around it (plane-flattened), and report the amplitude signal-to-noise
ratio. Label-free: reads the Z channel and fitted lobe positions only.

Usage:
  python3 test/refit_missed_lobes.py FEATURES_TSV NPY_DIR OUT_TSV
"""

import csv
import os
import sys

import numpy as np
from scipy.optimize import least_squares

PX = 10.0 / 512
WIN = 0.55  # half-window in nm


def read_z(path):
    return np.fromfile(path, dtype=np.float64).reshape(512, 512)


def lobe_positions(features_path, fname):
    out = []
    with open(features_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["file"] == fname:
                out.append((float(r["x_nm"]), float(r["y_nm"]), float(r["t_nm"])))
    return sorted(out, key=lambda p: p[2])


def chain_spacing(lobes):
    d = [np.hypot(lobes[i + 1][0] - lobes[i][0], lobes[i + 1][1] - lobes[i][1])
         for i in range(len(lobes) - 1)]
    return float(np.mean(d)) if d else 0.6


def window(img, cx, cy):
    x0 = max(0, int((cx - WIN) / PX))
    x1 = min(512, int((cx + WIN) / PX) + 1)
    y0 = max(0, int((cy - WIN) / PX))
    y1 = min(512, int((cy + WIN) / PX) + 1)
    if x1 - x0 < 10 or y1 - y0 < 10:
        return None, None, None, None
    sub = img[y0:y1, x0:x1].astype(float)
    h, w = sub.shape
    ys, xs = np.mgrid[0:h, 0:w]
    A = np.column_stack([xs.ravel(), ys.ravel(), np.ones(h * w)])
    coef, *_ = np.linalg.lstsq(A, sub.ravel(), rcond=None)
    flat = (sub - (coef[0] * xs + coef[1] * ys + coef[2])).ravel()
    return flat.reshape(h, w), x0, y0, (w, h)


def fit_gauss(img, cx, cy, angle, win=WIN):
    """fit a 2D Gaussian (amplitude, center, sigma_t, sigma_u) in the window;
    returns (A, cx, cy, st, su, snr, resid_std) or None if the window is out."""
    wdata = window(img, cx, cy)
    if wdata is None:
        return None
    sub, x0, y0, (w, h) = wdata
    ys, xs = np.mgrid[0:h, 0:w]
    th = np.radians(angle)
    ct, st = np.cos(th), np.sin(th)
    cx0, cy0 = (cx / PX - x0), (cy / PX - y0)

    def model(p):
        A, dx, dy, st_, su_ = p
        t = (xs - cx0 - dx) * ct + (ys - cy0 - dy) * st
        u = -(xs - cx0 - dx) * st + (ys - cy0 - dy) * ct
        return A * np.exp(-0.5 * (t / max(st_, 0.05))**2 - 0.5 * (u / max(su_, 0.05))**2)

    def resid(p):
        return (model(p) - sub).ravel()

    best = None
    for A0 in (np.percentile(sub, 99) - np.percentile(sub, 50),):
        for st0, su0 in ((0.14, 0.10), (0.10, 0.14), (0.12, 0.12)):
            p0 = (A0, 0.0, 0.0, st0, su0)
            try:
                r = least_squares(resid, p0, method="lm", max_nfev=200)
            except Exception:
                continue
            A, dx, dy, st_, su_ = r.x
            if st_ < 0.03 or su_ < 0.03 or abs(dx) > 0.5 * win / PX or abs(dy) > 0.5 * win / PX:
                continue
            resid_std = float(np.std(r.fun))
            snr = abs(A) / max(resid_std, 1e-12)
            if best is None or snr > best[0]:
                best = (snr, A, dx, dy, st_, su_, resid_std)
    return best


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    features_path, npy_dir, out_tsv = sys.argv[1:4]

    with open(out_tsv, "w") as g:
        g.write("file\tn_lobes\tspacing\texpected_dist\tA\tdx_nm\tdy_nm\tst_nm\tsu_nm\tsnr\tresid_std\n")
        for fname in sorted(os.listdir(npy_dir)):
            if not fname.endswith(".npy"):
                continue
            base = fname[:-4]
            img = read_z(os.path.join(npy_dir, fname))
            lobes = lobe_positions(features_path, base + ".sxm")
            n = len(lobes)
            if n < 3:
                continue
            xl, yl, _ = lobes[-1]
            xp, yp, _ = lobes[-2]
            xq, yq, _ = lobes[-3]
            d1 = np.array([xl - xp, yl - yp])
            d2 = np.array([xp - xq, yp - yq])
            dx, dy = d1 + 0.5 * d2
            L = np.hypot(dx, dy)
            if L < 1e-9:
                continue
            dx, dy = dx / L, dy / L
            spacing = chain_spacing(lobes)
            angle = np.degrees(np.arctan2(dy, dx))
            ex, ey = xl + dx * spacing, yl + dy * spacing
            fit = fit_gauss(img, ex, ey, angle)
            if fit is None:
                g.write(f"{base}\t{n}\t{spacing:.3f}\t{spacing:.3f}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
                print(f"{base}: fenetre hors image")
                continue
            snr, A, ddx, ddy, st_, su_, rstd = fit
            g.write(f"{base}\t{n}\t{spacing:.3f}\t{spacing:.3f}\t{A:.5f}\t{ddx*PX:.3f}\t"
                    f"{ddy*PX:.3f}\t{st_*PX:.3f}\t{su_*PX:.3f}\t{snr:.2f}\t{rstd:.6f}\n")
            print(f"{base}: A={A:.4f} snr={snr:.1f} dx={ddx*PX:.2f} st={st_*PX:.2f} su={su_*PX:.2f}")
    print(f"rapport: {out_tsv}")


if __name__ == "__main__":
    main()
