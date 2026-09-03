#!/usr/bin/env python3
"""Chantier: mold-assisted detection of missed terminal lobes (step 1 of the
end-to-end pipeline). Slides the cc GlcN/GlcNAc mold along the chain beyond
the last fitted lobe (after plane flattening) and reports the correlation
peak. Label-free: reads only the Z channel and the fitted lobe positions.

Usage:
  python3 test/detect_missed_lobes.py FEATURES_TSV TEMPLATES_TSV NPY_DIR OUT_TSV
"""

import csv
import os
import sys

import numpy as np

SIDE = 17
COORDS = np.arange(-0.32, 0.321, 0.04)
PX = 10.0 / 512


def load_template(path, typ):
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if int(r["type"]) == typ and int(r["parity"]) == 0 and int(r["mirror"]) == 0:
                t = np.array([float(r[f"p{i:03d}"]) for i in range(1, 290)]).reshape(17, 17)
                return t - t.mean()
    raise SystemExit(f"template type {typ} not found")


def read_z(path):
    return np.fromfile(path, dtype=np.float64).reshape(512, 512)


def lobe_positions(features_path, fname):
    out = []
    with open(features_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["file"] == fname:
                out.append((float(r["x_nm"]), float(r["y_nm"]), float(r["t_nm"])))
    return sorted(out, key=lambda p: p[2])


def ncc_at(img, cx, cy, tmpl, dx, dy):
    nx, ny = -dy, dx
    patch = np.zeros((17, 17))
    ok = np.zeros((17, 17), dtype=bool)
    for i, tt in enumerate(COORDS):
        for j, uu in enumerate(COORDS):
            x = cx + dx * tt + nx * uu
            y = cy + dy * tt + ny * uu
            px = int(round(x / PX))
            py = int(round(y / PX))
            if 0 <= px < 512 and 0 <= py < 512:
                patch[i, j] = img[py, px]
                ok[i, j] = True
    m = ok
    if m.sum() < 100:
        return float("nan")
    rows, cols = np.where(m)
    vals = patch[m]
    A = np.column_stack([rows.astype(float), cols.astype(float), np.ones(len(rows))])
    coef, *_ = np.linalg.lstsq(A, vals, rcond=None)
    flat = patch.copy()
    flat[m] = vals - A @ coef
    a = flat[m] - flat[m].mean()
    b = tmpl[m]
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / den) if den > 0 else 0.0


def main():
    if len(sys.argv) < 5:
        print(__doc__)
        sys.exit(1)
    features_path, templates_path, npy_dir, out_tsv = sys.argv[1:5]
    tmpl0 = load_template(templates_path, 0)
    tmpl1 = load_template(templates_path, 1)

    with open(out_tsv, "w") as g:
        g.write("file\tn_lobes\tchain_dir\tlast_lobe\tdist_peak\tpeak_glcn\tpeak_glcnac\tmean_ref_glcn\n")
        for fname in sorted(os.listdir(npy_dir)):
            if not fname.endswith(".npy"):
                continue
            base = fname[:-4]
            img = read_z(os.path.join(npy_dir, fname))
            lobes = lobe_positions(features_path, base + ".sxm")
            n = len(lobes)
            if n < 3:
                continue
            x5, y5, _ = lobes[-1]
            x4, y4, _ = lobes[-2]
            x3, y3, _ = lobes[-3]
            d1 = np.array([x5 - x4, y5 - y4])
            d2 = np.array([x4 - x3, y4 - y3])
            dx, dy = d1 + 0.5 * d2
            L = np.hypot(dx, dy)
            if L < 1e-9:
                continue
            dx, dy = dx / L, dy / L
            best = (0.0, None, None)
            for dist in np.arange(0.25, 1.35, 0.02):
                cx, cy = x5 + dx * dist, y5 + dy * dist
                n0 = ncc_at(img, cx, cy, tmpl0, dx, dy)
                n1 = ncc_at(img, cx, cy, tmpl1, dx, dy)
                if n0 > best[0]:
                    best = (n0, dist, n1)
            ref = [ncc_at(img, x, y, tmpl0, dx, dy) for (x, y, _) in lobes]
            ref_mean = float(np.mean(ref)) if ref else 0.0
            if best[1] is None:
                g.write(f"{base}\t{n}\t({dx:.3f},{dy:.3f})\t({x5:.2f},{y5:.2f})\t"
                        f"NA\tNA\tNA\t{ref_mean:.4f}\n")
                print(f"{base}: n={n} scan hors image")
                continue
            g.write(f"{base}\t{n}\t({dx:.3f},{dy:.3f})\t({x5:.2f},{y5:.2f})\t"
                    f"{best[1]:.3f}\t{best[0]:.4f}\t{best[2]:.4f}\t{ref_mean:.4f}\n")
            print(f"{base}: n={n} pic@dist={best[1]:.2f} NCC_GlcN={best[0]:.3f} "
                  f"NCC_GlcNAc={best[2]:.3f} ref_GlcN={ref_mean:.3f}")
    print(f"rapport: {out_tsv}")


if __name__ == "__main__":
    main()
