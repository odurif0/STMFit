#!/usr/bin/env python3
"""Adaptive-contour (constant-current) DFT-STM mold builder.

Label-free: reads only QE LDOS cubes, relaxed-geometry frames, and fit
feature/patch TSVs produced by the STM pipeline. Never reads benchmark truth,
sequences, or composition counts.

Pipeline (journal 2026-08-01/02, sections 8h-8r):
  1. parse the QE Gaussian cube (production cubes are orthogonal 240x180x250
     bohr grids; the cube may be written with non-orthogonal axes in general),
  2. trilinearly sample the LDOS onto a regular grid aligned with the slab
     frame (origin/t/u axes from the frame TSV),
  3. find the isovalue whose first-vacuum isosurface mean height equals the
     target height (0.50 nm) with >= 30% valid columns,
  4. z-score the height map into a 17x17 template (NA -> 0 after z-score),
     plus the 4 parity/mirror variants.

The resulting template TSV is scored against experimental residual patches by
test/score_connected_mold_templates.jl; the per-lobe cost_margin is the mold
feature used by the label-free GMM/k-means pipeline (see
test/build_cc_soft_champion.py).

Usage:
  python3 test/lib/cc_mold_builder.py CUBE0 CUBE1 FRAME0 FRAME1 OUT_TEMPLATE \
      [--height 0.50] [--half-nm 0.32] [--step-nm 0.04]
"""

import numpy as np
import sys

BOHR = 0.05291772109  # nm per bohr
SIDE = 17
COORDS = np.arange(-0.32, 0.321, 0.04)
TT, UU = np.meshgrid(COORDS, COORDS, indexing="ij")
DISK = (TT**2 + UU**2) <= (0.32**2)


def read_cube(path):
    """Parse a Gaussian cube (bohr units) into (origin_nm, axes, values)."""
    lines = open(path).read().splitlines()
    nat = abs(int(lines[2].split()[0]))
    origin = np.array([float(x) for x in lines[2].split()[1:4]]) * BOHR
    axes = []
    for ax in range(3):
        p = lines[3 + ax].split()
        axes.append((abs(int(p[0])), np.array([float(x) for x in p[1:4]]) * BOHR))
    vals = []
    for line in lines[6 + nat:]:
        for tok in line.split():
            try:
                vals.append(float(tok.replace("D", "E")))
            except ValueError:
                pass
    (n1, v1), (n2, v2), (n3, v3) = axes
    return origin, axes, np.array(vals[: n1 * n2 * n3])


def sample_volume(origin, axes, vals, pts):
    """Trilinear sample of the cube at Nx3 points (nm)."""
    (n1, v1), (n2, v2), (n3, v3) = axes
    M = np.column_stack([v1, v2, v3])
    Minv = np.linalg.inv(M)
    c = (pts - origin) @ Minv.T
    i0 = np.floor(c).astype(int)
    f = c - i0
    out = np.zeros(len(pts))
    wsum = np.zeros(len(pts))
    for di in (0, 1):
        for dj in (0, 1):
            for dk in (0, 1):
                ii = i0[:, 0] + di
                jj = i0[:, 1] + dj
                kk = i0[:, 2] + dk
                ok = (ii >= 0) & (ii < n1) & (jj >= 0) & (jj < n2) & (kk >= 0) & (kk < n3)
                w = (f[:, 0] if di else 1 - f[:, 0]) * (f[:, 1] if dj else 1 - f[:, 1]) * (f[:, 2] if dk else 1 - f[:, 2])
                idx = ii + n1 * jj + n1 * n2 * kk
                out[ok] += w[ok] * vals[idx[ok]]
                wsum[ok] += w[ok]
    out[wsum == 0] = np.nan
    return out


def surface_grid(cube_path, frame, t_range, u_range, z_step=0.005):
    """Sample LDOS columns along the frame normal; return (zs, NxNz matrix)."""
    origin, axes, vals = read_cube(cube_path)
    a1, a2, a3 = axes
    o = np.array(frame[0])
    th = np.array(frame[1])
    uh = np.array(frame[2])
    nh = np.cross(th, uh)
    nh /= np.linalg.norm(nh)
    tt, uu = np.meshgrid(t_range, u_range)
    tt = tt.ravel()
    uu = uu.ravel()
    zs = np.arange(-0.5, 2.6, z_step)
    base = o[None, :] + tt[:, None] * th[None, :] + uu[:, None] * uh[None, :]
    col = base[:, None, :] + zs[None, :, None] * nh[None, None, :]
    vc = sample_volume(origin, axes, vals, col.reshape(-1, 3)).reshape(len(tt), len(zs))
    return zs, vc


def iso_for_mean_height(vc, zs, target, min_support=0.3):
    """Isovalue whose first-vacuum surface mean height is closest to target."""
    isos = 10 ** np.linspace(-5.5, np.log10(np.nanmax(vc) * 0.8), 80)
    best = None
    best_err = np.inf
    for iso in isos:
        above = vc > iso
        rev = above[:, ::-1]
        first = rev.argmax(axis=1)
        found = rev.max(axis=1)
        z = zs[::-1][first]
        z[~found] = np.nan
        nv = found.sum()
        if nv < min_support * len(z):
            continue
        mh = np.nanmean(z)
        err = abs(mh - target)
        if err < best_err:
            best_err = err
            best = (iso, mh, nv, z)
    return best


def iso_for_mean_height_legacy(vc, zs, target):
    """Champion calibration: first isovalue whose mean height drops below the
    target (no support requirement). Reproduces the 2026-08-02 templates."""
    isos = 10 ** np.linspace(-5.5, np.log10(np.nanmax(vc) * 0.8), 80)
    for iso in isos:
        above = vc > iso
        rev = above[:, ::-1]
        first = rev.argmax(axis=1)
        found = rev.max(axis=1)
        z = zs[::-1][first]
        z[~found] = np.nan
        mh = np.nanmean(z)
        if mh < target:
            return (iso, mh, found.sum(), z)
    return None


def build_templates(cube_path, frame, label, out_path, target=0.50, t_half=0.32,
                    u_half=0.32, step=0.04, legacy=False):
    """Build the z-scored 17x17 (or rectangular) cc template for one type."""
    t_n = round(2 * t_half / step) + 1
    u_n = round(2 * u_half / step) + 1
    t_range = np.arange(-t_half, t_half + 0.001, step)
    u_range = np.arange(-u_half, u_half + 0.001, step)
    zs, vc = surface_grid(cube_path, frame, t_range, u_range)
    if legacy:
        iso, mh, nv, z = iso_for_mean_height_legacy(vc, zs, target)
    else:
        iso, mh, nv, z = iso_for_mean_height(vc, zs, target)
    if z is None:
        raise RuntimeError(f"{label}: no isovalue reaches mean height {target}")
    print(f"{label}: iso={iso:.5g} mean_h={mh:.4f} valid={nv}/{len(z)}")
    z = z.reshape(t_n, u_n)
    zf = np.where(np.isfinite(z), z, np.nan)
    mu = np.nanmean(zf)
    sd = np.nanstd(zf)
    zc = np.where(np.isfinite(zf), (zf - mu) / sd, 0.0)
    variants = [zc, zc[:, ::-1], zc[::-1, :], zc[::-1, ::-1]]
    first = not os.path.exists(out_path)
    with open(out_path, "a" if not first else "w") as f:
        if first:
            cols = ["name", "type", "parity", "mirror"] + [f"p{i:03d}" for i in range(1, t_n * u_n + 1)]
            f.write("\t".join(cols) + "\n")
        for vi, (pv, mv) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
            f.write(f"{label}_p{pv}_m{mv}\t{'0' if label == 'GlcN' else '1'}\t{pv}\t{mv}\t"
                    + "\t".join(f"{v:.7g}" for v in variants[vi].ravel()) + "\n")
    print(f"  template rows written for {label} -> {out_path}")
    return iso


if __name__ == "__main__":
    import os
    import numpy as np  # noqa: F401 (numpy used in module scope)
    if len(sys.argv) < 7:
        print(__doc__)
        sys.exit(1)
    cube0, cube1, frame0, frame1, out = sys.argv[1:6]
    args = sys.argv[6:]
    height = 0.50
    t_half = u_half = 0.32
    step = 0.04
    legacy = False
    i = 0
    while i < len(args):
        if args[i] == "--height":
            height = float(args[i + 1])
            i += 2
        elif args[i] == "--half-nm":
            t_half = u_half = float(args[i + 1])
            i += 2
        elif args[i] == "--half-u-nm":
            u_half = float(args[i + 1])
            i += 2
        elif args[i] == "--legacy":
            legacy = True
            i += 1
        elif args[i] == "--step-nm":
            step = float(args[i + 1])
            i += 2
        else:
            raise SystemExit(f"Unknown arg: {args[i]}")

    def read_frame(p):
        d = {}
        for line in open(p):
            parts = line.strip().split("\t")
            if len(parts) == 2 and parts[0] in ("origin_nm", "t_axis", "u_axis"):
                d[parts[0]] = [float(x) for x in parts[1].split(",")]
        return (d["origin_nm"], d["t_axis"], d["u_axis"])

    if os.path.exists(out):
        os.remove(out)
    build_templates(cube0, read_frame(frame0), "GlcN", out, height, t_half, u_half, step, legacy=legacy)
    build_templates(cube1, read_frame(frame1), "GlcNAc", out, height, t_half, u_half, step, legacy=legacy)
    print(f"templates: {out}")
