#!/usr/bin/env python3

"""Render simple static XYZ previews for QE mold structures.

This is a visualization helper only. It does not modify geometries and does not
enter fitting, scoring, or selection.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt


COLORS = {
    "H": "#f2f2f2",
    "C": "#303030",
    "N": "#2f5eff",
    "O": "#e3342f",
    "Cu": "#b87333",
}

RADII = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "Cu": 1.32,
}

SIZES = {
    "H": 18,
    "C": 42,
    "N": 48,
    "O": 48,
    "Cu": 18,
}


def read_xyz(path: Path):
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"XYZ too short: {path}")
    n = int(lines[0].strip())
    atoms = []
    for i, line in enumerate(lines[2 : 2 + n], start=1):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Bad XYZ atom line {i}: {line}")
        label = parts[4] if len(parts) >= 5 else f"{parts[0]}{i}"
        atoms.append(
            {
                "element": parts[0],
                "x": float(parts[1]),
                "y": float(parts[2]),
                "z": float(parts[3]),
                "label": label,
            }
        )
    return atoms


def dist(a, b):
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def infer_bonds(atoms):
    bonds = []
    for i, a in enumerate(atoms[:-1]):
        for j, b in enumerate(atoms[i + 1 :], start=i + 1):
            if a["element"] == "Cu" and b["element"] == "Cu":
                continue
            cutoff = 1.22 * (RADII.get(a["element"], 0.75) + RADII.get(b["element"], 0.75))
            d = dist(a, b)
            if 0.45 <= d <= min(cutoff, 2.15):
                bonds.append((i, j))
    return bonds


def projected(atom, axes):
    return atom[axes[0]], atom[axes[1]]


def render_panel(ax, atoms, bonds, axes, title, molecule_only=False):
    draw_atoms = [a for a in atoms if not molecule_only or a["element"] != "Cu"]
    draw_set = {id(a) for a in draw_atoms}
    for i, j in bonds:
        a, b = atoms[i], atoms[j]
        if id(a) not in draw_set or id(b) not in draw_set:
            continue
        x1, y1 = projected(a, axes)
        x2, y2 = projected(b, axes)
        ax.plot([x1, x2], [y1, y2], color="#6b7280", linewidth=0.8, alpha=0.75, zorder=1)

    for elem in sorted({a["element"] for a in draw_atoms}):
        group = [a for a in draw_atoms if a["element"] == elem]
        xs = [a[axes[0]] for a in group]
        ys = [a[axes[1]] for a in group]
        ax.scatter(
            xs,
            ys,
            s=SIZES.get(elem, 35),
            c=COLORS.get(elem, "#999999"),
            edgecolors="#111827",
            linewidths=0.3,
            label=elem,
            zorder=2,
        )
    ax.set_title(title)
    ax.set_xlabel(f"{axes[0]} / Å")
    ax.set_ylabel(f"{axes[1]} / Å")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.18)


def render(path: Path, out: Path, title: str, molecule_only: bool):
    atoms = read_xyz(path)
    bonds = infer_bonds(atoms)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    render_panel(axes[0], atoms, bonds, ("x", "y"), "top: x/y", molecule_only)
    render_panel(axes[1], atoms, bonds, ("x", "z"), "side: x/z", molecule_only)
    render_panel(axes[2], atoms, bonds, ("y", "z"), "side: y/z", molecule_only)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=max(1, len(labels)))
    fig.suptitle(title, fontsize=14)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xyz", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--title", default="XYZ preview")
    parser.add_argument("--molecule-only", action="store_true")
    args = parser.parse_args()
    render(args.xyz, args.out, args.title, args.molecule_only)


if __name__ == "__main__":
    main()
