#!/usr/bin/env python3
"""Plot unit-assignment predictions overlaid on STM chain fits.

Creates two outputs per file:
  1. Standalone chain diagram (x_nm, y_nm scatter, 0/1 colored)
  2. Annotated overlay on the existing best-fit PNG

Usage:
  python3 test/plot_unit_assignment.py \\
      --features results/unit_separability/lobe_features_selectedN_primary.tsv \\
      --predictions results/unit_assignment/stm_prelim_h050_half032_predictions.tsv \\
      --plots-dir results/best_plots \\
      --out-dir results/unit_assignment/plots_h050

Label-free: reads only predictions + positions. Truth is never used.
"""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import numpy as np

COLORS = {0: "#2166ac", 1: "#b2182b", "?": "#8c8c8c"}
LABELS = {0: "GlcN (0)", 1: "GlcNAc (1)", "?": "uncertain (?)"}


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_features(path):
    """file -> [{lobe, x, y, amp, axis_x, axis_y, origin_x, origin_y}]"""
    feats = {}
    for row in read_tsv(path):
        f = row["file"]
        if f not in feats:
            feats[f] = {"lobes": [], "axis": None}
        feats[f]["lobes"].append({
            "lobe": int(row["lobe"]),
            "x": float(row["x_nm"]),
            "y": float(row["y_nm"]),
            "amp": float(row["amplitude"]),
        })
        ax = float(row.get("axis_x", 0) or 0)
        ay = float(row.get("axis_y", 0) or 0)
        ox = float(row.get("origin_x_nm", 0) or 0)
        oy = float(row.get("origin_y_nm", 0) or 0)
        feats[f]["axis"] = (ax, ay, ox, oy)
    return feats


def load_predictions(path):
    preds = {}
    for row in read_tsv(path):
        f = row["file"]
        if f not in preds:
            preds[f] = {}
        pred = row["predicted"].strip()
        if pred == "?":
            parsed = "?"
        elif pred in ("0", "1"):
            parsed = int(pred)
        else:
            raise ValueError(f"Invalid prediction {pred!r} in {path}: {f} lobe {row['lobe']}")
        conf_raw = row.get("confidence", "").strip()
        confidence = float(conf_raw) if conf_raw and conf_raw not in ("NA", "NaN", "nan") else None
        preds[f][int(row["lobe"])] = {"predicted": parsed, "confidence": confidence}
    return preds


def prediction_value(preds, f_name, lobe):
    entry = preds[f_name].get(lobe)
    return "?" if entry is None else entry["predicted"]


def label_color(pred):
    return "black" if pred == "?" else "white"


def legend_handles():
    return [mpatches.Patch(color=COLORS[k], label=LABELS[k]) for k in (0, 1, "?")]


def write_summary(preds, files, out_path):
    counts = {0: 0, 1: 0, "?": 0}
    confidences = []
    for f_name in files:
        for entry in preds[f_name].values():
            counts[entry["predicted"]] += 1
            if entry["confidence"] is not None:
                confidences.append(entry["confidence"])
    lobes = counts[0] + counts[1] + counts["?"]
    mean_conf = sum(confidences) / len(confidences) if confidences else "NA"
    review_flags = "has_uncertain" if counts["?"] else "ok"
    if out_path.is_symlink():
        raise RuntimeError(f"Refusing to overwrite symlink: {out_path}")
    with open(out_path, "w") as f:
        f.write("files\tlobes\tpredicted_0\tpredicted_1\tuncertain\tmean_confidence\treview_flags\n")
        f.write(f"{len(files)}\t{lobes}\t{counts[0]}\t{counts[1]}\t{counts['?']}\t{mean_conf}\t{review_flags}\n")


def plot_chain_standalone(feats, preds, f_name, out_path):
    """Scatter plot of chain lobes with 0/1 color coding."""
    if f_name not in feats or f_name not in preds:
        return
    lobes = sorted(feats[f_name]["lobes"], key=lambda l: l["lobe"])
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))

    xs = [l["x"] for l in lobes]
    ys = [l["y"] for l in lobes]
    ax.plot(xs, ys, "k-", alpha=0.3, linewidth=1)

    for l in lobes:
        p = prediction_value(preds, f_name, l["lobe"])
        c = COLORS.get(p, "#999999")
        size = 50 + 500 * l["amp"] / max(x["amp"] for x in lobes)
        ax.scatter(l["x"], l["y"], c=c, s=size, zorder=5, edgecolors="white",
                   linewidths=1.5)
        ax.annotate(str(p), (l["x"], l["y"]), fontsize=9, fontweight="bold",
                    color=label_color(p), ha="center", va="center", zorder=6)

    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_title(f"{f_name}  N={len(lobes)}", fontsize=10)
    ax.set_aspect("equal")
    ax.invert_yaxis()  # STM convention: y increases downward
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_overlay(feats, preds, f_name, png_path, out_path, scan_nm=10.0):
    """Overlay 0/1 labels on existing best-fit PNG."""
    if f_name not in feats or f_name not in preds:
        return
    if not Path(png_path).exists():
        return
    lobes = sorted(feats[f_name]["lobes"], key=lambda l: l["lobe"])
    img = mpimg.imread(png_path)
    h, w = img.shape[:2]
    n_panels = w // h  # typically 3 panels of 800x800
    panel_w = h  # square panel

    fig, ax = plt.subplots(1, 1, figsize=(12, 4))
    ax.imshow(img, extent=(0.0, float(w), float(h), 0.0))  # y flipped for image
    ax.axis("off")

    # Map (x_nm, y_nm) to pixel coords in left panel
    # Assume scan_nm maps to panel_w pixels, origin at top-left of left panel
    for l in lobes:
        p = prediction_value(preds, f_name, l["lobe"])
        c = COLORS.get(p, "#999999")
        px = (l["x"] / scan_nm) * panel_w
        py = (l["y"] / scan_nm) * panel_w
        size = 80 + 400 * l["amp"] / max(x["amp"] for x in lobes)
        ax.scatter(px, py, c=c, s=size, zorder=5, edgecolors="yellow",
                   linewidths=2)
        ax.annotate(str(p), (px, py), fontsize=10, fontweight="bold",
                    color=label_color(p), ha="center", va="center", zorder=6)

    # Legend
    ax.legend(handles=legend_handles(), loc="upper right", fontsize=9, framealpha=0.9)
    ax.set_title(f"{f_name}", fontsize=10, loc="left")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_summary_grid(feats, preds, out_path, n_cols=5):
    """Grid of all chains with 0/1 predictions."""
    files = sorted(set(feats.keys()) & set(preds.keys()))
    n = len(files)
    n_rows = (n + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 2.8))
    axes = np.array(axes).reshape(-1) if n > 1 else np.array([axes])

    for idx, f_name in enumerate(files):
        ax = axes[idx]
        lobes = sorted(feats[f_name]["lobes"], key=lambda l: l["lobe"])
        if not lobes:
            ax.axis("off")
            continue
        xs = [l["x"] for l in lobes]
        ys = [l["y"] for l in lobes]
        ax.plot(xs, ys, "k-", alpha=0.2, linewidth=0.5)
        max_amp = max(l["amp"] for l in lobes)
        for l in lobes:
            p = prediction_value(preds, f_name, l["lobe"])
            c = COLORS.get(p, "#999")
            size = 15 + 80 * l["amp"] / max_amp
            ax.scatter(l["x"], l["y"], c=c, s=size, zorder=5, edgecolors="white",
                       linewidths=0.5)
        ax.set_title(f_name.replace("240817_", "").replace(".sxm", ""),
                     fontsize=7)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.tick_params(labelsize=5)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

    for idx in range(n, len(axes)):
        axes[idx].axis("off")

    fig.legend(handles=legend_handles(), loc="lower center", ncol=3, fontsize=9)
    fig.suptitle("Unit assignment — preliminary DFT-STM molds h=0.50 nm", fontsize=11)
    plt.tight_layout(rect=(0.0, 0.03, 1.0, 0.97))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--features", required=True)
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--plots-dir", default="results/best_plots")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--scan-nm", type=float, default=10.0)
    ap.add_argument("--mode", choices=["all", "grid", "standalone", "overlay"],
                    default="all")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    feats = load_features(args.features)
    preds = load_predictions(args.predictions)
    files = sorted(set(feats.keys()) & set(preds.keys()))
    print(f"Features: {len(feats)} files, Predictions: {len(preds)} files")
    print(f"Common: {len(files)} files")
    summary_path = out / "summary.tsv"
    write_summary(preds, files, summary_path)
    print(f"Summary: {summary_path}")

    if args.mode in ("all", "grid"):
        grid_path = out / "summary_grid.png"
        plot_summary_grid(feats, preds, grid_path)
        print(f"Grid: {grid_path}")

    if args.mode in ("all", "standalone"):
        sd = out / "standalone"
        sd.mkdir(exist_ok=True)
        for f in files:
            plot_chain_standalone(feats, preds, f, sd / f.replace(".sxm", "_chain.png"))
        print(f"Standalone: {sd}/ ({len(files)} files)")

    if args.mode in ("all", "overlay"):
        od = out / "overlay"
        od.mkdir(exist_ok=True)
        for f in files:
            png = Path(args.plots_dir) / f.replace(".sxm", "_best.png")
            plot_overlay(feats, preds, f, png, od / f.replace(".sxm", "_overlay.png"),
                         scan_nm=args.scan_nm)
        print(f"Overlay: {od}/ ({len(files)} files)")


if __name__ == "__main__":
    main()
