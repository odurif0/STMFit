#!/usr/bin/env python3
"""Reproduce the 2026-08-02 label-free unit-assignment champion.

Champion: soft vote (mean of probabilities) of
  A) k-means 4-view (+ interactions, 20 seeds) and
  B) GMM 1-view + per-channel constant-current margins + Fisher-mold
     margin + self-training 2.
The empirical mold is the Fisher discriminant: GMM-on-PCA10 cluster means
and the regularized noise covariance (label-free, amplitude convention),
half-split cross-validated (see test/lib/empirical_fisher_mold.py).
Score: 79.3% classified physical accuracy / 36 exact chains / 677 of 854
honest - the promotion bar (78.9% / 18 / 677) is MET (cross-validated).

Label-free: no truth, sequence, or composition prior anywhere (audited:
zero label references in all construction scripts). All grading happens only
in the final post-hoc report step.

Inputs (extracted on Viper, see docs/src/unit_assignment.md):
  --features      per-lobe feature TSV (full145, e.g. full_features_origN_bwd_uasym.tsv)
  --patches-fwd   forward residual 17x17 patch TSV (step 0.04, half 0.32)
  --patches-bwd   backward residual 17x17 patch TSV (same grid)
  --cube0/--cube1 converged GlcN / GlcNAc LDOS cubes (bohr)
  --frame0/--frame1 relaxed-geometry frame TSVs
  --out           final prediction TSV
  --workdir       scratch dir [results/champion_cc_soft]

Steps: cc templates (test/lib/cc_mold_builder.py --legacy) -> score both
channels (test/score_connected_mold_templates.jl) -> Fisher mold
(test/lib/empirical_fisher_mold.py, half-split CV margins) -> feature table
with mold_cc_fwd/bwd + emp_fisher margins -> GMM 1-view chan st2 ->
k-means 4-view -> soft vote -> post-hoc grade.
"""

import argparse
import csv
import math
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUILDER = os.path.join(ROOT, "test", "lib", "cc_mold_builder.py")
SCORER = os.path.join(ROOT, "test", "score_connected_mold_templates.jl")
GMM = os.path.join(ROOT, "test", "build_labelfree_gmm_predictions.jl")
KM = os.path.join(ROOT, "test", "build_labelfree_unit_predictions.jl")
GRADER = os.path.join(ROOT, "test", "report_unit_assignment_benchmark.jl")
BASE4 = "amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"


def run(cmd):
    print("+", " ".join(str(c) for c in cmd))
    subprocess.run([str(c) for c in cmd], check=True, cwd=ROOT)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features", required=True)
    ap.add_argument("--patches-fwd", required=True)
    ap.add_argument("--patches-bwd", required=True)
    ap.add_argument("--cube0", required=True)
    ap.add_argument("--cube1", required=True)
    ap.add_argument("--frame0", required=True)
    ap.add_argument("--frame1", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workdir", default=os.path.join(ROOT, "results", "champion_cc_soft"))
    args = ap.parse_args()

    wd = args.workdir
    os.makedirs(wd, exist_ok=True)
    templates = os.path.join(wd, "templates_cc_17.tsv")
    score_fwd = os.path.join(wd, "score_fwd.tsv")
    score_bwd = os.path.join(wd, "score_bwd.tsv")
    table = os.path.join(wd, "features_cc.tsv")
    pred_gmm = os.path.join(wd, "pred_gmm.tsv")
    pred_km = os.path.join(wd, "pred_km.tsv")

    # 1. adaptive-contour templates (legacy calibration = champion calibration)
    run(["python3", BUILDER, args.cube0, args.cube1, args.frame0, args.frame1,
         templates, "--height", "0.50", "--half-nm", "0.32", "--step-nm", "0.04",
         "--legacy"])

    # 2. score both channels against the templates (contrast mode)
    run(["julia", "--project=.", SCORER, "--patches", args.patches_fwd, "--templates",
         templates, "--template-mode", "contrast", "--prefix", "res", "--out", score_fwd])
    run(["julia", "--project=.", SCORER, "--patches", args.patches_bwd, "--templates",
         templates, "--template-mode", "contrast", "--prefix", "bwd_res", "--out", score_bwd])

    # 2.5 Fisher empirical mold on the fwd residual patches (half-split CV
    # margins: the mold is trained on one half and applied to the other).
    fisher = os.path.join(wd, "emp_fisher")
    run(["python3", os.path.join(ROOT, "test", "lib", "empirical_fisher_mold.py"),
         args.patches_fwd, "res", fisher])

    # 3. feature table: base features + per-channel cc margins + Fisher margin
    margins = {}
    for path, tag in ((score_fwd, "fwd"), (score_bwd, "bwd")):
        with open(path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                margins.setdefault((r["file"], int(r["lobe"])), {})[tag] = float(r["cost_margin"])
    fisher_m = {}
    with open(fisher + "_cv.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            fisher_m[(r["file"], int(r["lobe"]))] = -float(r["score"])
    with open(args.features) as f:
        feats = list(csv.DictReader(f, delimiter="\t"))
    cols = list(feats[0].keys()) + ["mold_cc_fwd", "mold_cc_bwd", "emp_fisher"]
    if "skew_ratio" in feats[0] and "split_log_skew" not in feats[0]:
        cols.append("split_log_skew")
    with open(table, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in feats:
            m = margins.get((r["file"], int(r["lobe"])))
            if m:
                r["mold_cc_fwd"] = f"{m['fwd']:.8g}"
                r["mold_cc_bwd"] = f"{m['bwd']:.8g}"
            fm = fisher_m.get((r["file"], int(r["lobe"])))
            if fm is not None:
                r["emp_fisher"] = f"{fm:.8g}"
            if "split_log_skew" in cols:
                try:
                    sk = float(r.get("skew_ratio", "nan"))
                except ValueError:
                    sk = float("nan")
                r["split_log_skew"] = f"{math.log(sk):.8g}" if sk > 0 else "NA"
            w.writerow(r)

    # 4. GMM 1-view + per-channel cc margins + Fisher margin + self-training 2
    run(["julia", "--project=.", GMM, "--features", table, "--out", pred_gmm,
         "--view", f"v_cc={BASE4},patch_u_asym,mold_cc_fwd,mold_cc_bwd,emp_fisher",
         "--seeds", "10", "--interactions", "--selftrain", "2"])

    # 5. k-means 4-view (base, base+split, base+com_t, base+diag45) + interactions
    run(["julia", "--project=.", KM, "--features", table, "--out", pred_km,
         "--view", f"v_base={BASE4}",
         "--view", f"v_split={BASE4},split_log_skew",
         "--view", f"v_comt={BASE4},bwd_neg_com_t",
         "--view", f"v_diag45={BASE4},bwd_neg_diag45",
         "--seeds", "20", "--interactions"])

    # 6. soft vote (mean of probabilities), drop non-manifest file
    def load(p):
        out = {}
        with open(p) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                try:
                    p1 = float(r["probability_1"])
                except (ValueError, KeyError):
                    p1 = 0.5
                out[(r["file"], int(r["lobe"]))] = p1
        return out

    km_p = load(pred_km)
    gmm_p = load(pred_gmm)
    with open(args.out, "w") as g:
        g.write("file\tlobe\tpredicted\tconfidence\n")
        for (file, lobe) in sorted(set(km_p) & set(gmm_p)):
            if file == "240310_Cu100009.sxm":
                continue
            p = (km_p[(file, lobe)] + gmm_p[(file, lobe)]) / 2
            g.write(f"{file}\t{lobe}\t{1 if p >= 0.5 else 0}\t{abs(p-0.5)*2:.8f}\n")
    print(f"champion predictions: {args.out}")

    # 7. post-hoc grade
    run(["julia", "--project=.", GRADER, "--full145-own-n",
         "--profile", f"champion={args.out}",
         "--outdir", os.path.join(wd, "grade")])


if __name__ == "__main__":
    main()
