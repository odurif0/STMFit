#!/usr/bin/env python3
"""Build the frozen 2026-08-02 label-free unit-assignment champion.

Champion: soft vote (mean of probabilities) of
  A) k-means 4-view (+ interactions, 20 seeds) and
  B) GMM 1-view + per-channel constant-current margins + Fisher-mold
     margin + self-training 2.
The empirical mold is the Fisher discriminant: GMM-on-PCA10 cluster means
and the regularized noise covariance (label-free, amplitude convention),
half-split cross-validated (see test/lib/empirical_fisher_mold.py).

Prediction only: no grading, control inputs, file exclusions, or composition
prior. Use the separate benchmark report command after predictions exist.
All feature keys are retained; unavailable component predictions emit '?'.
Confidence is the frozen soft-vote margin, not a calibrated probability.

Inputs:
  --features      selected-N per-lobe features, including split/bwd descriptors
  --patches-fwd   forward residual 17x17 patch TSV (step 0.04, half 0.32)
  --patches-bwd   backward residual 17x17 patch TSV (same grid)
  --cube0/--cube1 converged GlcN / GlcNAc LDOS cubes (bohr)
  --frame0/--frame1 relaxed-geometry frame TSVs
  --out           final prediction TSV (must not exist)
  --workdir       new directory for intermediate outputs (must not exist)
  --julia         Julia executable [julia; project runtime: 1.13]

Steps: frozen constant-current templates -> score both channels -> Fisher
half-split CV margins -> feature table -> GMM and k-means -> soft vote.
Scientific settings are unchanged from the frozen champion.
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
BASE4 = "amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"


def run(cmd):
    print("+", " ".join(str(c) for c in cmd))
    subprocess.run([str(c) for c in cmd], check=True, cwd=ROOT)


def read_lobe_table(path, required=("file", "lobe")):
    """Check production columns and keys without using an expected count."""
    with open(path, newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        columns = reader.fieldnames or []
        if len(columns) != len(set(columns)):
            raise ValueError(f"Duplicate columns in {path}")
        for column in columns:
            normalized = "".join(c for c in column.lower() if c.isalnum())
            if (normalized in {"sequence", "controlsequence", "expectedn", "targetn",
                               "manifest", "grade", "report"}
                    or any(token in normalized for token in ("benchmark", "truth"))
                    or normalized.endswith("sequence")):
                raise ValueError(f"Forbidden benchmark/control column: {column}")
        missing = set(required) - set(columns)
        if missing:
            raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty lobe table: {path}")
    keyed = {}
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Malformed TSV row in {path}")
        file = os.path.basename(row["file"].strip())
        lobe = int(row["lobe"])
        key = (file, lobe)
        if not file or lobe < 1:
            raise ValueError(f"Invalid lobe key in {path}: {key}")
        if key in keyed:
            raise ValueError(f"Duplicate lobe key in {path}: {key}")
        row["file"] = file
        keyed[key] = row
    return columns, keyed


def write_soft_vote(features, pred_km, pred_gmm, out):
    """Keep the frozen vote for valid pairs and retain unavailable lobes as '?'."""
    _, reference = read_lobe_table(features)
    _, km_rows = read_lobe_table(pred_km, ("file", "lobe", "predicted", "probability_1"))
    _, gmm_rows = read_lobe_table(pred_gmm, ("file", "lobe", "predicted", "probability_1"))
    for component in (km_rows, gmm_rows):
        if set(component) - set(reference):
            raise ValueError("Component predictions contain keys absent from features")

    def probability(rows, key):
        row = rows.get(key)
        if row is None or row["predicted"] == "?":
            return None
        if row["predicted"] not in ("0", "1"):
            raise ValueError(f"Invalid component prediction for {key}")
        try:
            value = float(row["probability_1"])
        except ValueError:
            return None
        if not math.isfinite(value):
            return None
        if not 0 <= value <= 1:
            raise ValueError(f"Component probability outside [0, 1] for {key}")
        return value

    output_rows = []
    for key in sorted(reference):
        pk, pg = probability(km_rows, key), probability(gmm_rows, key)
        missing = [name for name, value in (("kmeans", pk), ("gmm", pg)) if value is None]
        if missing:
            output_rows.append((*key, "?", "0.00000000", "NA",
                                "missing_component:" + ",".join(missing)))
        else:
            p = (pk + pg) / 2
            output_rows.append((*key, str(int(p >= 0.5)), f"{abs(p-0.5)*2:.8f}",
                                f"{p:.8f}", "ok"))
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("file", "lobe", "predicted", "confidence", "probability_1",
                         "invalid_reason"))
        writer.writerows(output_rows)


def main(argv=None):
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
    ap.add_argument("--workdir", required=True, help="New directory for intermediate outputs")
    ap.add_argument("--julia", default="julia", help="Julia executable (project runtime: 1.13)")
    args = ap.parse_args(argv)

    # Fail before computation, and never overwrite an earlier run.
    for path in (args.features, args.patches_fwd, args.patches_bwd):
        read_lobe_table(path)
    for path in (args.cube0, args.cube1, args.frame0, args.frame1):
        if not os.path.isfile(path):
            raise ValueError(f"Missing mold input: {path}")
    if os.path.lexists(args.out):
        raise FileExistsError(f"Output already exists: {args.out}")
    wd = args.workdir
    os.makedirs(wd, exist_ok=False)
    templates = os.path.join(wd, "templates_cc_17.tsv")
    score_fwd = os.path.join(wd, "score_fwd.tsv")
    score_bwd = os.path.join(wd, "score_bwd.tsv")
    table = os.path.join(wd, "features_cc.tsv")
    pred_gmm = os.path.join(wd, "pred_gmm.tsv")
    pred_km = os.path.join(wd, "pred_km.tsv")

    # 1. adaptive-contour templates (legacy calibration = champion calibration)
    run([sys.executable, BUILDER, args.cube0, args.cube1, args.frame0, args.frame1,
         templates, "--height", "0.50", "--half-nm", "0.32", "--step-nm", "0.04",
         "--legacy"])

    # 2. score both channels against the templates (contrast mode)
    run([args.julia, "--project=.", SCORER, "--patches", args.patches_fwd, "--templates",
         templates, "--template-mode", "contrast", "--prefix", "res", "--out", score_fwd])
    run([args.julia, "--project=.", SCORER, "--patches", args.patches_bwd, "--templates",
         templates, "--template-mode", "contrast", "--prefix", "bwd_res", "--out", score_bwd])

    # 2.5 Fisher empirical mold on the fwd residual patches (half-split CV
    # margins: the mold is trained on one half and applied to the other).
    fisher = os.path.join(wd, "emp_fisher")
    run([sys.executable, os.path.join(ROOT, "test", "lib", "empirical_fisher_mold.py"),
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
    run([args.julia, "--project=.", GMM, "--features", table, "--out", pred_gmm,
         "--view", f"v_cc={BASE4},patch_u_asym,mold_cc_fwd,mold_cc_bwd,emp_fisher",
         "--seeds", "10", "--interactions", "--selftrain", "2"])

    # 5. k-means 4-view (base, base+split, base+com_t, base+diag45) + interactions
    run([args.julia, "--project=.", KM, "--features", table, "--out", pred_km,
         "--view", f"v_base={BASE4}",
         "--view", f"v_split={BASE4},split_log_skew",
         "--view", f"v_comt={BASE4},bwd_neg_com_t",
         "--view", f"v_diag45={BASE4},bwd_neg_diag45",
         "--seeds", "20", "--interactions"])

    # 6. Soft vote on the full selected-N key set; grading is a separate command.
    write_soft_vote(args.features, pred_km, pred_gmm, args.out)
    print(f"champion predictions: {args.out}")

if __name__ == "__main__":
    main()
