#!/usr/bin/env python3
"""Build a production 10–20mer multi-pass report.

Inputs are the standard GCV run plus optional support-rescue runs.  All images
are retained.  The script never uses an expected N; it selects between passes by
support completeness and ell/circ consistency, and marks low-confidence cases for
review instead of excluding them.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path
from typing import Any

from PIL import Image, ImageDraw, ImageFont


PASSES = ("standard", "rescue", "aggressive")


def _font_display(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    """Display font for the large N label – Montserrat Bold (geometric, modern)."""
    for path in (
        "/usr/share/fonts/texlive-montserrat/Montserrat-Bold.otf",
        "/usr/share/fonts/texlive-montserrat/Montserrat-SemiBold.otf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    ):
        if Path(path).is_file():
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


def _font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    """Body font for metadata – DejaVu Sans (clean, professional)."""
    for path in (
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/texlive-dejavu/DejaVuSans.ttf",
    ):
        if Path(path).is_file():
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


def _read_summary(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        return {row["filepath"]: row for row in csv.DictReader(handle, delimiter="\t")}


def _as_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "NaN")
    try:
        return float(value)
    except ValueError:
        return math.nan


def _as_int(row: dict[str, str], key: str, default: int | None = None) -> int | None:
    """Parse an integer cell; returns *default* on NA/ERR/unparseable values so
    a single failed file does not abort the whole report."""
    try:
        return int(float(row[key]))
    except (ValueError, TypeError, KeyError):
        return default


def _metrics(row: dict[str, str]) -> dict[str, Any]:
    support_1d = _as_float(row, "support_1D_nm")
    support_2d = _as_float(row, "support_2D_ell_nm")
    # Sentinel -1 for unparseable N (NA/ERR rows from a fully failed fit):
    # downstream selection treats it as an invalid, never-favoured pass.
    n_ell = _as_int(row, "N_ell", -1)
    n_circ = _as_int(row, "N_circ", -1)
    n_sel = _as_int(row, "N_selected", -1)
    ratio = support_2d / support_1d if support_1d > 0 else math.nan
    return {
        "N": n_sel,
        "N_ell": n_ell,
        "N_circ": n_circ,
        "support_1d": support_1d,
        "support_2d": support_2d,
        "support_ratio": ratio,
        "ell_circ_delta": abs(n_ell - n_circ),
        "ambiguous": row.get("ambiguous_eff", "false").lower() == "true",
        "plot": row.get("best_plot", ""),
    }


def _score(metrics: dict[str, Any]) -> float:
    ratio = float(metrics["support_ratio"])
    disagreement = int(metrics["ell_circ_delta"])
    ambiguity = 0.3 if bool(metrics["ambiguous"]) else 0.0
    if not math.isfinite(ratio):
        return -10.0
    # Reward support completeness up to parity with the 1D reference (ratio 1.0),
    # then penalize over-expansion above parity, so the diagnostic "aggressive"
    # pass cannot win purely by inflating support beyond the 1D extent.
    ratio_score = ratio if ratio <= 1.0 else 1.0 - 3.0 * (ratio - 1.0)
    return 3.0 * ratio_score - 0.8 * disagreement - ambiguity


def _choose_pass(pass_metrics: dict[str, dict[str, Any]]) -> str:
    """Pick a pass using label-free support/comparison diagnostics.

    Standard GCV is canonical.  Rescue can override only when standard support is
    objectively truncated and the alternative substantially improves support
    without creating a large ell/circ disagreement.  This prevents global support
    expansion from replacing good standard fits.
    """
    standard = pass_metrics["standard"]
    choice = "standard"
    standard_ratio = float(standard["support_ratio"])

    if math.isfinite(standard_ratio) and standard_ratio < 0.65:
        viable = [
            name
            for name, metrics in pass_metrics.items()
            if float(metrics["support_ratio"]) > standard_ratio + 0.15
            and int(metrics["ell_circ_delta"]) <= 1
        ]
        if viable:
            choice = max(viable, key=lambda name: _score(pass_metrics[name]))

    # Safety: never select a rescue pass whose support collapsed relative to the
    # standard run.
    if choice != "standard":
        if float(pass_metrics[choice]["support_2d"]) < 0.8 * float(standard["support_2d"]):
            choice = "standard"

    return choice


def _confidence(chosen: dict[str, Any], all_n: list[int], chosen_pass: str) -> tuple[str, str]:
    flags: list[str] = []
    support_ratio = float(chosen["support_ratio"])
    n_selected = int(chosen["N"])
    disagreement = int(chosen["ell_circ_delta"])
    spread = max(all_n) - min(all_n)

    if math.isfinite(support_ratio) and support_ratio < 0.65:
        flags.append("support_truncated")
    if n_selected < 8:
        flags.append("low_N")
    if disagreement >= 2:
        flags.append("ell_circ_disagree")
    recovered_by_rescue = chosen_pass != "standard" and n_selected >= 8 and disagreement <= 1 and (
        not math.isfinite(support_ratio) or support_ratio >= 0.65
    )

    if spread >= 4 and not recovered_by_rescue:
        flags.append("support_sensitive_N")
    elif spread >= 4 and recovered_by_rescue:
        flags.append("support_recovered")
    if bool(chosen["ambiguous"]):
        flags.append("ambiguous_gcv")

    if recovered_by_rescue:
        confidence = "medium"
    elif not flags and spread <= 1 and disagreement == 0:
        confidence = "high"
    elif disagreement <= 1 and (not math.isfinite(support_ratio) or support_ratio >= 0.65) and spread <= 3:
        confidence = "medium"
    else:
        confidence = "review"

    return confidence, ";".join(flags) if flags else "ok"


def _guard_effect(standard_n: int, guard_n: int | None, final_n: int) -> tuple[str, str, str, str]:
    if guard_n is None:
        return "NA", "NA", "NA", "guard unavailable"
    if guard_n == standard_n:
        vs_standard = "same"
    else:
        vs_standard = "down" if guard_n < standard_n else "up"
    if guard_n == final_n:
        vs_final = "same_final"
    else:
        vs_final = "diff_final"
    note = f"guard vs standard: {standard_n}->{guard_n} ({vs_standard}); guard vs final: {final_n}->{guard_n} ({vs_final})"
    compact = f"{vs_standard};{vs_final}"
    return vs_standard, vs_final, compact, note


def _fit_text(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.FreeTypeFont | ImageFont.ImageFont, max_w: int) -> str:
    """Truncate *text* with '…' so its rendered width ≤ *max_w* pixels."""
    if draw.textlength(text, font=font) <= max_w:
        return text
    lo, hi = 0, len(text)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        candidate = text[:mid] + "…"
        if draw.textlength(candidate, font=font) <= max_w:
            lo = mid
        else:
            hi = mid - 1
    return text[:lo] + "…" if lo > 0 else "…"


def _annotate_plot(source: Path, destination: Path, lines: list[str], confidence: str, guard_vs_final: str) -> None:
    """Draw a minimal scientific-figure footer below the plot.

    The annotation is placed below the original image, not above it, so it never
    collides with the plot title, axes, legends, or colour bars generated by the
    plotting backend.  Text is clipped to the image width so nothing overflows.

        ┌──────────────────────────────────────────────────────────────┐
        │                      (original plot)                        │
        ╞═════════════ thin accent line (confidence color) ════════════╡
        │  N final = 14  │  251206_013.sxm                           │
        │                │  ● pass standard · confidence HIGH         │
        │                │    GCV 14/14/16 · guard 14 same/diff       │
        │                                                             │
    """
    image = Image.open(source).convert("RGB")
    width, height = image.size
    footer_h = 190
    annotated = Image.new("RGB", (width, height + footer_h), (255, 255, 255))
    annotated.paste(image, (0, 0))
    draw = ImageDraw.Draw(annotated)
    footer_y = height

    # Confidence accent colour (muted, scientific tones)
    accent_map = {"high": (46, 139, 87), "medium": (196, 155, 29), "review": (178, 60, 62)}
    accent = accent_map.get(confidence, (120, 120, 120))

    # Fonts
    n_font = _font_display(64)
    title_font = _font(38)
    meta_font = _font(29)

    # --- Left column: big N label, vertically centred ---
    pad_x = 44
    n_text = lines[0]
    n_bbox = draw.textbbox((0, 0), n_text, font=n_font)
    n_w = n_bbox[2] - n_bbox[0]
    n_ink_h = n_bbox[3] - n_bbox[1]
    n_y = footer_y + (footer_h - n_ink_h) // 2 - n_bbox[1]
    draw.text((pad_x, n_y), n_text, fill=(30, 30, 30), font=n_font)

    # --- Vertical separator ---
    sep_x = pad_x + n_w + 34
    draw.line((sep_x, footer_y + 24, sep_x, footer_y + footer_h - 24), fill=(210, 210, 210), width=2)

    # --- Right column: filename + metadata ---
    text_x = sep_x + 28
    right_w = int(width - text_x - 24)     # pixels available for right-column text

    # Filename (line 1)
    title_y = footer_y + 24
    draw.text((text_x, title_y), _fit_text(draw, lines[1], title_font, right_w), fill=(20, 20, 20), font=title_font)

    # Confidence dot + status (line 2)
    meta_y = title_y + 58
    dot_r = 6
    dot_cy = meta_y + 14
    draw.ellipse(
        (text_x, dot_cy - dot_r, text_x + 2 * dot_r, dot_cy + dot_r), fill=accent
    )
    meta_text_x = text_x + 2 * dot_r + 10
    meta_w = int(width - meta_text_x - 24)
    draw.text((meta_text_x, meta_y), _fit_text(draw, lines[2], meta_font, meta_w), fill=(100, 100, 100), font=meta_font)

    # GCV / guard comparison (line 3)
    draw.text((meta_text_x, meta_y + 44), _fit_text(draw, lines[3], meta_font, meta_w), fill=(100, 100, 100), font=meta_font)

    # --- Footer separator/accent line ---
    draw.line((0, footer_y, width, footer_y), fill=accent, width=3)

    annotated.save(destination)


def build_report(args: argparse.Namespace) -> list[dict[str, str]]:
    summaries = {
        "standard": _read_summary(Path(args.standard)),
        "rescue": _read_summary(Path(args.rescue)),
        "aggressive": _read_summary(Path(args.aggressive)),
    }
    guard_rows = _read_summary(Path(args.guard)) if args.guard else {}
    files = sorted(summaries["standard"])
    for pass_name, rows in summaries.items():
        missing = set(files) - set(rows)
        if missing:
            raise SystemExit(f"{pass_name} summary is missing files: {sorted(missing)}")

    output_dir = Path(args.output_dir)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    for stale_plot in plots_dir.glob("*.png"):
        stale_plot.unlink()

    final_rows: list[dict[str, str]] = []
    for file in files:
        pass_metrics = {name: _metrics(rows[file]) for name, rows in summaries.items()}
        # A fully failed standard fit (N_selected == "ERR"/"NA", coerced to the
        # -1 sentinel) is retained as an explicit error row instead of aborting
        # the whole report.
        if int(pass_metrics["standard"]["N"]) < 0:
            final_rows.append({
                "file": file,
                "N_final": "ERR",
                "chosen_pass": "standard",
                "confidence": "error",
                "flags": "standard_fit_failed",
                "guard_N": "NA",
                "guard_vs_standard": "NA",
                "guard_vs_final": "NA",
                "guard_effect": "NA",
                "guard_note": "standard pass produced no fit",
                "standard_N": "ERR",
                "rescue_N": "ERR",
                "aggressive_N": "ERR",
                "standard_support_nm": "NA",
                "rescue_support_nm": "NA",
                "aggressive_support_nm": "NA",
                "support_ratio_final": "NA",
                "N_ell_final": "NA",
                "N_circ_final": "NA",
                "plot": "",
            })
            continue
        choice = _choose_pass(pass_metrics)
        chosen = pass_metrics[choice]
        all_n = [int(pass_metrics[name]["N"]) for name in PASSES]
        confidence, flags = _confidence(chosen, all_n, choice)
        final_n = int(chosen["N"])
        standard_n = int(pass_metrics["standard"]["N"])
        guard_n = _as_int(guard_rows[file], "N_selected") if file in guard_rows else None
        guard_vs_standard, guard_vs_final, guard_effect, guard_note = _guard_effect(standard_n, guard_n, final_n)

        source_plot = Path(str(chosen["plot"]))
        final_plot = ""
        if source_plot.is_file():
            final_plot_path = plots_dir / source_plot.name
            guard_str = "NA" if guard_n is None else str(guard_n)
            guard_label = "same" if guard_vs_final == "same_final" else "DIFF"
            trio = (
                f"{pass_metrics['standard']['N']}/"
                f"{pass_metrics['rescue']['N']}/"
                f"{pass_metrics['aggressive']['N']}"
            )
            big_n = f"N final = {final_n}"
            title = file
            status = f"pass {choice} · confidence {confidence.upper()} · flags {flags}"
            model_line = f"GCV standard/rescue/aggressive = {trio} · guard = {guard_str} ({guard_label})"
            _annotate_plot(source_plot, final_plot_path, [big_n, title, status, model_line], confidence, guard_vs_final)
            final_plot = str(final_plot_path)

        final_rows.append(
            {
                "file": file,
                "N_final": str(chosen["N"]),
                "chosen_pass": choice,
                "confidence": confidence,
                "flags": flags,
                "guard_N": "NA" if guard_n is None else str(guard_n),
                "guard_vs_standard": guard_vs_standard,
                "guard_vs_final": guard_vs_final,
                "guard_effect": guard_effect,
                "guard_note": guard_note,
                "standard_N": str(pass_metrics["standard"]["N"]),
                "rescue_N": str(pass_metrics["rescue"]["N"]),
                "aggressive_N": str(pass_metrics["aggressive"]["N"]),
                "standard_support_nm": f"{float(pass_metrics['standard']['support_2d']):.3f}",
                "rescue_support_nm": f"{float(pass_metrics['rescue']['support_2d']):.3f}",
                "aggressive_support_nm": f"{float(pass_metrics['aggressive']['support_2d']):.3f}",
                "support_ratio_final": f"{float(chosen['support_ratio']):.3f}",
                "N_ell_final": str(chosen["N_ell"]),
                "N_circ_final": str(chosen["N_circ"]),
                "plot": final_plot,
            }
        )

    return final_rows


def write_outputs(rows: list[dict[str, str]], output_dir: Path) -> None:
    fields = [
        "file",
        "N_final",
        "chosen_pass",
        "confidence",
        "flags",
        "guard_N",
        "guard_vs_standard",
        "guard_vs_final",
        "guard_effect",
        "guard_note",
        "standard_N",
        "rescue_N",
        "aggressive_N",
        "standard_support_nm",
        "rescue_support_nm",
        "aggressive_support_nm",
        "support_ratio_final",
        "N_ell_final",
        "N_circ_final",
        "plot",
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "final_results.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with (output_dir / "final_results.md").open("w") as handle:
        handle.write("# 10–20mer final multi-pass results\n\n")
        handle.write(
            "Policy: GCV/N_eff with n_max=24. Robust down-guard disabled for final selection. "
            "All images retained. Support-rescue passes are used only when the standard support "
            "is objectively truncated. The robust guard is reported as an audit-only diagnostic. "
            "No expected N is used.\n\n"
        )
        handle.write("| file | N_final | pass | confidence | flags | guard vs final | guard vs standard | std/rescue/aggr N | plot |\n")
        handle.write("|---|---:|---|---|---|---|---|---|---|\n")
        for row in rows:
            trio = f"{row['standard_N']}/{row['rescue_N']}/{row['aggressive_N']}"
            guard_final = f"{row['guard_N']} ({row['guard_vs_final']})"
            guard_standard = f"{row['guard_N']} ({row['guard_vs_standard']})"
            handle.write(
                f"| {row['file']} | {row['N_final']} | {row['chosen_pass']} | "
                f"{row['confidence']} | {row['flags']} | {guard_final} | {guard_standard} | {trio} | `{row['plot']}` |\n"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--standard", default="results/10_20mer_analysis_nmax24_gcv/summary_overlap060_hard.tsv")
    parser.add_argument("--rescue", default="results/10_20mer_analysis_rescue/summary_overlap060_hard.tsv")
    parser.add_argument("--aggressive", default="results/10_20mer_analysis_rescue_aggressive/summary_overlap060_hard.tsv")
    parser.add_argument("--guard", default="results/10_20mer_analysis_guard_audit/summary_overlap060_hard.tsv")
    parser.add_argument("--output-dir", default="results/10_20mer_analysis_final")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = build_report(args)
    output_dir = Path(args.output_dir)
    write_outputs(rows, output_dir)
    print("n_rows", len(rows))
    print("N_final_distribution", dict(sorted(Counter(int(row["N_final"]) for row in rows if row["N_final"].lstrip("-").isdigit()).items())))
    print("confidence", dict(Counter(row["confidence"] for row in rows)))
    print(output_dir / "final_results.tsv")
    print(output_dir / "final_results.md")


if __name__ == "__main__":
    main()
