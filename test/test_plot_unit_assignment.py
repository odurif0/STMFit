#!/usr/bin/env python3
"""Focused label-free plotting checks; no raw images or external labels needed.

Run: MPLBACKEND=Agg python3 test/test_plot_unit_assignment.py
"""

import csv
import importlib.util
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

os.environ["MPLBACKEND"] = "Agg"

SCRIPT = Path(__file__).with_name("plot_unit_assignment.py")
SPEC = importlib.util.spec_from_file_location("plot_unit_assignment", SCRIPT)
plots = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(plots)


class UnitAssignmentPlotTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.features = self.root / "features.tsv"
        self.predictions = self.root / "predictions.tsv"
        self.features.write_text(
            "file\tlobe\tx_nm\ty_nm\tamplitude\n"
            "scan_A.sxm\t3\t1.2\t0.6\t0.07\n"
            "scan_A.sxm\t1\t0.2\t0.5\t0.05\n"
            "scan_A.sxm\t2\t0.7\t0.6\t0.10\n",
            encoding="utf-8",
        )
        self.predictions.write_text(
            "file\tlobe\tpredicted\tconfidence\n"
            "scan_A.sxm\t1\t0\t0.9\n"
            "scan_A.sxm\t2\t1\t0.8\n"
            "scan_A.sxm\t3\t?\t0.2\n",
            encoding="utf-8",
        )

    def tearDown(self):
        plots.plt.close("all")
        self.tmp.cleanup()

    def run_cli(self, mode):
        output = self.root / mode
        args = [
            str(SCRIPT), "--features", str(self.features),
            "--predictions", str(self.predictions),
            "--out-dir", str(output), "--mode", mode,
        ]
        with mock.patch.object(sys, "argv", args):
            plots.main()
        return output

    def assert_summary(self, output):
        with (output / "summary.tsv").open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["files"], "1")
        self.assertEqual(row["lobes"], "3")
        self.assertEqual(row["predicted_0"], "1")
        self.assertEqual(row["predicted_1"], "1")
        self.assertEqual(row["uncertain"], "1")
        self.assertEqual(row["review_flags"], "has_uncertain")
        self.assertAlmostEqual(float(row["mean_confidence"]), (0.9 + 0.8 + 0.2) / 3)

    def test_single_file_grid_and_generic_title(self):
        with mock.patch.object(plots.plt, "close") as close:
            output = self.run_cli("grid")
        figure = close.call_args.args[0]
        self.assertGreater((output / "summary_grid.png").stat().st_size, 0)
        self.assertEqual(figure._suptitle.get_text(), "Label-free unit-assignment map")
        self.assertEqual(len(figure.axes), 5)
        self.assertEqual(
            [text.get_text() for text in figure.legends[0].get_texts()],
            ["GlcN (0)", "GlcNAc (1)", "uncertain (?)"],
        )
        self.assert_summary(output)

    def test_standalone_preserves_uncertainty_confidence_and_convention(self):
        predictions = plots.load_predictions(self.predictions)
        self.assertEqual(predictions["scan_A.sxm"][3]["predicted"], "?")
        self.assertEqual(predictions["scan_A.sxm"][3]["confidence"], 0.2)
        self.assertEqual(plots.prediction_value(predictions, "scan_A.sxm", 4), "?")
        self.assertEqual(plots.COLORS["?"], "#8c8c8c")
        self.assertEqual(plots.label_color("?"), "black")
        with mock.patch.object(plots.plt, "close") as close:
            output = self.run_cli("standalone")
        figure = close.call_args.args[0]
        axis = figure.axes[0]
        self.assertGreater((output / "standalone" / "scan_A_chain.png").stat().st_size, 0)
        self.assertEqual([text.get_text() for text in axis.texts], ["0", "1", "?"])
        self.assertEqual(
            [text.get_text() for text in axis.get_legend().get_texts()],
            ["GlcN (0)", "GlcNAc (1)", "uncertain (?)"],
        )
        self.assert_summary(output)

    def test_invalid_prediction_is_rejected(self):
        self.predictions.write_text(
            "file\tlobe\tpredicted\tconfidence\n"
            "scan_A.sxm\t1\t2\t0.9\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "Invalid prediction '2'"):
            self.run_cli("standalone")


if __name__ == "__main__":
    unittest.main()
