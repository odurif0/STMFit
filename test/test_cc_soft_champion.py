#!/usr/bin/env python3
"""Focused tests for prediction-only frozen champion wiring (no model refits)."""
import csv
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).with_name("build_cc_soft_champion.py")
spec = importlib.util.spec_from_file_location("cc_soft_champion", SCRIPT)
champion = importlib.util.module_from_spec(spec)
spec.loader.exec_module(champion)


class ChampionProductionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.features = self.root / "features_full145_allowed_name.tsv"
        self.km = self.root / "km.tsv"
        self.gmm = self.root / "gmm.tsv"
        self.out = self.root / "predictions.tsv"

    def table(self, path, columns, rows):
        with open(path, "w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(rows)
        return path

    def components(self, keys, km_values, gmm_values):
        self.table(self.features, ("file", "lobe"), keys)
        for path, values in ((self.km, km_values), (self.gmm, gmm_values)):
            self.table(path, ("file", "lobe", "predicted", "probability_1"),
                       [(*key, "1", value) for key, value in zip(keys, values)])

    def vote(self):
        champion.write_soft_vote(self.features, self.km, self.gmm, self.out)
        with open(self.out, newline="") as stream:
            return list(csv.DictReader(stream, delimiter="\t"))

    def test_valid_vote_matches_frozen_arithmetic_including_tie(self):
        keys = [("unknown.sxm", i) for i in range(1, 6)]
        a, b = [0, .1, .4, 1, .7], [.2, .8, .6, .9, 1]
        self.components(keys, a, b)
        rows = self.vote()
        for row, ka, gb in zip(rows, a, b):
            p = (ka + gb) / 2
            self.assertEqual(row["predicted"], str(int(p >= .5)))
            self.assertEqual(row["confidence"], f"{abs(p-.5)*2:.8f}")
            self.assertEqual(row["probability_1"], f"{p:.8f}")
            self.assertEqual(row["invalid_reason"], "ok")

    def test_no_named_file_exclusion(self):
        self.components([("240310_Cu100009.sxm", 1)], [0], [0])
        self.assertEqual(self.vote()[0]["file"], "240310_Cu100009.sxm")

    def test_missing_component_keys_are_not_dropped(self):
        self.components([("long.sxm", 1), ("long.sxm", 2)], [.1, .8], [.2])
        rows = self.vote()
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[1]["predicted"], "?")
        self.assertEqual(rows[1]["confidence"], "0.00000000")
        self.assertEqual(rows[1]["invalid_reason"], "missing_component:gmm")

    def test_unavailable_component_probability_abstains(self):
        for unavailable in ("NA", "NaN", "Inf", ""):
            with self.subTest(unavailable=unavailable):
                self.components([("long.sxm", 1)], [unavailable], [.8])
                self.assertEqual(self.vote()[0]["predicted"], "?")
                self.out.unlink()

    def test_explicit_component_abstention_is_preserved(self):
        self.components([("long.sxm", 1)], [.5], [.8])
        self.table(self.km, ("file", "lobe", "predicted", "probability_1"),
                   [("long.sxm", 1, "?", .5)])
        self.assertEqual(self.vote()[0]["predicted"], "?")

    def test_invalid_probability_fails_before_output(self):
        self.components([("long.sxm", 1)], [1.2], [.8])
        with self.assertRaisesRegex(ValueError, "outside"):
            self.vote()
        self.assertFalse(self.out.exists())

    def test_extra_component_key_is_rejected(self):
        self.components([("long.sxm", 1)], [.1], [.2])
        self.table(self.gmm, ("file", "lobe", "predicted", "probability_1"),
                   [("other.sxm", 1, "0", .2)])
        with self.assertRaisesRegex(ValueError, "absent from features"):
            self.vote()

    def test_duplicate_keys_are_rejected(self):
        self.table(self.features, ("file", "lobe"), [("long.sxm", 1)] * 2)
        with self.assertRaisesRegex(ValueError, "Duplicate lobe"):
            champion.read_lobe_table(self.features)

    def test_control_columns_are_rejected(self):
        for column in ("sequence", "expected_N", "target_N", "control_sequence",
                       "benchmark_label", "truth"):
            with self.subTest(column=column):
                self.table(self.features, ("file", "lobe", column), [("long.sxm", 1, "x")])
                with self.assertRaisesRegex(ValueError, "Forbidden"):
                    champion.read_lobe_table(self.features)

    def test_existing_output_is_not_overwritten(self):
        self.components([("long.sxm", 1)], [.1], [.2])
        self.out.write_text("existing user output")
        with self.assertRaises(FileExistsError):
            self.vote()
        self.assertEqual(self.out.read_text(), "existing user output")

    def test_benchmark_flag_cannot_start_work(self):
        with patch.object(champion, "run") as runner, patch("sys.stderr", new_callable=io.StringIO):
            with self.assertRaises(SystemExit) as error:
                champion.main(["--truth", "forbidden.tsv"])
        self.assertNotEqual(error.exception.code, 0)
        runner.assert_not_called()

    def test_main_only_invokes_frozen_builders_and_chosen_runtimes(self):
        self.table(self.features, ("file", "lobe", "skew_ratio"), [("long.sxm", 1, 2)])
        patches = self.table(self.root / "patches.tsv", ("file", "lobe"), [("long.sxm", 1)])
        cube = self.root / "cube"
        cube.write_text("stub")
        work = self.root / "work"
        commands = []

        def stub_run(cmd):
            commands.append(cmd)
            script = str(cmd[2]) if cmd[1] == "--project=." else str(cmd[1])
            if script == champion.SCORER:
                self.table(Path(cmd[cmd.index("--out") + 1]),
                           ("file", "lobe", "cost_margin"), [("long.sxm", 1, .2)])
            elif script.endswith("empirical_fisher_mold.py"):
                self.table(Path(str(cmd[-1]) + "_cv.tsv"),
                           ("file", "lobe", "score"), [("long.sxm", 1, -.4)])
            elif script in (champion.GMM, champion.KM):
                self.table(Path(cmd[cmd.index("--out") + 1]),
                           ("file", "lobe", "predicted", "probability_1"),
                           [("long.sxm", 1, "0", .2)])

        args = ["--features", str(self.features), "--patches-fwd", str(patches),
                "--patches-bwd", str(patches), "--cube0", str(cube), "--cube1", str(cube),
                "--frame0", str(cube), "--frame1", str(cube), "--out", str(self.out),
                "--workdir", str(work), "--julia", "/chosen/julia-1.13"]
        with patch.object(champion, "run", side_effect=stub_run):
            champion.main(args)
        self.assertEqual(len(commands), 6)
        self.assertFalse(any("report_unit_assignment" in str(cmd) for cmd in commands))
        for cmd in commands:
            self.assertEqual(cmd[0], "/chosen/julia-1.13" if cmd[1] == "--project=."
                             else champion.sys.executable)
        self.assertTrue(self.out.is_file())


if __name__ == "__main__":
    unittest.main()
