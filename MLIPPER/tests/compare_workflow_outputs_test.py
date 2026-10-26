#!/usr/bin/env python3

import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compare_workflow_outputs import compare_jplace, compare_newick


class WorkflowOutputComparisonTest(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.directory = Path(self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def write_json(self, name, value):
        path = self.directory / name
        path.write_text(json.dumps(value))
        return path

    def write_text(self, name, value):
        path = self.directory / name
        path.write_text(value)
        return path

    def test_jplace_ignores_invocation_and_tolerates_small_float_drift(self):
        expected = {
            "tree": "(A:0.1{0},B:0.1{1});",
            "placements": [{"n": ["Q"], "p": [[0, -10.0, 1.0, 0.1, 0.2]]}],
            "metadata": {"invocation": "old"},
            "version": 3,
            "fields": ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"],
        }
        actual = copy.deepcopy(expected)
        actual["metadata"]["invocation"] = "new"
        actual["placements"][0]["p"][0][1] += 1.0e-7
        compare_jplace(
            self.write_json("expected.jplace", expected),
            self.write_json("actual.jplace", actual),
        )

    def test_jplace_rejects_changed_edge_assignment(self):
        expected = {
            "placements": [{"n": ["Q"], "p": [[0, -10.0]]}],
            "metadata": {},
        }
        actual = copy.deepcopy(expected)
        actual["placements"][0]["p"][0][0] = 1
        with self.assertRaises(AssertionError):
            compare_jplace(
                self.write_json("expected.jplace", expected),
                self.write_json("actual.jplace", actual),
            )

    def test_newick_accepts_reordering_but_rejects_topology_change(self):
        expected = self.write_text(
            "expected.nwk", "((A:0.1,B:0.2):0.3,C:0.4,D:0.5);"
        )
        reordered = self.write_text(
            "reordered.nwk", "(D:0.5,C:0.4,(B:0.2,A:0.1):0.3);"
        )
        changed = self.write_text(
            "changed.nwk", "((A:0.1,C:0.4):0.3,B:0.2,D:0.5);"
        )
        compare_newick(expected, reordered)
        with self.assertRaises(AssertionError):
            compare_newick(expected, changed)

    def test_newick_rejects_large_branch_length_change(self):
        expected = self.write_text("expected.nwk", "(A:0.1,B:0.2,C:0.3);")
        actual = self.write_text("actual.nwk", "(A:0.1,B:0.25,C:0.3);")
        with self.assertRaises(AssertionError):
            compare_newick(expected, actual)


if __name__ == "__main__":
    unittest.main()
