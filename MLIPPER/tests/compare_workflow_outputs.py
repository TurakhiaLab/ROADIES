#!/usr/bin/env python3

import argparse
import json
import math
from pathlib import Path


FLOAT_REL_TOLERANCE = 1.0e-7
FLOAT_ABS_TOLERANCE = 1.0e-6


def fail(message):
    raise AssertionError(message)


def compare_json(expected, actual, path="root"):
    if isinstance(expected, bool) or expected is None or isinstance(expected, str):
        if actual != expected:
            fail(f"{path}: expected {expected!r}, got {actual!r}")
        return
    if isinstance(expected, (int, float)):
        if not isinstance(actual, (int, float)):
            fail(f"{path}: expected a number, got {type(actual).__name__}")
        if isinstance(expected, int) and isinstance(actual, int):
            if actual != expected:
                fail(f"{path}: expected {expected}, got {actual}")
            return
        if not math.isclose(
            float(expected),
            float(actual),
            rel_tol=FLOAT_REL_TOLERANCE,
            abs_tol=FLOAT_ABS_TOLERANCE,
        ):
            fail(f"{path}: expected {expected:.17g}, got {actual:.17g}")
        return
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            fail(f"{path}: list length or type differs")
        for index, (expected_item, actual_item) in enumerate(zip(expected, actual)):
            compare_json(expected_item, actual_item, f"{path}[{index}]")
        return
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual) != set(expected):
            fail(f"{path}: object keys differ")
        for key in expected:
            compare_json(expected[key], actual[key], f"{path}.{key}")
        return
    fail(f"{path}: unsupported expected value type {type(expected).__name__}")


class NewickNode:
    def __init__(self, children=None, name="", length=None):
        self.children = children or []
        self.name = name
        self.length = length


class NewickParser:
    def __init__(self, text):
        self.text = text.strip()
        self.index = 0

    def parse(self):
        root = self.parse_node()
        self.skip_space()
        if self.peek() == ";":
            self.index += 1
        self.skip_space()
        if self.index != len(self.text):
            fail(f"Newick: unexpected text at offset {self.index}")
        return root

    def peek(self):
        return self.text[self.index] if self.index < len(self.text) else ""

    def skip_space(self):
        while self.peek().isspace():
            self.index += 1

    def parse_node(self):
        self.skip_space()
        children = []
        if self.peek() == "(":
            self.index += 1
            while True:
                children.append(self.parse_node())
                self.skip_space()
                if self.peek() == ",":
                    self.index += 1
                    continue
                if self.peek() != ")":
                    fail(f"Newick: expected ')' at offset {self.index}")
                self.index += 1
                break
        name = self.parse_token(":,();")
        length = None
        if self.peek() == ":":
            self.index += 1
            token = self.parse_token(",();")
            try:
                length = float(token)
            except ValueError:
                fail(f"Newick: invalid branch length {token!r}")
        return NewickNode(children, name, length)

    def parse_token(self, delimiters):
        self.skip_space()
        start = self.index
        while self.peek() and self.peek() not in delimiters:
            self.index += 1
        return self.text[start:self.index].strip()


def tree_edge_map(root):
    leaves = set()

    def collect(node):
        if not node.children:
            if not node.name:
                fail("Newick: unnamed leaf")
            if node.name in leaves:
                fail(f"Newick: duplicate leaf {node.name!r}")
            leaves.add(node.name)
            return {node.name}
        result = set()
        for child in node.children:
            result.update(collect(child))
        return result

    collect(root)
    edges = {}

    def visit(node):
        if not node.children:
            return {node.name}
        descendants = set()
        for child in node.children:
            child_leaves = visit(child)
            descendants.update(child_leaves)
            if child.length is None:
                fail("Newick: every non-root edge must have a branch length")
            other = leaves - child_leaves
            left = tuple(sorted(child_leaves))
            right = tuple(sorted(other))
            key = min(left, right)
            if key in edges:
                fail(f"Newick: duplicate split {key}")
            edges[key] = child.length
        return descendants

    visit(root)
    return leaves, edges


def compare_jplace(expected_path, actual_path):
    expected = json.loads(expected_path.read_text())
    actual = json.loads(actual_path.read_text())
    expected.get("metadata", {}).pop("invocation", None)
    actual.get("metadata", {}).pop("invocation", None)
    compare_json(expected, actual)


def compare_newick(expected_path, actual_path):
    expected_leaves, expected_edges = tree_edge_map(
        NewickParser(expected_path.read_text()).parse())
    actual_leaves, actual_edges = tree_edge_map(
        NewickParser(actual_path.read_text()).parse())
    if actual_leaves != expected_leaves:
        fail(
            "Newick taxa differ: "
            f"expected={sorted(expected_leaves)}, actual={sorted(actual_leaves)}"
        )
    if set(actual_edges) != set(expected_edges):
        fail("Newick topology differs: canonical split sets are not equal")
    for split, expected_length in expected_edges.items():
        actual_length = actual_edges[split]
        if not math.isclose(
            expected_length,
            actual_length,
            rel_tol=FLOAT_REL_TOLERANCE,
            abs_tol=FLOAT_ABS_TOLERANCE,
        ):
            fail(
                f"Newick branch {split}: expected {expected_length:.17g}, "
                f"got {actual_length:.17g}"
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected-jplace", type=Path, required=True)
    parser.add_argument("--actual-jplace", type=Path, required=True)
    parser.add_argument("--expected-tree", type=Path, required=True)
    parser.add_argument("--actual-tree", type=Path, required=True)
    args = parser.parse_args()

    compare_jplace(args.expected_jplace, args.actual_jplace)
    compare_newick(args.expected_tree, args.actual_tree)
    print("workflow_output_comparison: PASS")


if __name__ == "__main__":
    main()
