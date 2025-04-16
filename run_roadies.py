#!/usr/bin/env python3

import subprocess
import argparse

parser = argparse.ArgumentParser(description="Script to run ROADIES.")

parser.add_argument(
    "--mode",
    default="accurate",
    help="specify mode of operation - <accurate>, <balanced> OR <fast>",
)
parser.add_argument(
    "--noconverge", action="store_true", help="specify ROADIES to run in non convergence mode"
)
parser.add_argument("--cores", type=int, default=32, help="specify number of cores")
parser.add_argument(
    "--config", default="config/config.yaml", help="specify path of config file"
)
parser.add_argument(
    "--deep", default="False", help="specify if ROADIES will run in deep mode - to capture deeper phylogenetic timescales"
)

args = parser.parse_args()

script = "noconverge.py" if args.noconverge else "converge.py"

command = [
    "python",
    f"workflow/scripts/{script}",
    "--cores",
    str(args.cores),
    "--mode",
    args.mode,
    "--config",
    args.config,
    "--deep",
    args.deep,
]

subprocess.run(command, check=True)