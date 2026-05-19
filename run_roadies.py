#!/usr/bin/env python3

import subprocess
import argparse

parser = argparse.ArgumentParser(description="Script to run ROADIES.")

parser.add_argument("--mode", default="accurate",
                    help="accurate | balanced | fast | placement")

parser.add_argument("--noconverge", action="store_true",
                    help="run in non-convergence mode")

parser.add_argument("--cores", type=int, default=32,
                    help="number of CPU cores")

parser.add_argument("--config", default="config/config.yaml",
                    help="config file path")

parser.add_argument("--deep", action="store_true",
                    help="enable deep phylogeny mode")

parser.add_argument("--gpu", type=int, default=0,
                    help="number of GPUs to use (0 = CPU mode)")

parser.add_argument("--grow", action="store_true",
                    help="Specify if you want to grow your tree or if you want to update your tree in placement mode")

parser.add_argument("--no-clean", action="store_true",
                    help="Skip deleting the output directory before running (enables Snakemake resume)")

args = parser.parse_args()

# Pick script
script = "noconverge.py" if args.noconverge else "converge.py"

# Build command
command = [
    "python", f"workflow/scripts/{script}",
    "--cores", str(args.cores),
    "--mode", args.mode,
    "--config", args.config,
    "--gpu", str(args.gpu),
]

if args.deep:
    command.append("--deep")

if args.grow:
    command.append("--grow")

if args.no_clean:
    command.append("--no-clean")

print("Running:", " ".join(command))
subprocess.run(command, check=True)
