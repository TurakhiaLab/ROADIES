#!/usr/bin/env python3

import subprocess
import argparse
import os
from pathlib import Path

ROADIES_ROOT = Path(__file__).resolve().parent

parser = argparse.ArgumentParser(description="Script to run ROADIES.")

parser.add_argument("--mode", default="accurate",
                    help="accurate | balanced | fast | placement")

parser.add_argument("--noconverge", action="store_true",
                    help="run in non-convergence mode")

parser.add_argument("--cores", type=int, default=32,
                    help="number of CPU cores")

parser.add_argument("--config", default=None,
                    help="config file path (default: config/config.yaml bundled with ROADIES)")

parser.add_argument("--deep", action="store_true",
                    help="enable deep phylogeny mode")

parser.add_argument("--gpu", type=int, default=0,
                    help="number of GPUs to use (0 = CPU mode)")

parser.add_argument("--grow", action="store_true",
                    help="Specify if you want to grow your tree or if you want to update your tree in placement mode")

parser.add_argument("--clean", action="store_true",
                    help="Delete the output directory before running, for a genuine fresh start "
                         "(default is to leave existing output alone - safer against accidental double-launches)")

parser.add_argument("--cluster", action="store_true",
                    help="Submit Snakemake rule jobs to SLURM via sbatch for multi-node execution, "
                         "instead of running everything locally on this machine")

args = parser.parse_args()

# Resolve the config path against the caller's cwd (not the repo root) before
# anything changes cwd, so a relative --config keeps meaning what the user typed.
if args.config is None:
    config_path = str(ROADIES_ROOT / "config" / "config.yaml")
else:
    config_path = os.path.abspath(args.config)

# Pick script
script = "noconverge.py" if args.noconverge else "converge.py"

# Build command
command = [
    "python", f"workflow/scripts/{script}",
    "--cores", str(args.cores),
    "--mode", args.mode,
    "--config", config_path,
    "--gpu", str(args.gpu),
]

if args.deep:
    command.append("--deep")

if args.grow:
    command.append("--grow")

if args.clean:
    command.append("--clean")

if args.cluster:
    command.append("--cluster")

print("Running:", " ".join(command))
# Run with cwd pinned to the ROADIES repo root so every relative path used
# internally by converge.py/noconverge.py/Snakemake resolves correctly,
# regardless of the directory this script was invoked from.
subprocess.run(command, check=True, cwd=str(ROADIES_ROOT))
