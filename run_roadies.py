#!/usr/bin/python3

import subprocess
import argparse

parser = argparse.ArgumentParser(description='Script to run ROADIES.')

parser.add_argument('--mode', default='accurate', help='specify mode of operation - <accurate>, <balanced> OR <fast>')
parser.add_argument('--converge', action='store_true', help='specify ROADIES to run in convergence mode')
parser.add_argument('--cores', type=int, default=32, help='specify number of cores')
parser.add_argument('--config', default='config/config.yaml', help='specify path of config file')

args = parser.parse_args()

script = "converge.py" if args.converge else "noconverge.py"

subprocess.run(
    [
        "python",
        f"workflow/scripts/{script}",
        "--cores",
        str(args.cores),
        "--mode",
        args.mode,
        "--config",
        args.config,
    ],
    check=True
)