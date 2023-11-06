import sys
import subprocess
import argparse

if len(sys.argv) < 2 or "--help" in sys.argv:
    print(
        "Usage: run_roadies.py --cores <num_cores> --mode <mode_option> --converge --config <path of config file>"
    )
    sys.exit(1)

# Default values
mode_option = "accurate"
converge = False
num_cores = 32
config_path = "config/config.yaml"

# Create the parser
parser = argparse.ArgumentParser(description='Script to run ROADIES.')

# Add arguments
parser.add_argument('--mode', help='specify mode of operation - <accurate>, <balanced> OR <fast>')
parser.add_argument('--converge', action='store_true', help='specify ROADIES to run in convergence mode')
parser.add_argument('--cores', type=int, help='specify number of cores')
parser.add_argument('--config', help='specify path of config file')

# Parse the arguments
args = parser.parse_args()

# Set values
if args.mode:
    mode_option = args.mode
if args.converge:
    converge = True
if args.cores:
    num_cores = args.cores
if args.config:
    config_path = args.config

# Define the scripts to be executed based on the arguments
if converge:
    if mode_option == "accurate" or mode_option == "balanced" or mode_option == "fast":
        subprocess.run(
            [
                "python",
                "workflow/scripts/converge.py",
                "--cores",
                str(num_cores),
                "--mode",
                mode_option,
                "--config",
                config_path,
            ]
        )
    else:
        print("Invalid mode option.")
else:
    if mode_option == "accurate" or mode_option == "balanced" or mode_option == "fast":
        subprocess.run(
            [
                "python",
                "workflow/scripts/noconverge.py",
                "--cores",
                str(num_cores),
                "--mode",
                mode_option,
                "--config",
                config_path,
            ]
        )
    else:
        print("Invalid mode option.")
