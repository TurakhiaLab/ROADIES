import sys
import subprocess

# Check for the arguments
if len(sys.argv) < 2:
    print(
        "Usage: run_roadies.py --cores <num_cores> --mode <mode_option> --converge --cores <num_cores>"
    )
    sys.exit(1)

# Default values
mode_option = "accurate"
converge = False
num_cores = 32

# Parse command line arguments
for i in range(1, len(sys.argv)):
    if sys.argv[i] == "--mode" and i + 1 < len(sys.argv):
        mode_option = sys.argv[i + 1]
    elif sys.argv[i] == "--converge":
        converge = True
    elif sys.argv[i] == "--cores" and i + 1 < len(sys.argv):
        num_cores = int(sys.argv[i + 1])

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
            ]
        )
    else:
        print("Invalid mode option.")
