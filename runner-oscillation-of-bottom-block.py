#!/usr/bin/env python3
import os
from datetime import datetime
import sys

ARGS_NAMES = ['n', 'k', 'ang_freq', 'mu','n_osc', 'iters_per_osc']

# Generate timestamp
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
outputs_dir = f"outputs/{timestamp}"
source_file = "main-oscillation-of-bottom-block-copy.py"
input_file = "inputs_oscillation_of_bottom_block.txt"

# Check if the input file exists
if not os.path.isfile(input_file):
    print(f"Input file '{input_file}' not found.")
    exit(1)

# Create directory

# Read inputs from the file and run the script


def do_run(args, args_list):
    print(f"Running with args: {[(f'{ARGS_NAMES[i]}: {arg}') for i, arg in enumerate(args_list)]}")
    run_output_dir = f"{outputs_dir}/{'_'.join([(f'{ARGS_NAMES[i]}{arg}') for i, arg in enumerate(args_list)])}"
    os.makedirs(run_output_dir, exist_ok=True)
    os.system(
        f"python3 {source_file} -o {run_output_dir} {args} > {run_output_dir}/stdout.out")


if len(sys.argv) == 1:
    print("Reading from inputs.txt")
    with open(input_file, 'r') as f:
        for line in f:
            args = line.strip()
            args_arr = args.split()
            if (len(args_arr) != len(ARGS_NAMES)):
                raise Exception(f"inputs.txt input doesn't match expected args {ARGS_NAMES}")
            do_run(args, args_arr)
else:
    print("Using command line args")
    do_run(' '.join(sys.argv[1:]))
