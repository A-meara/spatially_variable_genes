#!/usr/bin/env python3
"""
scGCO wrapper for omnibenchmark.

This script runs in Python 3.12 (base omnibenchmark environment) and calls
the actual scGCO implementation inside a Python 3.9 Singularity container.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name
  --data.dataset: Input dataset h5ad file
"""

import subprocess
import argparse
import os
import sys

# Parse command line arguments
parser = argparse.ArgumentParser(description='scGCO wrapper for spatially variable gene detection')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Output directory where files will be saved')
parser.add_argument('--name', type=str, required=True,
                    help='Name of the dataset')
parser.add_argument('--data.dataset', type=str, required=True,
                    help='Path to input dataset h5ad file')

args = parser.parse_args()

# Path to the scgco singularity image
# This assumes the scgco.sif is available at a known location
SCGCO_SIF = os.environ.get('SCGCO_SIF', '/path/to/scgco.sif')

if not os.path.exists(SCGCO_SIF):
    print(f"ERROR: scGCO Singularity image not found at {SCGCO_SIF}", file=sys.stderr)
    print(f"Set SCGCO_SIF environment variable to the correct path", file=sys.stderr)
    sys.exit(1)

# Path to the actual scGCO script inside the container
SCGCO_SCRIPT = '/opt/scGCO_script.py'

print(f"Running scGCO via Singularity container: {SCGCO_SIF}", flush=True)
print(f"Input: {getattr(args, 'data.dataset')}", flush=True)
print(f"Output: {args.output_dir}/{args.name}.predictions.h5ad", flush=True)

# Build the singularity exec command
cmd = [
    'singularity', 'exec',
    '--bind', f'{os.path.dirname(getattr(args, "data.dataset"))}:/input:ro',  # Bind input dir as read-only
    '--bind', f'{args.output_dir}:/output:rw',  # Bind output dir as read-write
    SCGCO_SIF,
    'python3', SCGCO_SCRIPT,
    '--input_data', getattr(args, 'data.dataset'),
    '--output', os.path.join(args.output_dir, f"{args.name}.predictions.h5ad"),
    '--name', args.name
]

# Execute the command
try:
    result = subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True
    )
    print(result.stdout, flush=True)
    if result.stderr:
        print(result.stderr, file=sys.stderr, flush=True)

except subprocess.CalledProcessError as e:
    print(f"ERROR: scGCO failed with exit code {e.returncode}", file=sys.stderr)
    print(f"STDOUT:\n{e.stdout}", file=sys.stderr)
    print(f"STDERR:\n{e.stderr}", file=sys.stderr)
    sys.exit(e.returncode)

print("scGCO wrapper completed successfully!", flush=True)
