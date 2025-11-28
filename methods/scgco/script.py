#!/usr/bin/env python3
"""
scGCO wrapper for omnibenchmark - runs in Python 3.12.

This wrapper handles I/O in Python 3.12 (required by omnibenchmark),
but calls the actual scGCO implementation in a Python 3.9 Singularity container.

Omnibenchmark standard arguments:
  --output_dir: Directory where outputs will be saved
  --name: Dataset name
  --data.dataset: Input dataset h5ad file
"""

import subprocess
import argparse
import os
import sys

# Parse command line arguments (Python 3.12)
parser = argparse.ArgumentParser(description='scGCO wrapper for spatially variable gene detection')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Output directory where files will be saved')
parser.add_argument('--name', type=str, required=True,
                    help='Name of the dataset')
parser.add_argument('--data.dataset', type=str, required=True,
                    help='Path to input dataset h5ad file')

args = parser.parse_args()

# Get paths
input_data_path = getattr(args, 'data.dataset')
output_path = os.path.join(args.output_dir, f"{args.name}.predictions.h5ad")

# Create output directory
os.makedirs(args.output_dir, exist_ok=True)

# Path to scgco singularity image
# Check common locations
SCGCO_SIF_LOCATIONS = [
    os.environ.get('SCGCO_SIF'),
    '/home/ameara/porting/spatially_variable_genes/envs/scgco.sif',
    'envs/scgco.sif',
    '../envs/scgco.sif',
    '../../envs/scgco.sif',
]

SCGCO_SIF = None
for location in SCGCO_SIF_LOCATIONS:
    if location and os.path.exists(location):
        SCGCO_SIF = location
        break

if not SCGCO_SIF:
    print(f"ERROR: scGCO Singularity image not found!", file=sys.stderr)
    print(f"Searched locations:", file=sys.stderr)
    for loc in SCGCO_SIF_LOCATIONS:
        if loc:
            print(f"  - {loc}", file=sys.stderr)
    print(f"\nSet SCGCO_SIF environment variable to the correct path", file=sys.stderr)
    sys.exit(1)

print(f"Using scGCO container: {SCGCO_SIF}", flush=True)
print(f"Input: {input_data_path}", flush=True)
print(f"Output: {output_path}", flush=True)

# Get absolute paths for bind mounts
input_abs = os.path.abspath(input_data_path)
output_abs = os.path.abspath(output_path)
input_dir = os.path.dirname(input_abs)
output_dir = os.path.dirname(output_abs)

# Path to the implementation script inside the container
SCGCO_SCRIPT = '/opt/scgco_impl.py'

# Build singularity exec command
cmd = [
    'singularity', 'exec',
    '--bind', f'{input_dir}:{input_dir}:ro',  # Bind input dir as read-only
    '--bind', f'{output_dir}:{output_dir}:rw',  # Bind output dir as read-write
    SCGCO_SIF,
    'python3', SCGCO_SCRIPT,
    '--input_data', input_abs,
    '--output', output_abs,
    '--name', args.name
]

print(f"Running scGCO via Singularity...", flush=True)

# Execute the scgco container
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

    print("scGCO completed successfully!", flush=True)

except subprocess.CalledProcessError as e:
    print(f"ERROR: scGCO failed with exit code {e.returncode}", file=sys.stderr)
    if e.stdout:
        print(f"STDOUT:\n{e.stdout}", file=sys.stderr)
    if e.stderr:
        print(f"STDERR:\n{e.stderr}", file=sys.stderr)
    sys.exit(e.returncode)
except FileNotFoundError:
    print(f"ERROR: singularity command not found. Is Apptainer/Singularity installed?", file=sys.stderr)
    sys.exit(1)
