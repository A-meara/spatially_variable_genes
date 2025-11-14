#!/bin/bash
# Script to build all Apptainer/Singularity images for the spatially variable genes benchmark
# This builds images in the correct order (base images first, then method-specific images)

set -e  # Exit on any error

echo "========================================"
echo "Building Apptainer Images"
echo "========================================"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Step 1: Build base images
echo "Step 1/3: Building base images..."
echo "-----------------------------------"

echo "Building python_base_ob.sif..."
apptainer build python_base_ob.sif python_base_ob.def
echo "✓ python_base_ob.sif built successfully"
echo ""

echo "Building r_base_ob.sif..."
apptainer build r_base_ob.sif r_base_ob.def
echo "✓ r_base_ob.sif built successfully"
echo ""

# Step 2: Build R method-specific images
echo "Step 2/2: Building R method images..."
echo "-----------------------------------"

echo "Building boostgp.sif..."
apptainer build boostgp.sif boostgp.def
echo "✓ boostgp.sif built successfully"
echo ""

echo "Building spark_x.sif..."
apptainer build spark_x.sif spark_x.def
echo "✓ spark_x.sif built successfully"
echo ""

echo "========================================"
echo "Build Complete!"
echo "========================================"
echo ""
echo "Built images:"
ls -lh *.sif
echo ""
echo "Note: random_ranking, true_ranking, and correlation use python_base_ob.sif directly (no separate build needed)"
