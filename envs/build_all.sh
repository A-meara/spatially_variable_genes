#!/bin/bash
# Script to build all Apptainer/Singularity images for the spatially variable genes benchmark
# This builds images in the correct order (base images first, then method-specific images)
#
# Usage:
#   ./build_all.sh        # Interactive mode - prompts before rebuilding existing images
#   ./build_all.sh --yes  # Auto-rebuild mode - rebuilds all images without prompting

set -e  # Exit on any error

# Parse command-line arguments
AUTO_YES=false
if [[ "$1" == "--yes" || "$1" == "-y" ]]; then
    AUTO_YES=true
fi

echo "========================================"
echo "Building Apptainer Images"
echo "========================================"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Step 1: Build base images
echo "Step 1/2: Building base images..."
echo "-----------------------------------"

# Check if python_base_ob.sif exists
if [ -f "python_base_ob.sif" ]; then
    echo "python_base_ob.sif already exists."
    if [ "$AUTO_YES" = true ]; then
        echo "Auto-rebuild mode: rebuilding python_base_ob.sif..."
        apptainer build --force python_base_ob.sif python_base_ob.def
        echo "✓ python_base_ob.sif rebuilt successfully"
    else
        read -p "Do you want to rebuild it? (y/N): " rebuild_python
        rebuild_python=${rebuild_python:-N}
        if [[ $rebuild_python =~ ^[Yy]$ ]]; then
            echo "Rebuilding python_base_ob.sif..."
            apptainer build --force python_base_ob.sif python_base_ob.def
            echo "✓ python_base_ob.sif rebuilt successfully"
        else
            echo "⊙ Skipping python_base_ob.sif (using existing)"
        fi
    fi
else
    echo "Building python_base_ob.sif..."
    apptainer build python_base_ob.sif python_base_ob.def
    echo "✓ python_base_ob.sif built successfully"
fi
echo ""

# Check if r_base_ob.sif exists
if [ -f "r_base_ob.sif" ]; then
    echo "r_base_ob.sif already exists."
    if [ "$AUTO_YES" = true ]; then
        echo "Auto-rebuild mode: rebuilding r_base_ob.sif..."
        apptainer build --force r_base_ob.sif r_base_ob.def
        echo "✓ r_base_ob.sif rebuilt successfully"
    else
        read -p "Do you want to rebuild it? (y/N): " rebuild_r
        rebuild_r=${rebuild_r:-N}
        if [[ $rebuild_r =~ ^[Yy]$ ]]; then
            echo "Rebuilding r_base_ob.sif..."
            apptainer build --force r_base_ob.sif r_base_ob.def
            echo "✓ r_base_ob.sif rebuilt successfully"
        else
            echo "⊙ Skipping r_base_ob.sif (using existing)"
        fi
    fi
else
    echo "Building r_base_ob.sif..."
    apptainer build r_base_ob.sif r_base_ob.def
    echo "✓ r_base_ob.sif built successfully"
fi
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
