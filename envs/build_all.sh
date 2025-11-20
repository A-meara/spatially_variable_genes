#!/bin/bash
# Script to build all Apptainer/Singularity images for the spatially variable genes benchmark
# This builds images in the correct order (base images first, then method-specific images)
#
# Usage:
#   ./build_all.sh               # Interactive mode - prompts before rebuilding existing images
#   ./build_all.sh --yes         # Auto-rebuild mode - rebuilds all images without prompting
#   ./build_all.sh --new-only    # Build only images that don't exist yet (skip existing)

set -e  # Exit on any error

# Parse command-line arguments
AUTO_YES=false
NEW_ONLY=false
for arg in "$@"; do
    if [[ "$arg" == "--yes" || "$arg" == "-y" ]]; then
        AUTO_YES=true
    elif [[ "$arg" == "--new-only" || "$arg" == "-n" ]]; then
        NEW_ONLY=true
    fi
done

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

# Function to build a base image
build_base_image() {
    local base_name=$1
    local def_file="${base_name}.def"
    local sif_file="${base_name}.sif"

    if [ -f "$sif_file" ]; then
        echo "$sif_file already exists."
        if [ "$NEW_ONLY" = true ]; then
            echo "⊙ Skipping $sif_file (--new-only mode)"
        elif [ "$AUTO_YES" = true ]; then
            echo "Auto-rebuild mode: rebuilding $sif_file..."
            apptainer build --force "$sif_file" "$def_file"
            echo "✓ $sif_file rebuilt successfully"
        else
            read -p "Do you want to rebuild it? (y/N): " rebuild_response
            rebuild_response=${rebuild_response:-N}
            if [[ $rebuild_response =~ ^[Yy]$ ]]; then
                echo "Rebuilding $sif_file..."
                apptainer build --force "$sif_file" "$def_file"
                echo "✓ $sif_file rebuilt successfully"
            else
                echo "⊙ Skipping $sif_file (using existing)"
            fi
        fi
    else
        echo "Building $sif_file..."
        apptainer build "$sif_file" "$def_file"
        echo "✓ $sif_file built successfully"
    fi
    echo ""
}

# Build base images in order
build_base_image "python_base_ob"
build_base_image "r_base_ob"
build_base_image "base_pytorch_nvidia"
build_base_image "base_tensorflow_nvidia"

# Step 2: Build method-specific images (anything that depends on base images)
echo "Step 2/2: Building method-specific images..."
echo "-----------------------------------"

# Find all .def files that are NOT base images
method_defs=$(find . -maxdepth 1 -name "*.def" ! -name "python_base_ob.def" ! -name "r_base_ob.def" ! -name "base_pytorch_nvidia.def" ! -name "base_tensorflow_nvidia.def" -exec basename {} \; | sort)

if [ -z "$method_defs" ]; then
    echo "No method-specific .def files found"
else
    for def_file in $method_defs; do
        sif_file="${def_file%.def}.sif"

        # Check if sif file already exists
        if [ -f "$sif_file" ]; then
            if [ "$NEW_ONLY" = true ]; then
                echo "⊙ Skipping $sif_file (already exists, --new-only mode)"
                continue
            elif [ "$AUTO_YES" = true ]; then
                echo "Building $sif_file from $def_file..."
                apptainer build --force "$sif_file" "$def_file"
                echo "✓ $sif_file built successfully"
                echo ""
            else
                # Interactive mode - file exists, skip it
                echo "⊙ Skipping $sif_file (already exists)"
                continue
            fi
        else
            # File doesn't exist, build it
            echo "Building $sif_file from $def_file..."
            apptainer build "$sif_file" "$def_file"
            echo "✓ $sif_file built successfully"
            echo ""
        fi
    done
fi

echo "========================================"
echo "Build Complete!"
echo "========================================"
echo ""
echo "Built images:"
ls -lh *.sif 2>/dev/null || echo "No .sif files found"
echo ""
echo "Note: Some methods may use base images directly without separate .def files"
