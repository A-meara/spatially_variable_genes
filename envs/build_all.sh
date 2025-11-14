#!/bin/bash
# Script to build all Apptainer/Singularity images for the spatially variable genes benchmark
# This builds images in the correct order (base images first, then method-specific images)
#
# Usage:
#   ./build_all.sh               # Interactive mode - prompts before rebuilding existing images
#   ./build_all.sh --yes         # Auto-rebuild mode - rebuilds all images without prompting
#   ./build_all.sh --yes --quiet # Auto-rebuild with minimal output and progress spinner

set -e  # Exit on any error

# Parse command-line arguments
AUTO_YES=false
QUIET=false
for arg in "$@"; do
    if [[ "$arg" == "--yes" || "$arg" == "-y" ]]; then
        AUTO_YES=true
    elif [[ "$arg" == "--quiet" || "$arg" == "-q" ]]; then
        QUIET=true
    fi
done

# Spinner function for quiet mode
spinner() {
    local pid=$1
    local delay=0.1
    local spinstr='|/-\'
    while ps -p $pid > /dev/null 2>&1; do
        local temp=${spinstr#?}
        printf " [%c]  " "$spinstr"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "\b\b\b\b\b\b"
    done
    printf "    \b\b\b\b"
}

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
        if [ "$QUIET" = true ]; then
            echo -n "Rebuilding python_base_ob.sif..."
            apptainer build --force --silent python_base_ob.sif python_base_ob.def > /dev/null 2>&1 &
            spinner $!
            wait $!
            echo " ✓"
        else
            echo "Auto-rebuild mode: rebuilding python_base_ob.sif..."
            apptainer build --force python_base_ob.sif python_base_ob.def
            echo "✓ python_base_ob.sif rebuilt successfully"
        fi
    else
        read -p "Do you want to rebuild it? (y/N): " rebuild_python
        rebuild_python=${rebuild_python:-N}
        if [[ $rebuild_python =~ ^[Yy]$ ]]; then
            if [ "$QUIET" = true ]; then
                echo -n "Rebuilding python_base_ob.sif..."
                apptainer build --force --silent python_base_ob.sif python_base_ob.def > /dev/null 2>&1 &
                spinner $!
                wait $!
                echo " ✓"
            else
                echo "Rebuilding python_base_ob.sif..."
                apptainer build --force python_base_ob.sif python_base_ob.def
                echo "✓ python_base_ob.sif rebuilt successfully"
            fi
        else
            echo "⊙ Skipping python_base_ob.sif (using existing)"
        fi
    fi
else
    if [ "$QUIET" = true ]; then
        echo -n "Building python_base_ob.sif..."
        apptainer build --silent python_base_ob.sif python_base_ob.def > /dev/null 2>&1 &
        spinner $!
        wait $!
        echo " ✓"
    else
        echo "Building python_base_ob.sif..."
        apptainer build python_base_ob.sif python_base_ob.def
        echo "✓ python_base_ob.sif built successfully"
    fi
fi
echo ""

# Check if r_base_ob.sif exists
if [ -f "r_base_ob.sif" ]; then
    echo "r_base_ob.sif already exists."
    if [ "$AUTO_YES" = true ]; then
        if [ "$QUIET" = true ]; then
            echo -n "Rebuilding r_base_ob.sif..."
            apptainer build --force --silent r_base_ob.sif r_base_ob.def > /dev/null 2>&1 &
            spinner $!
            wait $!
            echo " ✓"
        else
            echo "Auto-rebuild mode: rebuilding r_base_ob.sif..."
            apptainer build --force r_base_ob.sif r_base_ob.def
            echo "✓ r_base_ob.sif rebuilt successfully"
        fi
    else
        read -p "Do you want to rebuild it? (y/N): " rebuild_r
        rebuild_r=${rebuild_r:-N}
        if [[ $rebuild_r =~ ^[Yy]$ ]]; then
            if [ "$QUIET" = true ]; then
                echo -n "Rebuilding r_base_ob.sif..."
                apptainer build --force --silent r_base_ob.sif r_base_ob.def > /dev/null 2>&1 &
                spinner $!
                wait $!
                echo " ✓"
            else
                echo "Rebuilding r_base_ob.sif..."
                apptainer build --force r_base_ob.sif r_base_ob.def
                echo "✓ r_base_ob.sif rebuilt successfully"
            fi
        else
            echo "⊙ Skipping r_base_ob.sif (using existing)"
        fi
    fi
else
    if [ "$QUIET" = true ]; then
        echo -n "Building r_base_ob.sif..."
        apptainer build --silent r_base_ob.sif r_base_ob.def > /dev/null 2>&1 &
        spinner $!
        wait $!
        echo " ✓"
    else
        echo "Building r_base_ob.sif..."
        apptainer build r_base_ob.sif r_base_ob.def
        echo "✓ r_base_ob.sif built successfully"
    fi
fi
echo ""

# Step 2: Build method-specific images (anything that depends on base images)
echo "Step 2/2: Building method-specific images..."
echo "-----------------------------------"

# Find all .def files that are NOT base images
method_defs=$(find . -maxdepth 1 -name "*.def" ! -name "python_base_ob.def" ! -name "r_base_ob.def" -exec basename {} \; | sort)

if [ -z "$method_defs" ]; then
    echo "No method-specific .def files found"
else
    for def_file in $method_defs; do
        sif_file="${def_file%.def}.sif"
        if [ "$QUIET" = true ]; then
            echo -n "Building $sif_file..."
            apptainer build --force --silent "$sif_file" "$def_file" > /dev/null 2>&1 &
            spinner $!
            wait $!
            echo " ✓"
        else
            echo "Building $sif_file from $def_file..."
            apptainer build --force "$sif_file" "$def_file"
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
