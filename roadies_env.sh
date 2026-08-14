#!/bin/bash
set -euo pipefail

# Required installations: (uncomment next 2 lines if you have sudo access, otherwise make sure following tools are installed before proceeding)
# sudo apt-get update
# sudo apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git vim screen default-jre libgomp1 libboost-all-dev cmake

# Define installation paths and check for directory
CONDA_PATH="${HOME}/conda"
ROADIES_ENV_SETUP="roadies_env.sh"

# Download and install Mambaforge if not already installed
if [ ! -d "${CONDA_PATH}" ]; then
    wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/download/24.11.3-2/Miniforge3-24.11.3-2-Linux-x86_64.sh"
    bash Miniforge3.sh -b -p "${CONDA_PATH}"
fi

# Create and setup the Conda environment if it doesn't exist
source "${CONDA_PATH}/etc/profile.d/conda.sh"  # Temporarily source conda for this script
source "${CONDA_PATH}/etc/profile.d/mamba.sh"  # Temporarily source mamba for this script

if ! conda env list | grep -q "roadies_env"; then
    mamba create -y -c conda-forge -c bioconda --name roadies_env snakemake alive-progress biopython iqtree=2.2.0.3 numpy lastz mashtree matplotlib seaborn treeswift=1.1.28 fasttree=2.1.11 python=3.11 raxml-ng ete3 lastz=1.04.52 aster=1.19 pyyaml seaborn epa-ng gappa libstdcxx-ng
fi
conda activate roadies_env

# Clone PASTA if not already done
if [ ! -d "pasta" ]; then
    git clone https://github.com/smirarab/pasta.git
fi

# Clone Sate-Tools if not already done
if [ ! -d "sate-tools-linux" ]; then
    git clone https://github.com/smirarab/sate-tools-linux.git
fi

# Setup PASTA if not already done
if [ -d "pasta" ]; then
    mafft_file="pasta/bin/mafft"
    if [ ! -f "$mafft_file" ]; then
        cd pasta
        python3 setup.py develop --user
        cd ..
    fi
fi

# Clone and build TWILIGHT (used for placement-mode MSA) from source, same as PASTA above
if [ ! -d "TWILIGHT" ]; then
    git clone https://github.com/TurakhiaLab/TWILIGHT.git
fi
if [ -d "TWILIGHT" ] && [ ! -f "TWILIGHT/bin/twilight" ]; then
    # TWILIGHT/CMakeLists.txt hardcodes -march=native, which ties the binary
    # to this specific host's CPU features and can also break nvcc parsing
    # newer GCC's AVX512BF16/AMX intrinsics headers. Portable baseline
    # instead (AVX2/FMA/BMI2, safe since Haswell/2013).
    sed -i 's/-march=native/-march=x86-64-v3/g' TWILIGHT/CMakeLists.txt
    cd TWILIGHT
    # Force TWILIGHT's vendored oneTBB build rather than trusting whatever
    # libtbb-dev a host's package manager provides (some don't ship a modern
    # CMake config), and route around two oneTBB build issues on newer
    # toolchains: cmake >=4 refusing its cmake_minimum_required(<3.5), and
    # its test suite hitting a GCC 12.4 -Wstringop-overflow false positive.
    CMAKE_WRAP_DIR="$(mktemp -d)"
    REAL_CMAKE="$(command -v cmake)"
    # CMake's check_language(CUDA) auto-detection doesn't reliably find nvcc
    # via PATH/CUDACXX in every environment - pass it explicitly when found.
    REAL_NVCC="$(command -v nvcc || true)"
    cat > "${CMAKE_WRAP_DIR}/cmake" <<EOF
#!/bin/bash
if [[ "\$1" == --* ]]; then
    exec "${REAL_CMAKE}" "\$@"
else
    exec "${REAL_CMAKE}" -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DTBB_TEST=OFF ${REAL_NVCC:+-DCMAKE_CUDA_COMPILER=${REAL_NVCC}} "\$@"
fi
EOF
    chmod +x "${CMAKE_WRAP_DIR}/cmake"
    cat > "${CMAKE_WRAP_DIR}/dpkg" <<'EOF'
#!/bin/bash
exit 1
EOF
    chmod +x "${CMAKE_WRAP_DIR}/dpkg"
    if command -v nvcc &>/dev/null; then
        PATH="${CMAKE_WRAP_DIR}:${PATH}" bash install/buildTWILIGHT.sh cuda
    else
        PATH="${CMAKE_WRAP_DIR}:${PATH}" bash install/buildTWILIGHT.sh
    fi
    rm -rf "${CMAKE_WRAP_DIR}"
    cd ..
fi

# Build MLIPPER (GPU placement tool) from source, best-effort: only needed for
# `--mode placement --gpu`, so skip without failing the rest of the setup if
# CUDA or libpll aren't available on this system.
if command -v nvcc &>/dev/null && { [ -f /usr/local/include/libpll/pll.h ] || [ -f /usr/include/libpll/pll.h ]; }; then
    echo "Building MLIPPER from source..."
    (cd MLIPPER && make) || echo "Warning: MLIPPER build failed; keeping the existing MLIPPER/MLIPPER binary, if any."
else
    echo "Warning: CUDA (nvcc) and/or libpll not found; skipping MLIPPER build. GPU-accelerated placement (--mode placement with --gpu) needs MLIPPER - once CUDA/libpll are available, build it with: bash MLIPPER/install/setup_host.sh"
fi

# Build sampling code
if [ ! -d "workflow/scripts/sampling/build" ]; then
    cd workflow/scripts/sampling
    mkdir build
    cd build
    cmake ..
    make
    cd ../../../..
fi

# Install ete3 for the user, only if it is not already installed
python3 -m pip show ete3 &>/dev/null || python3 -m pip install --user ete3

echo "Setup complete. Remember to source '${ROADIES_ENV_SETUP}' before running your projects."
