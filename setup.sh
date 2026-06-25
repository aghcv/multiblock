#!/usr/bin/env bash
set -euo pipefail

# Minimal setup for the 'light' branch
# - Creates conda env with cmake, ninja, vtk
# - Configures and builds the project with Ninja

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
ENV_NAME="${ENV_NAME:-cxxgeom}"
ENV_PREFIX="${ENV_PREFIX:-$HOME/envs}"
USE_CRUN="${USE_CRUN:-}"

say() { echo "[setup] $*"; }

ENV_PATH="${ENV_PREFIX}/${ENV_NAME}"

if [[ -n "${USE_CRUN}" ]]; then
  if ! command -v crun >/dev/null 2>&1; then
    say "USE_CRUN is set but 'crun' is not available."
    say "Load modules first: module load container_env python3"
    exit 1
  fi
  if command -v module >/dev/null 2>&1; then
    module load container_env python3 || true
  else
    say "Module command not found; ensure modules are loaded before running."
  fi
else
  if ! command -v conda >/dev/null 2>&1; then
    say "Conda not found. Please install Miniforge/Conda and re-run."
    exit 1
  fi
fi

# 1) Conda channel configuration
if [[ -n "${USE_CRUN}" ]]; then
  say "Configuring conda-forge channel (strict priority) via crun"
  crun conda config --add channels conda-forge || true
  crun conda config --set channel_priority strict || true
else
  say "Configuring conda-forge channel (strict priority)"
  conda config --add channels conda-forge || true
  conda config --set channel_priority strict || true
fi

# 2) Create or update environment
say "Creating environment '${ENV_NAME}' from environment.yml"
if [[ -n "${USE_CRUN}" ]]; then
  mkdir -p "${ENV_PREFIX}"
  crun -p "${ENV_PATH}" conda env remove -p "${ENV_PATH}" --yes || true
  crun conda env create -p "${ENV_PATH}" -f "${PROJECT_DIR}/environment.yml"
else
  conda env remove -n "${ENV_NAME}" --yes || true
  conda env create -f "${PROJECT_DIR}/environment.yml" -n "${ENV_NAME}"
fi

# 3) Activate and build
say "Configuring CMake (Release)"
rm -rf "${PROJECT_DIR}/build"
if [[ -n "${USE_CRUN}" ]]; then
  crun -p "${ENV_PATH}" cmake -S "${PROJECT_DIR}" -B "${PROJECT_DIR}/build" -G Ninja -DCMAKE_BUILD_TYPE=Release
else
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${ENV_NAME}"
  cmake -S "${PROJECT_DIR}" -B "${PROJECT_DIR}/build" -G Ninja -DCMAKE_BUILD_TYPE=Release
fi

say "Building with Ninja"
if [[ -n "${USE_CRUN}" ]]; then
  crun -p "${ENV_PATH}" cmake --build "${PROJECT_DIR}/build" --config Release
else
  cmake --build "${PROJECT_DIR}/build" --config Release
fi

say "Done. Binary is at: ${PROJECT_DIR}/build/multiblock.mbx"
