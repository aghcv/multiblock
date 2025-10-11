#!/usr/bin/env bash
set -e
LOG="$CONDA_PREFIX/install_log.txt"

echo "🚀 Starting post-install builds..." | tee "$LOG"

# ============================================================
# 1️⃣ Build and install OpenVDB
# ============================================================
SRC_DIR="$CONDA_PREFIX/../openvdb-src"
if [ ! -d "$CONDA_PREFIX/include/openvdb" ]; then
  echo "🔧 Building OpenVDB from source..." | tee -a "$LOG"
  rm -rf "$SRC_DIR"
  cd $CONDA_PREFIX/..
  git clone --depth=1 https://github.com/AcademySoftwareFoundation/openvdb.git openvdb-src
  cd openvdb-src
  mkdir build && cd build

  cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
    -DUSE_BLOSC=OFF \
    -DUSE_EXR=OFF \
    -DUSE_TBB=ON \
    -DOPENVDB_BUILD_PYTHON_MODULE=OFF \
    -DOPENVDB_BUILD_BINARIES=OFF \
    -DOPENVDB_INSTALL_CMAKE_MODULES=ON \
    -DOPENVDB_USE_DEPRECATED_ABI=OFF \
    -DOPENVDB_BUILD_CORE=ON

  cmake --build . --parallel $(sysctl -n hw.ncpu) 2>&1 | tee -a "$LOG"

  # Explicitly install everything, including CMake exports
  cmake --install . 2>&1 | tee -a "$LOG"
  cmake --build . --target install_cmake 2>&1 | tee -a "$LOG"
  cmake --build . --target install_config 2>&1 | tee -a "$LOG"

  echo "✅ OpenVDB installed successfully" | tee -a "$LOG"

  cd $CONDA_PREFIX/..
  rm -rf "$SRC_DIR"
else
  echo "ℹ️ OpenVDB already present, skipping build." | tee -a "$LOG"
fi


# ============================================================
# 2️⃣ Build and install VMTK
# ============================================================
SRC_DIR="$CONDA_PREFIX/../vmtk-src"
if [ ! -d "$CONDA_PREFIX/include/vmtk" ]; then
  echo "🔧 Building VMTK from source..." | tee -a "$LOG"
  rm -rf "$SRC_DIR"
  cd $CONDA_PREFIX/..
  git clone --depth=1 https://github.com/vmtk/vmtk.git vmtk-src
  cd vmtk-src
  mkdir build && cd build

  cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
    -DBUILD_SHARED_LIBS=ON \
    -DVMTK_MINIMAL_BUILD=ON \
    -DVMTK_WITH_PYTHON=OFF

  cmake --build . --parallel $(sysctl -n hw.ncpu) 2>&1 | tee -a "$LOG"
  cmake --install . 2>&1 | tee -a "$LOG"

  echo "✅ VMTK installed successfully" | tee -a "$LOG"

  cd $CONDA_PREFIX/..
  rm -rf "$SRC_DIR"
else
  echo "ℹ️ VMTK already present, skipping build." | tee -a "$LOG"
fi



# ============================================================
# 3️⃣ Environment verification
# ============================================================
echo "🔎 Verifying environment..." | tee -a "$LOG"

cmake --version | tee -a "$LOG"
conda list | egrep "vtk|openvdb|tbb|eigen" | tee -a "$LOG"

# Library checks
find $CONDA_PREFIX/lib -name "libvmtk*" | tee -a "$LOG" && echo "✅ VMTK libraries found" | tee -a "$LOG"
find $CONDA_PREFIX/lib -name "libopenvdb*" | tee -a "$LOG" && echo "✅ OpenVDB libraries found" | tee -a "$LOG"

echo "🎉 Post-install setup complete."
