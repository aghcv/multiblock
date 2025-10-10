#!/usr/bin/env bash
set -e
echo "🔧 Building VMTK from source..."

git clone --depth=1 https://github.com/vmtk/vmtk.git
cd vmtk
mkdir build && cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
  -DBUILD_SHARED_LIBS=ON
cmake --build . --parallel $(nproc)
cmake --install .

echo "✅ VMTK installed to $CONDA_PREFIX"
