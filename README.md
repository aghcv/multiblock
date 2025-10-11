# multiblock

`multiblock` provides tools to analyze and manipulate surface meshes using VTK multiblock structures.

```bash
multiblock/
├── CMakeLists.txt
├── cmake/                        # helper scripts (optional)
│   └── FindVMTK.cmake            # custom finder if not installed
├── src/
│   ├── main.cpp                  # example app
│   ├── centerline_extraction.cpp
│   ├── distance_field.cpp
│   ├── ownership_map.cpp
│   └── ...
├── include/
│   ├── fastvessels/
│   │   ├── centerline_extraction.hpp
│   │   ├── distance_field.hpp
│   │   ├── ownership_map.hpp
│   │   └── common.hpp
├── tests/
│   └── test_basic.cpp
└── README.md
```


# 🧩 Quickstart: Set up multiblock-core (C++ geometry environment)

# 1. Clone the repository
git clone https://github.com/<your-org>/<your-repo>.git
cd <your-repo>

# 2. Configure conda-forge channel (once)
conda config --add channels conda-forge
conda config --set channel_priority strict

# 3. Create a clean environment
conda env remove -n cxxgeom --yes || true
conda env create -f environment.yml
conda activate cxxgeom
bash .conda/post-link.sh

# 4. Verify installation
# already implemented in .coda/post_link.sh
# cmake --version
# conda list | egrep "vtk|openvdb|tbb|eigen"
# find $CONDA_PREFIX/lib -name "libvmtk*" && echo "✅ VMTK installed successfully"

# 5. Configure the CMake build
rm -rf build
mkdir -p build && cd build
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release

# 6. Build the project
ninja

# 7. Run example
./multiblock ../data/example.vtp
