# multiblock

`multiblock` provides tools to analyze and manipulate surface meshes using VTK multiblock structures.


multiblock-core/
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

# 4. Verify installation
cmake --version
conda list | egrep "vtk|openvdb|tbb|eigen"

# Verify VMTK installation after automated build
find $CONDA_PREFIX/lib -name "libvmtk*" && echo "✅ VMTK installed successfully"

# ✅ Environment 'cxxgeom' is now ready for building and automating VMTK
