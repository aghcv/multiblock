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

---

# 🧩 Quickstart: Set up multiblock-core (C++ geometry environment)

## 1. Clone the repository
```bash
git clone https://github.com/<your-org>/<your-repo>.git
cd <your-repo>
```

## 2. Configure conda-forge channel (once)
```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```

## 3. Create a clean environment
```bash
conda env remove -n cxxgeom --yes || true
conda env create -f environment.yml
conda activate cxxgeom
bash .conda/post-link.sh
```

## 4. Verify installation
```bash
# already implemented in .conda/post_link.sh
# cmake --version
# conda list | egrep "vtk|openvdb|tbb|eigen"
# find $CONDA_PREFIX/lib -name "libvmtk*" && echo "✅ VMTK installed successfully"
```

## 5. Configure the CMake build
```bash
rm -rf build
mkdir -p build && cd build
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release
```

## 6. Build the project
```bash
ninja
```

## 7. Run example
```bash
./multiblock ../data/example.vtp
```

---

# 🧠 Developer Setup (Cross-IDE)

To enable smooth development and header auto-completion across all IDEs (VS Code, CLion, Qt Creator, Visual Studio, etc.),  
`multiblock` exports a `compile_commands.json` file that describes all compiler flags and include paths.

## 1️⃣ CMake configuration
In the top-level `CMakeLists.txt`, make sure this line is present:
```cmake
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
```

This automatically generates a `build/compile_commands.json` file after configuration.

## 2️⃣ Clone + build
```bash
git clone https://github.com/<your-org>/<your-repo>.git
cd <your-repo>
mkdir build && cd build
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release
ninja
```

## 3️⃣ IDE support
Most IDEs can detect `compile_commands.json` automatically.  
If yours doesn’t, open your IDE’s **CMake** or **Project Settings** and point it to:

```
/path/to/multiblock/build/compile_commands.json
```

## 4️⃣ Optional (VS Code only)
If you use VS Code, create a file `.vscode/settings.json` at the repository root:

```json
{
    "C_Cpp.default.compileCommands": "${workspaceFolder}/build/compile_commands.json"
}
```

Reload the VS Code window (`Cmd ⇧ P → Developer: Reload Window`) and all VTK headers should resolve automatically 🎯
