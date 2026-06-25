# multiblock

Note: You are on the 'light' branch. This is a minimal VTK-only setup.

`multiblock` provides tools to analyze and manipulate surface meshes using VTK multiblock structures.

```bash
multiblock/
├── CMakeLists.txt
├── src/
│   ├── main.cpp                  # example app
│   ├── io_utils.cpp
│   └── ...
├── include/
│   ├── fastvessels/
│   │   ├── io_utils.hpp
│   │   └── ...
├── environment.yml
├── setup.sh
└── README.md
```

---

# 🧩 Quickstart: minimal VTK-only setup

## 1. Clone the repository
```bash
git clone https://github.com/<your-org>/<your-repo>.git
cd <your-repo>
```

## 2. One-step setup
Run the setup script:
```bash
bash setup.sh
```

HPC (ODU) example:
```bash
module load container_env python3

env USE_CRUN=1 ENV_NAME=multiblock ENV_PREFIX=~/envs bash setup.sh # EXAMPLE: env USE_CRUN=1 ENV_NAME=cxxgeom ENV_PREFIX=/home/aghorban/envs bash setup.sh
```

## 3. Run example
```bash
./build/multiblock.mbx
```
On first run, it generates a config file in output/ (default: output/solver_input.txt) and runs a small compute benchmark.

### Solver input
- Default config file: output/solver_input.txt (auto-generated if missing)
- Default geometry path: raw/geometry.vtp (VTP must exist before running the solver)
- Key fields:
    - `feature`: `default`, `cpu_accelerated`, `gpu_accelerated`, `cpu_gpu_accelerated`
    - `cpu_threads`: threads per rank (`0` = auto)
    - `num_gpus`: `-1` = auto-detect (env-based), `0` = CPU-only
    - `mode`: `run` or `preprocess`

### Preprocess (auto-benchmark and write optimized config)
Edit output/solver_input.txt:
```text
mode=preprocess
```
Then run:
```bash
./build/multiblock.mbx
```
It writes output/solver_input.txt (updated in-place by preprocess) and output/results.txt after a run.

---

# 🧠 Developer Setup (optional)

`multiblock` exports a `compile_commands.json` file for IDE integration.

## 1️⃣ CMake configuration
In the top-level `CMakeLists.txt`, make sure this line is present:
```cmake
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
```

This generates `build/compile_commands.json` after configuration.

## 2️⃣ Build
To rebuild manually:
```bash
mkdir -p build && cd build
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
