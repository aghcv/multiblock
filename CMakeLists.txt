cmake_minimum_required(VERSION 3.20)
project(multiblock LANGUAGES CXX)

# Add local cmake/ directory for custom find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Ensure CMake searches the conda environment's cmake path
if(DEFINED ENV{CONDA_PREFIX})
    list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}/lib/cmake")
endif()

# ------------------------------------------------------------
# Build configuration
# ------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------
find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    FiltersCore
    FiltersModeling
    IOXML
    IOLegacy
)
# Search in module mode since FindOpenVDB.cmake is installed under $CONDA_PREFIX/lib/cmake/OpenVDB
if(DEFINED ENV{CONDA_PREFIX})
    list(APPEND CMAKE_MODULE_PATH
        "$ENV{CONDA_PREFIX}/lib/cmake/OpenVDB"
        "$ENV{CONDA_PREFIX}/lib/cmake/openvdb"
    )
endif()

find_package(OpenVDB MODULE REQUIRED)
find_package(TBB REQUIRED)
find_package(Eigen3 REQUIRED)

# Optional VMTK (quietly)
find_package(VMTK QUIET)
if(VMTK_FOUND)
    message(STATUS "✅ Found VMTK at ${VMTK_DIR}")
    add_definitions(-DWITH_VMTK)
else()
    message(WARNING "⚠️ VMTK not found — continuing without it.")
endif()

# ------------------------------------------------------------
# Source files
# ------------------------------------------------------------
file(GLOB SRC_FILES src/*.cpp)

# ------------------------------------------------------------
# Executable target
# ------------------------------------------------------------
add_executable(multiblock ${SRC_FILES})

# Make headers visible to all targets
target_include_directories(multiblock
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Link external libraries
target_link_libraries(multiblock
    PRIVATE
        VTK::CommonCore
        VTK::CommonDataModel
        VTK::FiltersCore
        VTK::FiltersModeling
        VTK::IOXML
        VTK::IOLegacy
        OpenVDB::openvdb
        TBB::tbb
        Eigen3::Eigen
)

if(VMTK_FOUND)
    target_link_libraries(multiblock PRIVATE
        VMTK::vmtkCommon
        VMTK::vmtkComputationalGeometry
        VMTK::vmtkIO
    )
endif()

# Enable automatic VTK module initialization
vtk_module_autoinit(
    TARGETS multiblock
    MODULES
        VTK::CommonCore
        VTK::CommonDataModel
        VTK::FiltersCore
        VTK::FiltersModeling
        VTK::IOXML
        VTK::IOLegacy
)
