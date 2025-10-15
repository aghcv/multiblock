cmake_minimum_required(VERSION 3.20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
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
# Core library (used by tests and executable)
# ------------------------------------------------------------
file(GLOB LIB_SOURCES CONFIGURE_DEPENDS
    src/io_utils.cpp
    src/centerline_extraction.cpp
    src/distance_field.cpp
    src/ownership_map.cpp
)

add_library(multiblock_lib STATIC ${LIB_SOURCES})
target_include_directories(multiblock_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(multiblock_lib
    PUBLIC
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
    target_link_libraries(multiblock_lib PUBLIC
        VMTK::vmtkCommon
        VMTK::vmtkComputationalGeometry
        VMTK::vmtkIO
    )
endif()

# ------------------------------------------------------------
# Executable (links the library)
# ------------------------------------------------------------
add_executable(multiblock src/main.cpp)
target_link_libraries(multiblock PRIVATE multiblock_lib)

vtk_module_autoinit(
    TARGETS multiblock multiblock_lib
    MODULES
        VTK::CommonCore
        VTK::CommonDataModel
        VTK::FiltersCore
        VTK::FiltersModeling
        VTK::IOXML
        VTK::IOLegacy
)

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
        VTK::IOGeometry
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

# ------------------------------------------------------------
# ✅ Add Unit Tests
# ------------------------------------------------------------
enable_testing()

# Add GoogleTest (if installed via conda or system package)
# Add GoogleTest (try system package first, otherwise fetch googletest)
include(FetchContent)
find_package(GTest QUIET)
if(GTest_FOUND)
    message(STATUS "Found GTest: ${GTest_DIR}")
else()
    message(STATUS "GTest not found; Fetching googletest via FetchContent...")
    # Use a known googletest release
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/release-1.12.1.zip
    )
    # Recommended to avoid overriding MSVC runtime on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

# Now add the tests/ directory (it links to GTest::gtest_main)
add_subdirectory(tests)

# ------------------------------------------------------------
# ✅ (Optional) Add Examples
# ------------------------------------------------------------
#if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/CMakeLists.txt")
#    add_subdirectory(examples)
#endif()
