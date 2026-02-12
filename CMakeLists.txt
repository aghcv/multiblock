cmake_minimum_required(VERSION 3.20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(multiblock LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

if(DEFINED ENV{CONDA_PREFIX})
    list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}/lib/cmake")
endif()

find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    IOImport
    IOXML
    IOLegacy
    IOGeometry
    FiltersCore
    FiltersGeometry
    RenderingCore
    RenderingOpenGL2
)

add_library(multiblock_lib STATIC
    src/model/common.cpp
    src/model/env_info.cpp
    src/model/config.cpp
    src/geometry/centerline_detection.cpp
    src/geometry/mblock_util.cpp
    src/geometry/io_utils.cpp
    src/simulation/solver.cpp
    src/preprocessing/preprocess.cpp
    src/results/postprocess.cpp
)
target_include_directories(multiblock_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(multiblock_lib PUBLIC
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::IOImport
    VTK::IOXML
    VTK::IOLegacy
    VTK::IOGeometry
    VTK::FiltersCore
    VTK::FiltersGeometry
    VTK::RenderingCore
    VTK::RenderingOpenGL2
)

add_executable(multiblock src/app/app_main.cpp)
target_include_directories(multiblock PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(multiblock PRIVATE multiblock_lib
    VTK::CommonCore VTK::CommonDataModel VTK::IOImport VTK::IOXML VTK::IOLegacy VTK::IOGeometry
    VTK::FiltersCore VTK::FiltersGeometry VTK::RenderingCore VTK::RenderingOpenGL2)

set_target_properties(multiblock PROPERTIES
    SUFFIX ".mbx"
)

vtk_module_autoinit(TARGETS multiblock multiblock_lib MODULES
    VTK::CommonCore VTK::CommonDataModel VTK::IOImport VTK::IOXML VTK::IOLegacy VTK::IOGeometry
    VTK::FiltersCore VTK::FiltersGeometry VTK::RenderingCore VTK::RenderingOpenGL2)

# No unit tests on the 'light' branch
