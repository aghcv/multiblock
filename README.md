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





# quickstart
