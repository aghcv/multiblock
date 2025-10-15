file(GLOB TEST_SOURCES CONFIGURE_DEPENDS *.cpp)

add_executable(multiblock_tests ${TEST_SOURCES})
target_link_libraries(multiblock_tests
    PRIVATE
        multiblock_lib
        VTK::IOGeometry        # explicitly link it here too
)

include(GoogleTest)

# Link to gtest_main in a portable way: prefer the namespaced target from find_package,
# otherwise use the target provided by FetchContent (gtest_main).
if(TARGET GTest::gtest_main)
    target_link_libraries(multiblock_tests PRIVATE GTest::gtest_main)
elseif(TARGET gtest_main)
    target_link_libraries(multiblock_tests PRIVATE gtest_main)
else()
    message(FATAL_ERROR "gtest_main target not found. Ensure GTest is available or FetchContent succeeded.")
endif()

gtest_discover_tests(multiblock_tests)
