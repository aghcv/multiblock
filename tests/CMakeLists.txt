file(GLOB TEST_SOURCES CONFIGURE_DEPENDS *.cpp)

add_executable(multiblock_tests ${TEST_SOURCES})
target_link_libraries(multiblock_tests
    PRIVATE
        multiblock_lib
        GTest::gtest_main
        VTK::IOGeometry        # ✅ explicitly link it here too
)

include(GoogleTest)
gtest_discover_tests(multiblock_tests)
