#include <gtest/gtest.h>
#include "fastvessels/io_utils.hpp"

TEST(IOUtils, CanReadSimpleSTL)
{
    auto mesh = fastvessels::ReadPolyData("../../tests/data/cube_1x1x1.stl");
    ASSERT_TRUE(mesh != nullptr);
}
