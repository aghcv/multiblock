#include <gtest/gtest.h>
#include "fastvessels/io_utils.hpp"
#include "fastvessels/vtk_vdb_converter.hpp"

TEST(VTKVDBConverter, RoundTripStructured)
{
    auto vol = fastvessels::ReadDataSet("../../tests/data/heart.vtk");
    auto grid = fastvessels::ConvertVTKToVDB(vol);
    ASSERT_TRUE(grid != nullptr);

    auto back = fastvessels::ConvertVDBToVTK(grid);
    ASSERT_EQ(back->GetDimensions()[0], vol->GetDimensions()[0]);
}
