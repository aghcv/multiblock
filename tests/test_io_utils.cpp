#include <gtest/gtest.h>
#include "fastvessels/io_utils.hpp"
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSet.h>

// ------------------------------------------------------------------
// Test: CanReadSimpleSTL
// ------------------------------------------------------------------
TEST(IOUtils, CanReadSimpleSTL)
{
    // Read using the unified function
    auto dataset = fastvessels::ReadDataSet("../../tests/data/cube_1x1x1.stl");
    ASSERT_TRUE(dataset != nullptr);

    // Should be PolyData for STL files
    auto poly = vtkPolyData::SafeDownCast(dataset);
    ASSERT_TRUE(poly != nullptr);

    // Basic geometry sanity checks
    EXPECT_GT(poly->GetNumberOfPoints(), 0);
    EXPECT_GT(poly->GetNumberOfPolys(), 0);
}

// ------------------------------------------------------------------
// Optional: CanReadUnstructuredVTU
// ------------------------------------------------------------------
TEST(IOUtils, CanReadUnstructuredVTU)
{
    // This test will run only if you have a sample VTU file
    auto dataset = fastvessels::ReadDataSet("../../tests/data/heart_sample.vtu");
    ASSERT_TRUE(dataset != nullptr);

    // VTU should load as an UnstructuredGrid
    auto grid = vtkUnstructuredGrid::SafeDownCast(dataset);
    ASSERT_TRUE(grid != nullptr);

    // Basic sanity checks
    EXPECT_GT(grid->GetNumberOfPoints(), 0);
    EXPECT_GT(grid->GetNumberOfCells(), 0);
}

// ------------------------------------------------------------------
// Optional: CanHandleUppercaseExtensions
// ------------------------------------------------------------------
TEST(IOUtils, HandlesUppercaseExtensions)
{
    auto dataset = fastvessels::ReadDataSet("../../tests/data/CUBE_1X1X1.STL");
    ASSERT_TRUE(dataset != nullptr);

    // Case-insensitivity check
    auto poly = vtkPolyData::SafeDownCast(dataset);
    ASSERT_TRUE(poly != nullptr);
}
