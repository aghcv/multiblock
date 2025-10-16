#pragma once

#include <string>
#include <openvdb/openvdb.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>

namespace fastvessels {

/// Convert vtkImageData (structured volume) to OpenVDB FloatGrid
openvdb::FloatGrid::Ptr ConvertVTKToVDB(vtkSmartPointer<vtkImageData> image);

/// Convert OpenVDB FloatGrid back to vtkImageData
vtkSmartPointer<vtkImageData> ConvertVDBToVTK(const openvdb::FloatGrid::Ptr& grid);

/// Detect direction from extension (.vtk/.vti → VDB or .vdb → VTK) and perform conversion
void ConvertDataSet(const std::string& inputPath, const std::string& outputPath);

} // namespace fastvessels
