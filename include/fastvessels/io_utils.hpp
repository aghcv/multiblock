#pragma once
#include <string>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

namespace fastvessels {

vtkSmartPointer<vtkPolyData> ReadPolyData(const std::string& filename);

}  // namespace fastvessels
