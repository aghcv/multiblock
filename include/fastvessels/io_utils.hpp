#pragma once
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <string>

namespace fastvessels {

vtkSmartPointer<vtkPolyData> ReadPolyData(const std::string& filename);

}
