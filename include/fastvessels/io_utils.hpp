#pragma once
#include <string>
#include <vtkSmartPointer.h>
#include <vtkDataSet.h>

namespace fastvessels {

/**
 * @brief Read any VTK-compatible dataset (surface or volumetric).
 * 
 * Supported formats:
 *   - .stl → vtkPolyData
 *   - .vtp → vtkPolyData
 *   - .vtk → vtkPolyData or vtkUnstructuredGrid
 *   - .vtu → vtkUnstructuredGrid
 *
 * @param filename Path to the file.
 * @return vtkSmartPointer<vtkDataSet> (can be cast to vtkPolyData or vtkUnstructuredGrid)
 * @throws std::runtime_error if the file type is unsupported or reading fails.
 */
vtkSmartPointer<vtkDataSet> ReadDataSet(const std::string& filename);

}  // namespace fastvessels
