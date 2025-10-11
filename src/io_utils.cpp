#include "fastvessels/io_utils.hpp"
#include <vtkSTLReader.h>
#include <vtkPLYReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

namespace fastvessels {

vtkSmartPointer<vtkPolyData> ReadPolyData(const std::string& filename) {
    std::string ext = fs::path(filename).extension().string();
    vtkSmartPointer<vtkPolyDataAlgorithm> reader;

    if (ext == ".stl") {
        reader = vtkSmartPointer<vtkSTLReader>::New();
    } else if (ext == ".ply") {
        reader = vtkSmartPointer<vtkPLYReader>::New();
    } else if (ext == ".vtk") {
        reader = vtkSmartPointer<vtkPolyDataReader>::New();
    } else if (ext == ".vtp") {
        reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    } else {
        throw std::runtime_error("Unsupported file format: " + ext);
    }

    reader->SetFileName(filename.c_str());
    reader->Update();

    auto poly = reader->GetOutput();
    if (!poly || poly->GetNumberOfPoints() == 0)
        throw std::runtime_error("Failed to read geometry or empty mesh.");

    std::cout << "✅ Loaded " << filename << " with "
              << poly->GetNumberOfPoints() << " points and "
              << poly->GetNumberOfCells() << " cells.\n";
    return poly;
}

} // namespace fastvessels
