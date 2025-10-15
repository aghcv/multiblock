#include "fastvessels/io_utils.hpp"

#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <stdexcept>
#include <algorithm>

namespace fastvessels {

vtkSmartPointer<vtkDataSet> ReadDataSet(const std::string& filename) {
    // Get lowercase extension
    auto dotPos = filename.find_last_of('.');
    if (dotPos == std::string::npos)
        throw std::runtime_error("File has no extension: " + filename);

    std::string ext = filename.substr(dotPos + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    // Reader pointer
    vtkSmartPointer<vtkDataSet> dataset;

    if (ext == "stl") {
        auto reader = vtkSmartPointer<vtkSTLReader>::New();
        reader->SetFileName(filename.c_str());
        reader->Update();
        dataset = reader->GetOutput();
    } 
    else if (ext == "vtp") {
        auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName(filename.c_str());
        reader->Update();
        dataset = reader->GetOutput();
    } 
    else if (ext == "vtu") {
        auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        reader->SetFileName(filename.c_str());
        reader->Update();
        dataset = reader->GetOutput();
    } 
    else if (ext == "vtk") {
        // Legacy .vtk can be PolyData or UnstructuredGrid
        auto polyReader = vtkSmartPointer<vtkPolyDataReader>::New();
        polyReader->SetFileName(filename.c_str());
        polyReader->Update();
        if (polyReader->GetOutput() && polyReader->GetOutput()->GetNumberOfPoints() > 0) {
            dataset = polyReader->GetOutput();
        } else {
            auto gridReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
            gridReader->SetFileName(filename.c_str());
            gridReader->Update();
            dataset = gridReader->GetOutput();
        }
    }
    else {
        throw std::runtime_error("Unsupported file format: " + ext);
    }

    if (!dataset)
        throw std::runtime_error("Failed to read file: " + filename);

    return dataset;
}

}  // namespace fastvessels
