#include "fastvessels/io_utils.hpp"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <stdexcept>

namespace fastvessels {

vtkSmartPointer<vtkPolyData> ReadPolyData(const std::string& filename) {
    std::string extension = filename.substr(filename.find_last_of('.') + 1);
    vtkSmartPointer<vtkPolyData> polydata;

    if (extension == "vtk") {
        auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName(filename.c_str());
        reader->Update();
        polydata = reader->GetOutput();
    } else if (extension == "vtp") {
        auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName(filename.c_str());
        reader->Update();
        polydata = reader->GetOutput();
    } else {
        throw std::runtime_error("Unsupported file format: " + extension);
    }

    return polydata;
}

}  // namespace fastvessels
