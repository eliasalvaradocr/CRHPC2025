#include <vtkSphereSource.h>
#include <vtkPolyDataWriter.h>
#include <vtkNew.h>

int main()
{
    // Define the sphere geometry
    vtkNew<vtkSphereSource> sphere;
    sphere->SetRadius(5.0);
    sphere->SetPhiResolution(50);
    sphere->SetThetaResolution(50);
    sphere->Update();

    // Write to VTK file
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName("sphere.vtk");
    writer->SetInputData(sphere->GetOutput());
    writer->Write();

    return 0;
}