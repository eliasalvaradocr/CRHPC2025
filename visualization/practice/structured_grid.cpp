#include <vtkStructuredGrid.h>
#include <vtkPoints.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkSmartPointer.h>

int main()
{
    // Define structured grid
    vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();
    grid->SetDimensions(10, 10, 10);

    // Define points and add to grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            for (int k = 0; k < 10; k++) {
                points->InsertNextPoint(i, j, k);
            }
        }
    }
    grid->SetPoints(points);

    // Write to VTU file
    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName("structured_grid.vts");
    writer->SetInputData(grid);
    writer->Write();

    return 0;
}
