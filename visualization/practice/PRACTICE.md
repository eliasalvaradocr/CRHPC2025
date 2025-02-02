# Basic VTK File Writing Examples for ParaView (C++)

This document provides instructions for writing VTK files (VTS, VTU, VTK) using the **VTK library** in C++. These files can be used with ParaView for visualization.

---

## Example 1: Writing a Sphere to a `.vtk` File

### Objective

Create a sphere and write it to a `.vtk` file using the **VTK library** in C++.

### Instructions

1. **Include VTK Headers:**
   - Ensure that VTK is installed and configured for your C++ environment.

2. **Define Geometry:**
   - Use `vtkSphereSource` to define the sphere.
   - Set parameters like the radius and resolution.

3. **Create PolyData:**
   - Use the `GetOutput()` method of `vtkSphereSource` to get the sphere geometry in `vtkPolyData`.

4. **Write to File:**
   - Use `vtkPolyDataWriter` to write the geometry to a `.vtk` file.

5. **Code Snippet:**
   ```cpp
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

# Example 2: Writing a Structured Grid to a .vts File 

## Objective

Write a structured grid to a .vtu file using the VTK library (XML-based VTK format).

## Instructions

1. Include VTK Headers:
   Ensure that VTK is installed and configured for your C++ environment.

2. Define the Grid:
   * Use vtkStructuredGrid to create the grid.
   * Set dimensions using SetDimensions().

3. Define Points:
   * Use vtkPoints to add point coordinates for the grid.

4. Write to File:

   * Use vtkXMLStructuredGridWriter to write the grid to a .vts file.

* Code Snippet:

```cpp
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
```

# Example 3: Writing PolyData (Line) to a .vtk File

## Objective

Create a simple 3D line and write it to a .vtk file in ASCII format.

## Instructions

1. Include VTK Headers:
 * Ensure that VTK is installed and configured for your C++ environment.

2. Define Points:

 * Use vtkPoints to specify the coordinates of the line's endpoints.

3. Create a Line:

* Use vtkCellArray to define the connectivity of the points (create a line connecting two points).

4. Create PolyData:

* Use vtkPolyData to store the points and lines.

5. Write to File:

* Use vtkPolyDataWriter to write the line to a .vtk file.

6. Code Snippet:

```cpp 
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>

int main()
{
    // Define points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(0.0, 0.0, 0.0);
    points->InsertNextPoint(1.0, 1.0, 1.0);

    // Define a line
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType line[2] = {0, 1};
    lines->InsertNextCell(2, line);

    // Create PolyData
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);

    // Write to VTK file
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName("line.vtk");
    writer->SetInputData(polyData);
    writer->Write();

    return 0;
}

```