#include <vtkXMLPStructuredGridWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>


//VTK Library
#include <vtkXMLPStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkMPIController.h>
#include <vtkProgrammableFilter.h>
#include <vtkInformation.h>

struct Args {
  vtkProgrammableFilter* pf;
  int local_extent[6];
};

// function to operate on the point attribute data
void execute (void* arg) {
  Args* args = reinterpret_cast<Args*>(arg);
  auto info = args->pf->GetOutputInformation(0);
  auto output_tmp = args->pf->GetOutput();
  auto input_tmp  = args->pf->GetInput();
  vtkStructuredGrid* output = dynamic_cast<vtkStructuredGrid*>(output_tmp);
  vtkStructuredGrid* input  = dynamic_cast<vtkStructuredGrid*>(input_tmp);
  output->ShallowCopy(input);
  output->SetExtent(args->local_extent);
}

void write_vtk(const std::string &filename, int *grid, int N, int rank, 
                    int size, int coords[3], int local_extent[6], 
                    int global_extent[6], int physdims[3], vtkMPIController* contr) {

    std::ostringstream rank_fname;
    rank_fname << filename << ".pvts";

    // Calculate the origin based on the process coordinates
    double origin[3];
    double spacing[3] = {1.0, 1.0, 1.0}; // Spacing between points in each dimension

    // Create points for the structured grid
    vtkNew<vtkPoints> points;
    for (int z = 0; z < physdims[2]; ++z) {
        for (int y = 0; y < physdims[1]; ++y) {
            for (int x = 0; x < physdims[0]; ++x) {
                double px = local_extent[0] + x * spacing[0];
                double py = local_extent[2] + y * spacing[1];
                double pz = local_extent[4] + z * spacing[2];
                points->InsertNextPoint(px, py, pz);
            }
        }
    }

    // Assign points to the grid
    //structuredGrid->SetPoints(points);

    // Add a scalar field (no tuples)
    vtkNew<vtkIntArray> scalarField;

    scalarField->SetName("ScalarField"); // Name of the scalar field
    scalarField->SetNumberOfComponents(1); // Single scalar value per point
    scalarField->SetNumberOfTuples((physdims[0]-1)*(physdims[1]-1)*(physdims[2]-1));
    
    for (vtkIdType i = 0; i < (physdims[0]-1)*(physdims[1]-1)*(physdims[2]-1); ++i) {
        scalarField->SetValue(i,grid[i]); // Assign a unique value per point
    }
    
    vtkNew<vtkProgrammableFilter> pf;
    
    Args args;
    args.pf = pf;
    for(int i=0; i<6; ++i) args.local_extent[i] = local_extent[i];

    pf->SetExecuteMethod(execute, &args);

    // Create a structured grid and assign point data and cell data to it
    vtkNew<vtkStructuredGrid> structuredGrid;
    
    structuredGrid->SetExtent(global_extent);
    pf->SetInputData(structuredGrid);
    structuredGrid->SetPoints(points);
    structuredGrid->GetCellData()->AddArray(scalarField);

    auto parallel_writer = vtkSmartPointer<vtkXMLPStructuredGridWriter>::New();
    parallel_writer->SetInputConnection(pf->GetOutputPort());
    parallel_writer->SetController(contr);
    parallel_writer->SetFileName(rank_fname.str().c_str());
    parallel_writer->SetNumberOfPieces(size);
    parallel_writer->SetStartPiece(rank);
    parallel_writer->SetEndPiece(rank);
    parallel_writer->SetDataModeToBinary();
    parallel_writer->Update();
    parallel_writer->Write();

}