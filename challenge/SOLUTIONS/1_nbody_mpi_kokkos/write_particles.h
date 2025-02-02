#include <Kokkos_Core.hpp>

using View = Kokkos::View<double**, Kokkos::Device<Kokkos::DefaultExecutionSpace,Kokkos::SharedSpace>>;
typedef Kokkos::MDRangePolicy<Kokkos::Rank<2>> mdrange_policy;    // particle force calculation "loop"


#include <vtkXMLPPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPointData.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

//VTK Library
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkMPIController.h>
#include <vtkProgrammableFilter.h>
#include <vtkInformation.h>

struct Args {
  vtkProgrammableFilter* pf;
};

// function to operate on the point attribute data
void execute (void* arg) {
  Args* args = reinterpret_cast<Args*>(arg);
  auto info = args->pf->GetOutputInformation(0);
  auto output_tmp = args->pf->GetOutput();
  auto input_tmp  = args->pf->GetInput();
  vtkPolyData* output = dynamic_cast<vtkPolyData*>(output_tmp);
  vtkPolyData* input  = dynamic_cast<vtkPolyData*>(input_tmp);
  output->ShallowCopy(input);
}


void write_vtk_polydata(const std::string &filename, 
                        const View &particles, const int Np,
                        int rank, int size, 
                        vtkMPIController* contr) {
    std::ostringstream rank_fname;
    rank_fname << filename << ".pvtp";

    // Create points for the particles
    vtkNew<vtkPoints> points;
    for (int i = 0; i < Np; i++) {
        points->InsertNextPoint(particles(i,1), particles(i,2), particles(i,3));
    }

    // Create a scalar field for the particles
    vtkNew<vtkDoubleArray> ranks;
    vtkNew<vtkDoubleArray> velocities;
    vtkNew<vtkDoubleArray> forces;
    velocities->SetName("Velocity"); // Name of the scalar field
    velocities->SetNumberOfComponents(3); // Single scalar value per point
    forces->SetName("Force"); // Name of the scalar field
    forces->SetNumberOfComponents(3); // Single scalar value per point
    ranks->SetName("Rank"); // Name of the scalar field
    ranks->SetNumberOfComponents(1); // Single scalar value per point

    for (int i = 0; i < Np; i++) {
      velocities->InsertNextTuple3(particles(i,4), particles(i,5), particles(i,6));
      forces->InsertNextTuple3(particles(i,7), particles(i,8), particles(i,9));
      ranks->InsertNextValue(particles(i,10));
    }


    vtkNew<vtkProgrammableFilter> pf;
    
    Args args;
    args.pf = pf;

    pf->SetExecuteMethod(execute, &args);

    // Create polydata and add points and scalar field
    vtkNew<vtkPolyData> polyData;
    
    pf->SetInputData(polyData);

    polyData->SetPoints(points);
    polyData->GetPointData()->AddArray(velocities);
    polyData->GetPointData()->AddArray(forces);
    polyData->GetPointData()->AddArray(ranks);

    // Create and configure the parallel writer
    vtkNew<vtkXMLPPolyDataWriter> parallel_writer;
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
