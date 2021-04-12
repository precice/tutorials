#include "utilities.h"
#include "FluidComputeSolution.h"

#include "precice/SolverInterface.hpp"
#include <iostream>
#include <mpi.h>
#include <vector>

using namespace precice;
using namespace precice::constants;

int main(int argc, char** argv)
{

  std::cout << "Starting Fluid Solver..." << std::endl;
  if (argc != 5 && argc != 6) {
    std::cout << std::endl;
    std::cout << "Fluid: Usage: mpiexec -np <#procs> " << argv[0] << " <configurationFileName> <N> <tau> <kappa> -parallel" << std::endl;
    std::cout << "or" << std::endl;
    std::cout << "Usage: " << argv[0] << " <configurationFileName> <N> <tau> <kappa>" << std::endl;
    std::cout << std::endl;
    std::cout << "N:     Number of mesh elements, needs to be equal for fluid and Solid solver." << std::endl;
    std::cout << "tau:   Dimensionless time step size." << std::endl;
    std::cout << "kappa: Dimensionless structural stiffness." << std::endl;

    return -1;
  }

  std::string configFileName(argv[1]);
  int domainSize = atoi(argv[2]); //N
  int chunkLength = domainSize + 1; //serial run
  double tau = atof(argv[3]);
  double kappa = atof(argv[4]);

  std::string solverName = "Fluid";

  std::string outputFilePrefix = "../postproc/out_fluid"; //extra

    
  int gridOffset, rank = 0, size = 1;


  if (argc == 6){
	  MPI_Init(&argc, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  	MPI_Comm_size(MPI_COMM_WORLD, &size);
    outputFilePrefix += std::to_string(rank); //extra

    if ((domainSize + 1) % size == 0) {
        chunkLength = (domainSize + 1) / size;
        gridOffset = rank * chunkLength;
    } else if (rank < (domainSize + 1) % size) {
        chunkLength = (domainSize + 1) / size + 1;
        gridOffset = rank * chunkLength;
      } else {
        chunkLength = (domainSize + 1) / size;
        gridOffset = ((domainSize + 1) % size) * ((domainSize + 1) / size + 1) + (rank - ((domainSize + 1) % size)) * (domainSize + 1) / size;
      }
  }



  SolverInterface interface(solverName, configFileName, rank , size);
  std::cout << "preCICE configured..." << std::endl;


  int dimensions = interface.getDimensions();
  int meshID = interface.getMeshID("Fluid-Nodes-Mesh");
  int pressureID = interface.getDataID("Pressure", meshID);
  int crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID);

  std::vector<double> pressure, crossSectionLength, velocity, crossSectionLength_old;
  std::vector<int> vertexIDs;
  std::vector<double> grid;

  grid.resize(dimensions * chunkLength);
  pressure.resize(chunkLength);
  vertexIDs.resize(chunkLength);
  crossSectionLength.resize(chunkLength);
  crossSectionLength_old.resize(chunkLength);
  velocity.resize(chunkLength);
  
  for (int i = 0; i < chunkLength; i++) {
    pressure[i] = 0.0;
    crossSectionLength[i] = 1.0;
    crossSectionLength_old[i] = 1.0;
    velocity[i] = 10.0;
  }

  if (argc==6){
    for (int i = 0; i < chunkLength; i++) {
      for (int j = 0; j < dimensions; j++) {
          grid[i * dimensions + j] = j == 0 ? gridOffset + (double)i : 0.0;
      }
    
    }
  } else {
    for (int i = 0; i < chunkLength; i++) {
      for (int j = 0; j < dimensions; j++) {
          grid[i * dimensions + j] = i * (1 - j) *0.1;
      }
    }
  }


  interface.setMeshVertices(meshID, chunkLength, grid.data(), vertexIDs.data());

  std::cout << "Initialize preCICE..." << std::endl;
  interface.initialize();

  double t = 0.0;
  double dt = 0.01;

  if (interface.isActionRequired(actionWriteInitialData())) {
    interface.writeBlockScalarData(pressureID, chunkLength, vertexIDs.data(), pressure.data());
    interface.markActionFulfilled(actionWriteInitialData());
  }

  interface.initializeData();

  if (interface.isReadDataAvailable()) {
    interface.readBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs.data(), crossSectionLength.data());
  }

  int out_counter = 0;

  while (interface.isCouplingOngoing()) {
    int convergenceCounter = 0;
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      interface.markActionFulfilled(actionWriteIterationCheckpoint());
    }

    if (argc == 6){
	    fluidComputeSolutionParallel(rank, size, domainSize, chunkLength, kappa, tau, 0.0, t+dt,
      pressure.data(), 
      crossSectionLength.data(), 
      velocity.data());
    } else {
      fluidComputeSolutionSerial(crossSectionLength.data(),
      crossSectionLength_old.data(),   
	    velocity.data(),           
	    pressure.data(),     
	    t, domainSize, kappa, tau); 
    }

    interface.writeBlockScalarData(pressureID, chunkLength, vertexIDs.data(), pressure.data());

    interface.advance(dt);

    interface.readBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs.data(), crossSectionLength.data());

    if (interface.isActionRequired(actionReadIterationCheckpoint())) { // i.e. not yet converged
      interface.markActionFulfilled(actionReadIterationCheckpoint());
      convergenceCounter++;
    } else {
      t += dt;
      write_vtk(t, out_counter, outputFilePrefix.c_str(), chunkLength, grid.data(), velocity.data(), pressure.data(), crossSectionLength.data());
      for (int i = 0; i < chunkLength; i++) {
        crossSectionLength_old[i] = crossSectionLength[i];
      }
      out_counter++;
    }
  }

  std::cout << "Exiting FluidSolver" << std::endl;
  interface.finalize();
  if (argc == 6){
    MPI_Finalize();

  }

  return 0;
}
