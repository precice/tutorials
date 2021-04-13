#include "SolidSolver.h"
#include <iostream>
#include <stdlib.h>
#include "mpi.h"
#include "precice/SolverInterface.hpp"

int main(int argc, char **argv)
{
  std::cout << "Starting Solid Solver..." << std::endl;
  using namespace precice;
  using namespace precice::constants;

  if (argc != 3 && argc != 4) {
    std::cout << "Fluid: Usage: mpiexec -np <#procs> " << argv[0] << " <configurationFileName> <N> -parallel" << std::endl;
    std::cout << "or" << std::endl;
    std::cout << "Usage: " << argv[0] << " configurationFileName> <N>" << std::endl;
    std::cout << std::endl;
    std::cout << "N:     Number of mesh elements, needs to be equal for fluid and Solid solver." << std::endl;

    return -1;
  }

  const bool parallel = (argc == 4);

  std::string configFileName(argv[1]);
  int         domainSize  = atoi(argv[2]); // N
  int         chunkLength = domainSize + 1;

  std::cout << "N: " << domainSize << std::endl;
  std::cout << "inputs: " << argc << std::endl;

  const std::string solverName = "Solid";

  int gridOffset, rank = 0, size = 1;
  if (parallel) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if ((domainSize + 1) % size == 0) {
      chunkLength = (domainSize + 1) / size;
      gridOffset  = rank * chunkLength;
    } else if (rank < (domainSize + 1) % size) {
      chunkLength = (domainSize + 1) / size + 1;
      gridOffset  = rank * chunkLength;
    } else {
      chunkLength = (domainSize + 1) / size;
      gridOffset  = ((domainSize + 1) % size) * ((domainSize + 1) / size + 1) + (rank - ((domainSize + 1) % size)) * (domainSize + 1) / size;
    }
  }

  SolverInterface interface(solverName, configFileName, rank, size);
  std::cout << "preCICE configured..." << std::endl;

  int dimensions           = interface.getDimensions();
  int meshID               = interface.getMeshID("Solid-Nodes-Mesh");
  int crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID);
  int pressureID           = interface.getDataID("Pressure", meshID);

  std::vector<double> pressure(chunkLength, 0.0);
  std::vector<double> crossSectionLength(chunkLength, 1.0);
  std::vector<double> grid(dimensions * chunkLength);

  if (parallel) {
    for (int i = 0; i < chunkLength; i++) {
      for (int j = 0; j < dimensions; j++) {
        grid[i * dimensions + j] = j == 0 ? gridOffset + (double) i : 0.0;
      }
    }
  } else {
    for (int i = 0; i < chunkLength; i++) {
      for (int j = 0; j < dimensions; j++) {
        grid[i * dimensions + j] = i * (1 - j);
      }
    }
  }

  std::vector<int> vertexIDs(chunkLength);
  interface.setMeshVertices(meshID, chunkLength, grid.data(), vertexIDs.data());
  std::cout << "Initialize preCICE..." << std::endl;
  interface.initialize();

  double t  = 0;
  double dt = 0.01;

  if (interface.isActionRequired(actionWriteInitialData())) {
    interface.writeBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs.data(), crossSectionLength.data());
    interface.markActionFulfilled(actionWriteInitialData());
  }

  interface.initializeData();

  while (interface.isCouplingOngoing()) {
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      interface.markActionFulfilled(actionWriteIterationCheckpoint());
    }

    if (interface.isReadDataAvailable()) {
      interface.readBlockScalarData(pressureID, chunkLength, vertexIDs.data(), pressure.data());
    }

    SolidComputeSolution(rank, size, chunkLength, pressure.data(), crossSectionLength.data()); // Call Solver

    if (interface.isWriteDataRequired(dt)) {
      interface.writeBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs.data(), crossSectionLength.data());
    }

    interface.advance(dt);

    if (interface.isActionRequired(actionReadIterationCheckpoint())) { // i.e. fluid not yet converged
      interface.markActionFulfilled(actionReadIterationCheckpoint());
    } else {
      t += dt;
    }
  }

  std::cout << "Exiting SolidSolver" << std::endl;
  interface.finalize();
  if (parallel) {
    MPI_Finalize();
  }
  return 0;
}
