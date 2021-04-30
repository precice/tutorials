#include "SolidSolver.h"
#include <iostream>
#include <stdlib.h>
#include "precice/SolverInterface.hpp"

int main(int argc, char **argv)
{
  std::cout << "Starting Solid Solver..." << std::endl;
  using namespace precice;
  using namespace precice::constants;

  if (argc != 2) {
    std::cout << "Fluid: Usage: " << argv[0] << " configurationFileName>" << std::endl;
    std::cout << std::endl;

    return -1;
  }

  std::string configFileName(argv[1]);
  int         domainSize  = 100; // N
  int         chunkLength = domainSize + 1;

  std::cout << "N: " << domainSize << std::endl;
  std::cout << "inputs: " << argc << std::endl;

  const std::string solverName = "Solid";

  SolverInterface interface(solverName, configFileName, 0, 1);
  std::cout << "preCICE configured..." << std::endl;

  int dimensions           = interface.getDimensions();
  int meshID               = interface.getMeshID("Solid-Nodes-Mesh");
  int crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID);
  int pressureID           = interface.getDataID("Pressure", meshID);

  std::vector<double> pressure(chunkLength, 0.0);
  std::vector<double> crossSectionLength(chunkLength, 1.0);
  std::vector<double> grid(dimensions * chunkLength);

  for (int i = 0; i < chunkLength; i++) {
    for (int j = 0; j < dimensions; j++) {
      grid[i * dimensions + j] = i * (1 - j);
    }
  }

  std::vector<int> vertexIDs(chunkLength);
  interface.setMeshVertices(meshID, chunkLength, grid.data(), vertexIDs.data());

  double t  = 0;
  std::cout << "Initialize preCICE..." << std::endl;
  double dt = interface.initialize();

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

    SolidComputeSolution(chunkLength, pressure.data(), crossSectionLength.data()); // Call Solver

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
  return 0;
}
