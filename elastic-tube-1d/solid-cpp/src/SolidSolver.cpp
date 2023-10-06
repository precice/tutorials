#include "SolidSolver.h"
#include <iostream>
#include <stdlib.h>
#include "precice/SolverInterface.hpp"

int main(int argc, char **argv)
{
  std::cout << "Starting Solid Solver..." << std::endl;
  using namespace precice;

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
  auto meshName               = "Solid-Nodes-Mesh";
  auto crossSectionLengthName = "CrossSectionLength";
  auto pressureName           = "Pressure";

  std::vector<double> pressure(chunkLength, 0.0);
  std::vector<double> crossSectionLength(chunkLength, 1.0);
  std::vector<double> grid(dimensions * chunkLength);

  for (int i = 0; i < chunkLength; i++) {
    for (int j = 0; j < dimensions; j++) {
      grid[i * dimensions + j] = i * (1 - j);
    }
  }

  std::vector<int> vertexIDs(chunkLength);
  interface.setMeshVertices(meshName, chunkLength, grid.data(), vertexIDs.data());

  if (interface.requiresInitialData()) {
    interface.writeBlockScalarData(meshName, crossSectionLengthName, chunkLength, vertexIDs.data(), crossSectionLength.data());
  }

  std::cout << "Initialize preCICE..." << std::endl;
  double dt = interface.initialize();

  while (interface.isCouplingOngoing()) {
    if (interface.requiresWritingCheckpoint()) {
    }

    interface.readBlockScalarData(meshName, pressureName, chunkLength, vertexIDs.data(), pressure.data());

    SolidComputeSolution(chunkLength, pressure.data(), crossSectionLength.data()); // Call Solver

    interface.writeBlockScalarData(meshName, crossSectionLengthName, chunkLength, vertexIDs.data(), crossSectionLength.data());

    interface.advance(dt);

    if (interface.requiresReadingCheckpoint()) { // i.e. fluid not yet converged
    }
  }

  std::cout << "Exiting SolidSolver" << std::endl;
  interface.finalize();
  return 0;
}
