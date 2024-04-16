#include "SolidSolver.h"
#include <iostream>
#include <stdlib.h>
#include "precice/precice.hpp"

int main(int argc, char **argv)
{
  std::cout << "Starting Solid Solver..." << std::endl;
  using namespace precice;

  if (argc != 2) {
    std::cout << "Fluid: Usage: " << argv[0] << " configurationFileName>" << std::endl;
    std::cout << std::endl;

    return -1;
  }

  std::string  configFileName(argv[1]);
  const int    domainSize  = 100; // N
  const int    chunkLength = domainSize + 1;
  const double tubeLength  = 10;

  std::cout << "N: " << domainSize << std::endl;
  std::cout << "inputs: " << argc << std::endl;

  const std::string solverName = "Solid";

  precice::Participant interface(solverName, configFileName, 0, 1);
  std::cout << "preCICE configured..." << std::endl;

  auto      meshName               = "Solid-Nodes-Mesh";
  auto      crossSectionLengthName = "CrossSectionLength";
  auto      pressureName           = "Pressure";
  const int dimensions             = interface.getMeshDimensions(meshName);

  std::vector<double> pressure(chunkLength, 0.0);
  std::vector<double> crossSectionLength(chunkLength, 1.0);

  std::vector<double> grid(dimensions * chunkLength, 0.0);
  const double        dx = tubeLength / domainSize;
  for (int i = 0; i < chunkLength; ++i) {
    grid[i * dimensions] = dx * i;
  }

  std::vector<int> vertexIDs(chunkLength);
  interface.setMeshVertices(meshName, grid, vertexIDs);

  if (interface.requiresInitialData()) {
    interface.writeData(meshName, crossSectionLengthName, vertexIDs, crossSectionLength);
  }

  std::cout << "Initialize preCICE..." << std::endl;
  interface.initialize();

  while (interface.isCouplingOngoing()) {
    if (interface.requiresWritingCheckpoint()) {
    }
    double dt = interface.getMaxTimeStepSize();

    interface.readData(meshName, pressureName, vertexIDs, dt, pressure);

    SolidComputeSolution(chunkLength, pressure.data(), crossSectionLength.data()); // Call Solver

    interface.writeData(meshName, crossSectionLengthName, vertexIDs, crossSectionLength);

    interface.advance(dt);

    if (interface.requiresReadingCheckpoint()) { // i.e. fluid not yet converged
    }
  }

  std::cout << "Exiting SolidSolver" << std::endl;
  return 0;
}
