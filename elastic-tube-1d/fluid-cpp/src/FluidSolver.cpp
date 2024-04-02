#include "FluidComputeSolution.h"
#include "utilities.h"

#include <cmath>
#include <iostream>
#include <vector>
#include "precice/precice.hpp"

using namespace precice;

int main(int argc, char **argv)
{

  std::cout << "Starting Fluid Solver..." << std::endl;
  if (argc != 2) {
    std::cout << std::endl;
    std::cout << "Fluid: Usage: " << argv[0] << " <configurationFileName>" << std::endl;

    return -1;
  }

  std::string  configFileName(argv[1]);
  int          domainSize  = 100; // N
  int          chunkLength = domainSize + 1;
  const double kappa       = 100;
  const double L           = 10.0; // tube length

  const std::string solverName = "Fluid";

  std::string outputFilePrefix = "./output/out_fluid"; // extra

  precice::Participant interface(solverName, configFileName, 0, 1);
  std::cout << "preCICE configured..." << std::endl;

  auto      meshName               = "Fluid-Nodes-Mesh";
  auto      pressureName           = "Pressure";
  auto      crossSectionLengthName = "CrossSectionLength";
  const int dimensions             = interface.getMeshDimensions(meshName);

  std::vector<int> vertexIDs(chunkLength);

  const double PI = 3.141592653589793;

  const double r0        = 1 / sqrt(PI);         // radius of the tube
  const double a0        = std::pow(r0, 2) * PI; // cross sectional area
  const double u0        = 10;                   // mean velocity
  const double ampl      = 3;                    // amplitude of varying velocity
  const double frequency = 10;                   // frequency of variation
  const double t_shift   = 0;                    // temporal shift of variation
  const double p0        = 0;                    // pressure at outlet
  const double vel_in_0  = u0 + ampl * sin(frequency * (t_shift) *PI);

  std::vector<double> pressure(chunkLength, p0);
  std::vector<double> pressure_old(pressure);
  std::vector<double> crossSectionLength(chunkLength, a0);
  std::vector<double> crossSectionLength_old(crossSectionLength);
  std::vector<double> velocity(chunkLength, vel_in_0);
  std::vector<double> velocity_old(velocity);
  std::vector<double> grid(dimensions * chunkLength);

  const double cellwidth = (L / domainSize);
  for (int i = 0; i < chunkLength; i++) {
    for (int d = 0; d < dimensions; d++) {
      if (d == 0) {
        grid[i * dimensions] = i * cellwidth;
      } else {
        grid[i * dimensions + d] = 0.0;
      }
    }
  }

  interface.setMeshVertices(meshName, grid, vertexIDs);

  if (interface.requiresInitialData()) {
    interface.writeData(meshName, pressureName, vertexIDs, pressure);
  }

  double t = 0.0;
  std::cout << "Initialize preCICE..." << std::endl;
  interface.initialize();

  interface.readData(meshName, crossSectionLengthName, vertexIDs, 0, crossSectionLength);

  std::copy(crossSectionLength.begin(), crossSectionLength.end(), crossSectionLength_old.begin());

  // initialize such that mass conservation is fulfilled
  for (int i = 0; i < chunkLength; ++i) {
    velocity_old[i] = vel_in_0 * crossSectionLength_old[0] / crossSectionLength_old[i];
  }

  int out_counter = 0;

  while (interface.isCouplingOngoing()) {
    if (interface.requiresWritingCheckpoint()) {
    }

    auto dt = interface.getMaxTimeStepSize();

    fluidComputeSolutionSerial(
        // values from last time window
        velocity_old.data(), pressure_old.data(), crossSectionLength_old.data(),
        // last received crossSectionLength
        crossSectionLength.data(),
        t + dt, // used for inlet velocity
        domainSize,
        kappa,
        dt, // tau
        // resulting velocity pressure
        velocity.data(),
        pressure.data());

    interface.writeData(meshName, pressureName, vertexIDs, pressure);

    interface.advance(dt);

    interface.readData(meshName, crossSectionLengthName, vertexIDs, interface.getMaxTimeStepSize(), crossSectionLength);

    if (interface.requiresReadingCheckpoint()) {
    } else {
      t += dt;
      write_vtk(t, out_counter, outputFilePrefix.c_str(), chunkLength, grid.data(), velocity.data(), pressure.data(), crossSectionLength.data());
      for (int i = 0; i < chunkLength; i++) {
        crossSectionLength_old[i] = crossSectionLength[i];
        pressure_old[i]           = pressure[i];
        velocity_old[i]           = velocity[i];
      }
      out_counter++;
    }
  }

  std::cout << "Exiting FluidSolver" << std::endl;
  return 0;
}
