#include "FluidComputeSolution.h"
#include "utilities.h"

#include <iostream>
#include <vector>
#include <cmath>
#include "precice/SolverInterface.hpp"

using namespace precice;
using namespace precice::constants;

int main(int argc, char **argv)
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
  int         domainSize  = atoi(argv[2]);  //N
  int         chunkLength = domainSize + 1; //serial run
  const double tau        = atof(argv[3]);
  const double kappa      = atof(argv[4]);
  const double L = 10.0; // tube length

  const std::string solverName = "Fluid";

  std::string outputFilePrefix = "./output/out_fluid"; //extra

  int gridOffset, rank = 0, size = 1;

  SolverInterface interface(solverName, configFileName, 0, 1);
  std::cout << "preCICE configured..." << std::endl;

  const int dimensions           = interface.getDimensions();
  const int meshID               = interface.getMeshID("Fluid-Nodes-Mesh");
  const int pressureID           = interface.getDataID("Pressure", meshID);
  const int crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID);

  std::vector<int>    vertexIDs(chunkLength);

  const double PI = 3.141592653589793;

  const double r0 = 1 / sqrt(PI); // radius of the tube
  const double a0 = std::pow(r0, 2) * PI; // cross sectional area
  const double E = 10000; // elasticity module
  const double u0 = 10; // mean velocity
  const double ampl = 3; // amplitude of varying velocity
  const double frequency = 10; // frequency of variation
  const double t_shift = 0; // temporal shift of variation
  const double p0 = 0; // pressure at outlet
  const double vel_in_0 = u0 + ampl * sin(frequency * (t_shift) * PI);

  std::vector<double> pressure(chunkLength, p0);
  std::vector<double> pressure_old(pressure);
  std::vector<double> crossSectionLength(chunkLength, a0);
  std::vector<double> crossSectionLength_old(crossSectionLength);
  std::vector<double> velocity(chunkLength, vel_in_0);
  std::vector<double> velocity_old(velocity);
  std::vector<double> grid(dimensions * chunkLength);

  const double cellwidth =(L / domainSize) ;
  for (int i = 0; i < chunkLength; i++) {
    for (int d = 0; d < dimensions; d++) {
      if (d == 0) {
        grid[i * dimensions] = i * cellwidth;
      } else {
        grid[i * dimensions + d] = 0.0;
      }
    }
  }

  interface.setMeshVertices(meshID, chunkLength, grid.data(), vertexIDs.data());

  std::cout << "Initialize preCICE..." << std::endl;
  interface.initialize();

  double t  = 0.0;
  double dt = 0.01;

  if (interface.isActionRequired(actionWriteInitialData())) {
    interface.writeBlockScalarData(pressureID, chunkLength, vertexIDs.data(), pressure.data());
    interface.markActionFulfilled(actionWriteInitialData());
  }

  interface.initializeData();

  if (interface.isReadDataAvailable()) {
    interface.readBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs.data(), crossSectionLength.data());
  }

  std::copy(crossSectionLength.begin(), crossSectionLength.end(), crossSectionLength_old.begin());

  // initialize such that mass conservation is fulfilled
  for(int i = 0; i < chunkLength; ++i) {
    velocity_old[i] = vel_in_0 * crossSectionLength_old[0] / crossSectionLength_old[i];
  }

  int out_counter = 0;

  while (interface.isCouplingOngoing()) {
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      interface.markActionFulfilled(actionWriteIterationCheckpoint());
    }

    
    if (interface.isReadDataAvailable()) {
      interface.readBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs.data(), crossSectionLength.data());
    }

    fluidComputeSolutionSerial(
        // values from last time window
        velocity_old.data(), pressure_old.data(), crossSectionLength_old.data(),
        // last received crossSectionLength
        crossSectionLength.data(),
        t+dt, // used for inlet velocity
        domainSize, 
        kappa, 
        dt, // tau
        // resulting velocity pressure
        velocity.data(),
        pressure.data());
    
    if (interface.isWriteDataRequired(dt)) {
      interface.writeBlockScalarData(pressureID, chunkLength, vertexIDs.data(), pressure.data());
    }
    
    interface.advance(dt);

    //interface.readBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs.data(), crossSectionLength.data());
    if (interface.isReadDataAvailable()) {
      interface.readBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs.data(), crossSectionLength.data());
    }

    if (interface.isActionRequired(actionReadIterationCheckpoint())) { // i.e. not yet converged
      interface.markActionFulfilled(actionReadIterationCheckpoint());
    } else {
      t += dt;
      write_vtk(t, out_counter, outputFilePrefix.c_str(), chunkLength, grid.data(), velocity.data(), pressure.data(), crossSectionLength.data());
      for (int i = 0; i < chunkLength; i++) {
        crossSectionLength_old[i] = crossSectionLength[i];
        pressure_old[i] = pressure[i];
        velocity_old[i] = velocity[i];
      }
      out_counter++;
    }
  }

  std::cout << "Exiting FluidSolver" << std::endl;
  interface.finalize();
  return 0;
}
