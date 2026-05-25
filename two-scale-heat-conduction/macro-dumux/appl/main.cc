// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include "properties.hh"

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/method.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/loadsolution.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <list>
#include <map>

int main(int argc, char **argv)
{
  using namespace Dumux;

  // define the type tag for this problem
  using TypeTag = Properties::TTag::TYPETAG;

  // make some types available
  using Scalar = GetPropType<TypeTag, Properties::Scalar>;

  // initialize MPI, finalize is done automatically on exit
  const auto &mpiHelper = Dune::MPIHelper::instance(argc, argv);

  // print dumux start message
  if (mpiHelper.rank() == 0)
    DumuxMessage::print(/*firstCall=*/true);

  // parse command line arguments and input file
  Parameters::init(argc, argv);
  Parameters::print();

  // try to create a grid (from the given grid file or the input file)
  GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
  gridManager.init();

  // we compute on the leaf grid view
  const auto &leafGridView = gridManager.grid().leafGridView();

  // create the finite volume grid geometry
  using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
  auto gridGeometry  = std::make_shared<GridGeometry>(leafGridView);

  // the problem (initial and boundary conditions)
  using Problem = GetPropType<TypeTag, Properties::Problem>;
  auto problem  = std::make_shared<Problem>(gridGeometry);

  ////////////////////////
  // Initialize preCICE //
  ////////////////////////

  // Tell preCICE about:
  // - Name of solver
  // - What rank of how many ranks this instance is
  // Configure preCICE. For now the config file is hardcoded.
  std::string meshName;

  auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();

  const auto runWithCoupling = getParam<bool>("Precice.RunWithCoupling");

  if (runWithCoupling) {
    couplingParticipant.announceConfig(mpiHelper.rank(), mpiHelper.size());
    meshName = couplingParticipant.getMeshNames()[0];
    // verify that dimensions match
    const int preciceDim = couplingParticipant.getMeshDimensions(meshName);
    const int dim        = int(leafGridView.dimension);
    std::cout << "Coupling dims = " << preciceDim << " , leafgrid dims = " << dim << std::endl;
    if (preciceDim != dim)
      DUNE_THROW(Dune::InvalidStateException, "Dimensions do not match");
  }

  // get mesh coordinates
  std::vector<double> coords;
  std::vector<int>    coupledElementIdxs;

  // coordinate loop (created vectors are 1D)
  // these positions of cell centers are later communicated to precice
  std::cout << "Coordinates: " << std::endl;
  for (const auto &element : elements(leafGridView, Dune::Partitions::interior)) {
    auto fvGeometry = localView(*gridGeometry);
    fvGeometry.bindElement(element);
    for (const auto &scv : scvs(fvGeometry)) {
      coupledElementIdxs.push_back(scv.elementIndex());
      const auto &pos = scv.center();
      for (const auto p : pos) {
        coords.push_back(p);
        std::cout << p << "  ";
      }
      std::cout << " ;" << std::endl;
    }
  }

  std::cout << "Number of Coupled Cells:" << coupledElementIdxs.size() << std::endl;

  int numberOfElements;
  if (runWithCoupling) {
    numberOfElements = coords.size() / couplingParticipant.getMeshDimensions(meshName);
  } else {
    numberOfElements = coupledElementIdxs.size();
  }

  if (runWithCoupling) {
    couplingParticipant.setMesh(meshName, coords);

    // couples between dumux element indices and preciceIndices;
    couplingParticipant.createIndexMapping(coupledElementIdxs);
  }

  // initialize the coupling data
  std::string readDatak00;
  std::string readDatak01;
  std::string readDatak10;
  std::string readDatak11;
  std::string readDataPorosity;
  std::string writeDataConcentration;

  if (runWithCoupling) {
    readDatak00            = couplingParticipant.getReadDataNamesOnMesh(meshName)[0];
    readDatak01            = couplingParticipant.getReadDataNamesOnMesh(meshName)[1];
    readDatak10            = couplingParticipant.getReadDataNamesOnMesh(meshName)[2];
    readDatak11            = couplingParticipant.getReadDataNamesOnMesh(meshName)[3];
    readDataPorosity       = couplingParticipant.getReadDataNamesOnMesh(meshName)[4];
    writeDataConcentration = couplingParticipant.getWriteDataNamesOnMesh(meshName)[0];
  }

  // the solution vector (initialized with zeros) NElements x 2(pressure,
  // temperature)
  using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
  SolutionVector x(gridGeometry->numDofs());
  problem->applyInitialSolution(x);
  auto xOld = x;

  // initialize the coupling data
  std::vector<double> temperatures;
  for (int solIdx = 0; solIdx < numberOfElements; ++solIdx) {
    temperatures.push_back(x[solIdx][problem->returnTemperatureIdx()]);
  };

  if (runWithCoupling) {
    couplingParticipant.writeQuantityVector(meshName, writeDataConcentration,
                                            temperatures);
    if (couplingParticipant
            .requiresToWriteInitialData()) { // not called in our example
      //    couplingParticipant.writeQuantityToOtherSolver(meshName,
      //    writeDataTemperature);
    }

    // initialize other data since they are not allowed to be read from micro
    // side
    std::vector<double> kInitial(numberOfElements, 1.0);
    std::vector<double> porosityInitial(numberOfElements, 0.5);
    couplingParticipant.writeQuantityVector(meshName, readDatak00, kInitial);
    couplingParticipant.writeQuantityVector(meshName, readDatak01, kInitial);
    couplingParticipant.writeQuantityVector(meshName, readDatak10, kInitial);
    couplingParticipant.writeQuantityVector(meshName, readDatak11, kInitial);
    couplingParticipant.writeQuantityVector(meshName, readDataPorosity,
                                            porosityInitial);
  }

  // the grid variables
  using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
  auto gridVariables  = std::make_shared<GridVariables>(problem, gridGeometry);
  gridVariables->init(x);

  // initialize the vtk output module
  using IOFields = GetPropType<TypeTag, Properties::IOFields>;
  VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x,
                                                           problem->name());
  IOFields::initOutputModule(vtkWriter);
  // add model specific output fields
  vtkWriter.addField(problem->getPorosity(), "Porosity");
  vtkWriter.addField(problem->getK00(), "K00");
  vtkWriter.addField(problem->getK01(), "K01");
  vtkWriter.addField(problem->getK10(), "K10");
  vtkWriter.addField(problem->getK11(), "K11");
  problem->updateVtkOutput(x);
  vtkWriter.write(0.0);

  // output every vtkOutputInterval time step
  const int vtkOutputInterval = getParam<int>("TimeLoop.OutputInterval");

  // initialize preCICE
  if (runWithCoupling) {
    couplingParticipant.initialize();
  }

  // time loop parameters
  const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
  double     preciceDt;
  double     solverDt;
  double     dt;

  if (runWithCoupling) {
    preciceDt = couplingParticipant.getMaxTimeStepSize();
    solverDt  = getParam<Scalar>("TimeLoop.InitialDt");
    dt        = std::min(preciceDt, solverDt);
  } else {
    dt = getParam<Scalar>("TimeLoop.InitialDt");
  }

  // instantiate time loop
  auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
  timeLoop->setMaxTimeStepSize(getParam<Scalar>("TimeLoop.MaxDt"));

  // initialize adapter checkpointing
  if (runWithCoupling) {
    couplingParticipant.initializeCheckpoint(x, *gridVariables, *timeLoop);
  }

  // the assembler with time loop for instationary problem
  using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
  auto assembler  = std::make_shared<Assembler>(problem, gridGeometry,
                                               gridVariables, timeLoop, xOld);

  // the linear solver
  using LinearSolver =
      ILUBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>,
                            LinearAlgebraTraitsFromAssembler<Assembler>>;
  auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(),
                                                     gridGeometry->dofMapper());

  // the non-linear solver
  using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
  NewtonSolver nonLinearSolver(assembler, linearSolver);

  // time loop
  int n_out = 0; // counts timesteps for the output interval
  std::cout << "Time Loop starts" << std::endl;
  timeLoop->start();
  do {
    if (runWithCoupling) {
      if (couplingParticipant.isCouplingOngoing() == false)
        break;

      // write checkpoint
      couplingParticipant.writeCheckpointIfRequired();

      preciceDt = couplingParticipant.getMaxTimeStepSize();
      solverDt  = std::min(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()),
                           timeLoop->maxTimeStepSize());
      dt        = std::min(preciceDt, solverDt);

      // read porosity and conductivity data from other solver
      // TODO: data needs to be updated if Newton solver adapts time-step size
      // and coupling data is interpolated in time
      couplingParticipant.readQuantityFromOtherSolver(meshName, readDatak00,
                                                      dt);
      couplingParticipant.readQuantityFromOtherSolver(meshName, readDatak01,
                                                      dt);
      couplingParticipant.readQuantityFromOtherSolver(meshName, readDatak10,
                                                      dt);
      couplingParticipant.readQuantityFromOtherSolver(meshName, readDatak11,
                                                      dt);
      couplingParticipant.readQuantityFromOtherSolver(meshName,
                                                      readDataPorosity, dt);
      // store coupling data in spatial params, exchange with MPI
      problem->spatialParams().updateCouplingData();
    } else {
      dt = std::min(
          nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()),
          timeLoop->maxTimeStepSize());
    }
    // set new dt as suggested by the Newton solver or by preCICE
    timeLoop->setTimeStepSize(dt);

    std::cout << "nonLinearSolver starts with target dt: " << dt << std::endl;

    // linearize & solve
    nonLinearSolver.solve(x, *timeLoop);

    // save dt value that was actually used by the non-linear solver
    dt = timeLoop->timeStepSize();

    // DuMux advance + report
    gridVariables->advanceTimeStep();
    timeLoop->advanceTimeStep();
    timeLoop->reportTimeStep();
    xOld = x;

    if (runWithCoupling) {
      int solIdx = 0;
      for (const auto &element : elements(leafGridView, Dune::Partitions::interior)) {
        auto fvGeometry = localView(*gridGeometry);
        fvGeometry.bindElement(element);
        for (const auto &scv : scvs(fvGeometry)) {
          temperatures[solIdx++] =
              x[scv.elementIndex()][problem->returnTemperatureIdx()];
        }
      }

      couplingParticipant.writeQuantityVector(meshName,
                                              writeDataConcentration, temperatures);
      couplingParticipant.writeQuantityToOtherSolver(meshName,
                                                     writeDataConcentration);

      // advance preCICE
      if ((!fabs(preciceDt - dt)) < 1e-14) {
        std::cout << "dt from preCICE is different than dt from DuMuX."
                  << " preCICE dt = " << preciceDt
                  << " and DuMuX dt = " << solverDt
                  << " resulted in dt = " << dt
                  << std::endl;
      }
      std::flush(std::cout);
      couplingParticipant.advance(dt);

      // reset to checkpoint if not converged
      if (couplingParticipant.readCheckpointIfRequired()) {
        gridVariables->advanceTimeStep();
        xOld = x;
        continue;
      }
    }

    if (runWithCoupling) {
      // if coupling, write VTK output only when time window is complete
      if (couplingParticipant.isTimeWindowComplete()) {
        n_out += 1;
        if (n_out == vtkOutputInterval) {
          problem->updateVtkOutput(x);
          vtkWriter.write(timeLoop->time());
          n_out = 0;
        }
      }
    } else {
      n_out += 1;
      if (n_out == vtkOutputInterval) {
        problem->updateVtkOutput(x);
        vtkWriter.write(timeLoop->time());
        n_out = 0;
      }
    }

    std::cout << "Time: " << timeLoop->time() << std::endl;

  } while (!timeLoop->finished());

  timeLoop->finalize(leafGridView.comm());

  ////////////////////////////////////////////////////////////
  // finalize, print dumux message to say goodbye
  ////////////////////////////////////////////////////////////
  if (runWithCoupling) {
    couplingParticipant.finalize();
  }
  // print dumux end message
  if (mpiHelper.rank() == 0) {
    Parameters::print();
    DumuxMessage::print(/*firstCall=*/false);
  }

  return 0;
} // end main
