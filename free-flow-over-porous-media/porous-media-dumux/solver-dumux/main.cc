// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
/*!
 * \file
 *
 * \brief A test problem for the coupled Stokes/Darcy problem (1p)
 */
#include <config.h>

#include <ctime>
#include <fstream>
#include <iostream>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>

#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/method.hh>

#include <dumux/io/format.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include "properties.hh"

#include "dumux-precice/couplingadapter.hh"

/*!
 * \brief Returns the velocity at the interface using Darcy's law for reconstruction
 */
template <class FluxVariables,
          class Problem,
          class Element,
          class FVElementGeometry,
          class ElementVolumeVariables,
          class SubControlVolumeFace,
          class ElementFluxVariablesCache>
auto velocityAtInterface(const Problem                   &problem,
                         const Element                   &element,
                         const FVElementGeometry         &fvGeometry,
                         const ElementVolumeVariables    &elemVolVars,
                         const SubControlVolumeFace      &scvf,
                         const ElementFluxVariablesCache &elemFluxVarsCache)
{
  const int     phaseIdx = 0;
  FluxVariables fluxVars;
  fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf,
                elemFluxVarsCache);
  auto upwindTerm = [phaseIdx](const auto &volVars) {
    return volVars.mobility(phaseIdx);
  };
  const auto scalarVelocity =
      fluxVars.advectiveFlux(phaseIdx, upwindTerm) / scvf.area();
  return scalarVelocity;
}

template <class FluxVariables,
          class Problem,
          class GridVariables,
          class SolutionVector>
void setInterfaceVelocities(const Problem        &problem,
                            const GridVariables  &gridVars,
                            const SolutionVector &sol,
                            const std::string     meshName,
                            const std::string     dataName)
{
  const auto &gridGeometry      = problem.gridGeometry();
  auto        fvGeometry        = localView(gridGeometry);
  auto        elemVolVars       = localView(gridVars.curGridVolVars());
  auto        elemFluxVarsCache = localView(gridVars.gridFluxVarsCache());

  auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();

  for (const auto &element : elements(gridGeometry.gridView())) {
    fvGeometry.bind(element);
    elemVolVars.bind(element, fvGeometry, sol);
    elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

    for (const auto &scvf : scvfs(fvGeometry)) {
      if (couplingParticipant.isCoupledEntity(scvf.index())) {
        // TODO: What to do here?
        const double v = velocityAtInterface<FluxVariables>(
            problem, element, fvGeometry, elemVolVars, scvf,
            elemFluxVarsCache);
        couplingParticipant.writeScalarQuantityOnFace(
            meshName, dataName, scvf.index(), v);
      }
    }
  }
}

int main(int argc, char **argv)
{
  using namespace Dumux;

  // initialize MPI, finalize is done automatically on exit
  const auto &mpiHelper = Dune::MPIHelper::instance(argc, argv);

  // print dumux start message
  if (mpiHelper.rank() == 0)
    DumuxMessage::print(/*firstCall=*/true);

  // parse command line arguments and input file
  Parameters::init(argc, argv);

  using DarcyTypeTag = Properties::TTag::DarcyOneP;

  using DarcyGridManager =
      Dumux::GridManager<GetPropType<DarcyTypeTag, Properties::Grid>>;
  DarcyGridManager darcyGridManager;
  darcyGridManager.init("Darcy"); // pass parameter group

  // we compute on the leaf grid view
  const auto &darcyGridView = darcyGridManager.grid().leafGridView();

  // create the finite volume grid geometry
  using DarcyGridGeometry =
      GetPropType<DarcyTypeTag, Properties::GridGeometry>;
  auto darcyGridGeometry = std::make_shared<DarcyGridGeometry>(darcyGridView);
  darcyGridGeometry->update(darcyGridManager.grid().leafGridView());

  using DarcyProblem = GetPropType<DarcyTypeTag, Properties::Problem>;
  auto darcyProblem  = std::make_shared<DarcyProblem>(darcyGridGeometry);

  // the solution vector
  using SolutionVector =
      GetPropType<DarcyTypeTag, Properties::SolutionVector>;
  SolutionVector sol;
  sol.resize(darcyGridGeometry->numDofs());

  // Initialize preCICE.Tell preCICE about:
  std::string preciceConfigFilename = "../precice-config.xml";

  auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();
  couplingParticipant.announceConfig(mpiHelper.rank(), mpiHelper.size());

  const std::string meshName = couplingParticipant.getMeshNames()[0];
  const int         dim      = couplingParticipant.getMeshDimensions(meshName);
  std::cout << dim << "  " << int(DarcyGridGeometry::GridView::dimension)
            << std::endl;
  if (dim != int(DarcyGridGeometry::GridView::dimension))
    DUNE_THROW(Dune::InvalidStateException, "Dimensions do not match");

  // GET mesh coordinates
  const double xMin =
      getParam<std::vector<double>>("Grid.LowerLeft")[0];
  const double xMax =
      getParam<std::vector<double>>("Grid.UpperRight")[0];
  std::vector<double> coords; //( dim * vertexSize );
  std::vector<int>    coupledScvfIndices;

  for (const auto &element : elements(darcyGridView)) {
    auto fvGeometry = localView(*darcyGridGeometry);
    fvGeometry.bindElement(element);

    for (const auto &scvf : scvfs(fvGeometry)) {
      static constexpr auto eps = 1e-7;
      const auto           &pos = scvf.center();
      if (pos[1] > darcyGridGeometry->bBoxMax()[1] - eps) {
        if (pos[0] > xMin - eps && pos[0] < xMax + eps) {
          coupledScvfIndices.push_back(scvf.index());
          for (const auto p : pos)
            coords.push_back(p);
        }
      }
    }
  }

  couplingParticipant.setMesh(meshName, coords);
  couplingParticipant.createIndexMapping(coupledScvfIndices);

  const std::string dataNameV = couplingParticipant.getWriteDataNamesOnMesh(meshName)[0];
  const std::string dataNameP = couplingParticipant.getReadDataNamesOnMesh(meshName)[0];

  darcyProblem->applyInitialSolution(sol);

  // the grid variables
  using DarcyGridVariables =
      GetPropType<DarcyTypeTag, Properties::GridVariables>;
  auto darcyGridVariables =
      std::make_shared<DarcyGridVariables>(darcyProblem, darcyGridGeometry);
  darcyGridVariables->init(sol);

  // intialize the vtk output module

  VtkOutputModule<DarcyGridVariables,
                  GetPropType<DarcyTypeTag, Properties::SolutionVector>>
      darcyVtkWriter(*darcyGridVariables, sol, darcyProblem->name());
  using DarcyVelocityOutput =
      GetPropType<DarcyTypeTag, Properties::VelocityOutput>;
  darcyVtkWriter.addVelocityOutput(
      std::make_shared<DarcyVelocityOutput>(*darcyGridVariables));
  GetPropType<DarcyTypeTag, Properties::IOFields>::initOutputModule(
      darcyVtkWriter);
  darcyVtkWriter.write(0.0);

  using FluxVariables = GetPropType<DarcyTypeTag, Properties::FluxVariables>;
  if (couplingParticipant.requiresToWriteInitialData()) {
    setInterfaceVelocities<FluxVariables>(
        *darcyProblem, *darcyGridVariables, sol, meshName, dataNameV);
    couplingParticipant.writeQuantityToOtherSolver(meshName, dataNameV);
  }
  couplingParticipant.initialize();

  // initialize checkpointing
  couplingParticipant.initializeCheckpoint(sol, *darcyGridVariables);

  // the assembler for a stationary problem
  using Assembler = FVAssembler<DarcyTypeTag, DiffMethod::numeric>;
  auto assembler  = std::make_shared<Assembler>(
      darcyProblem, darcyGridGeometry, darcyGridVariables);

  // the linear solver
  using LinearSolver =
      UMFPackIstlSolver<SeqLinearSolverTraits,
                        LinearAlgebraTraitsFromAssembler<Assembler>>;
  auto linearSolver = std::make_shared<LinearSolver>();

  // the non-linear solver
  using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
  NewtonSolver nonLinearSolver(assembler, linearSolver);

  double preciceDt = couplingParticipant.getMaxTimeStepSize();
  auto   dt        = preciceDt;

  double vtkTime = 1.0;

  while (couplingParticipant.isCouplingOngoing()) {
    couplingParticipant.writeCheckpointIfRequired();

    couplingParticipant.readQuantityFromOtherSolver(meshName, dataNameP,
                                                    dt);

    // solve the non-linear system
    nonLinearSolver.solve(sol);

    setInterfaceVelocities<FluxVariables>(
        *darcyProblem, *darcyGridVariables, sol, meshName, dataNameV);
    couplingParticipant.writeQuantityToOtherSolver(meshName, dataNameV);

    couplingParticipant.advance(dt);
    preciceDt = couplingParticipant.getMaxTimeStepSize();
    dt        = std::min(preciceDt, dt);

    if (!couplingParticipant.readCheckpointIfRequired()) {
      vtkTime += 1.;
      darcyVtkWriter.write(vtkTime);
    }
  }
  couplingParticipant.finalize();

  ////////////////////////////////////////////////////////////
  // finalize, print dumux message to say goodbye
  ////////////////////////////////////////////////////////////

  // print dumux end message
  if (mpiHelper.rank() == 0) {
    Parameters::print();
    DumuxMessage::print(/*firstCall=*/false);
  }

  return 0;
} // end main
