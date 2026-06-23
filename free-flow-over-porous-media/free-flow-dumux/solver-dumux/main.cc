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
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/istl/io.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/partial.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/io/format.hh>
#include <dumux/io/grid/gridmanager.hh>

#include <dumux/freeflow/navierstokes/velocityoutput.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>

#include <dumux/multidomain/newtonsolver.hh>

#include "properties.hh"

#include "dumux-precice/couplingadapter.hh"

#include <iostream>
#include <memory>
#include <type_traits>

// TODO
//  Helper function to put pressure on interface

template <class ElementVolumeVariables, class SubControlVolumeFace>
auto velocityAtInterface(const ElementVolumeVariables &elemVolVars,
                         const SubControlVolumeFace   &scvf)
{
  assert(scvf.isFrontal());
  const double scalarVelocity = elemVolVars[scvf.insideScvIdx()].velocity();
  auto         velocity       = scvf.unitOuterNormal();
  velocity[scvf.normalAxis()] = scalarVelocity;
  return velocity;
}

template <typename TypeTag,
          class Problem,
          class Element,
          class SubControlVolumeFace,
          class FVElementGeometry,
          class ElementVolumeVariables,
          class ElementFluxVariablesCache>
auto pressureAtInterface(const Problem                   *problem,
                         const Element                   &element,
                         const SubControlVolumeFace      &scvf,
                         const FVElementGeometry         &fvGeometry,
                         const ElementVolumeVariables    &elemVolVars,
                         const ElementFluxVariablesCache &elemFluxVarsCache)
{
  using ElementBoundaryTypes =
      Dumux::GetPropType<TypeTag, Dumux::Properties::ElementBoundaryTypes>;
  using LocalResidual =
      Dumux::GetPropType<TypeTag, Dumux::Properties::LocalResidual>;
  using NumEqVector = Dumux::NumEqVector<
      Dumux::GetPropType<TypeTag, Dumux::Properties::PrimaryVariables>>;
#if DUMUX_VERSION_MAJOR >= 3 && DUMUX_VERSION_MINOR >= 9
  using SlipVelocityPolicy = Dumux::NavierStokesSlipVelocity<
      typename FVElementGeometry::GridGeometry::DiscretizationMethod,
      Dumux::NavierStokes::SlipConditions::BJ>;
  using FluxHelper = Dumux::NavierStokesMomentumBoundaryFlux<
      typename FVElementGeometry::GridGeometry::DiscretizationMethod,
      SlipVelocityPolicy>;
#else
  using FluxHelper = Dumux::NavierStokesMomentumBoundaryFluxHelper;
#endif

  ElementBoundaryTypes elemBcTypes;
  auto                 localResidual = LocalResidual(problem);
  elemBcTypes.update(*problem, element, fvGeometry);

  NumEqVector flux(0.0);
  NumEqVector neumannFlux(0.0);
  const auto &scv = fvGeometry.scv(scvf.insideScvIdx());
  for (const auto &otherScvf : scvfs(fvGeometry, scv)) {
    if (otherScvf.index() == scvf.index())
      continue;
    flux += localResidual.maybeHandleNeumannBoundary(
        *problem, element, fvGeometry, elemVolVars, elemBcTypes,
        elemFluxVarsCache, otherScvf);
  }
  flux += FluxHelper::fixedPressureMomentumFlux(
              *problem, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, 0.0,
              /*zeroNormalVelocityGradient=*/false)[scvf.normalAxis()] *
          scvf.area();
  return -1 * scvf.directionSign() * flux / scvf.area();
}

template <typename MomentumTypeTag,
          class Problem,
          class GridVariables,
          class SolutionVector>
void setInterfacePressures(const std::shared_ptr<Problem> problem,
                           const GridVariables           &gridVars,
                           const SolutionVector          &sol,
                           const std::string              meshName,
                           const std::string              dataName)
{
  const auto &gridGeometry      = problem->gridGeometry();
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
        const auto p = pressureAtInterface<MomentumTypeTag>(
            problem.get(), element, scvf, fvGeometry, elemVolVars,
            elemFluxVarsCache);
        couplingParticipant.writeScalarQuantityOnFace(
            meshName, dataName, scvf.index(), p);
      }
    }
  }
}

template <class MomentumProblem, class GridVariables, class SolutionVector>
void setInterfaceVelocities(const MomentumProblem &problem,
                            const GridVariables   &gridVars,
                            const SolutionVector  &sol,
                            const std::string      meshName,
                            const std::string      dataName)
{
  const auto &gridGeometry = problem.gridGeometry();
  auto        fvGeometry   = localView(gridGeometry);
  auto        elemVolVars  = localView(gridVars.curGridVolVars());

  auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();

  for (const auto &element : elements(gridGeometry.gridView())) {
    fvGeometry.bindElement(element);
    elemVolVars.bindElement(element, fvGeometry, sol);

    for (const auto &scvf : scvfs(fvGeometry)) {
      if (couplingParticipant.isCoupledEntity(scvf.index())) {
        // TODO: What to do here?
        const auto v =
            velocityAtInterface(elemVolVars, scvf)[scvf.normalAxis()];
        couplingParticipant.writeScalarQuantityOnFace(
            meshName, dataName, scvf.index(), v);
      }
    }
  }
}

template <class MomentumProblem, class GridVariables, class SolutionVector>
std::tuple<double, double, double> writeVelocitiesOnInterfaceToFile(
    const std::string                     &meshName,
    const std::string                     &filename,
    const std::shared_ptr<MomentumProblem> problem,
    const GridVariables                   &gridVars,
    const SolutionVector                  &sol)
{
  const auto &gridGeometry = problem->gridGeometry();
  auto        fvGeometry   = localView(gridGeometry);
  auto        elemVolVars  = localView(gridVars.curGridVolVars());

  const auto &couplingParticipant =
      Dumux::Precice::CouplingAdapter::getInstance();

  std::ofstream ofs(filename + ".csv",
                    std::ofstream::out | std::ofstream::trunc);
  ofs << "x,y,";
  if (couplingParticipant.getMeshDimensions(meshName) == 3)
    ofs << "z,";
  ofs << "velocityY"
      << "\n";

  double min = std::numeric_limits<double>::max();
  double max = std::numeric_limits<double>::min();
  double sum = 0.;
  for (const auto &element : elements(gridGeometry.gridView())) {
    fvGeometry.bind(element);
    elemVolVars.bind(element, fvGeometry, sol);

    for (const auto &scvf : scvfs(fvGeometry)) {
      if (couplingParticipant.isCoupledEntity(scvf.index())) {
        const auto &pos = scvf.center();
        for (int i = 0;
             i < couplingParticipant.getMeshDimensions(meshName); ++i) {
          ofs << pos[i] << ",";
        }
        const double v =
            velocityAtInterface(elemVolVars, scvf)[scvf.normalAxis()];
        max = std::max(v, max);
        min = std::min(v, min);
        sum += v;
        const int prec = ofs.precision();
        ofs << std::setprecision(std::numeric_limits<double>::digits10 +
                                 1)
            << v << "\n";
        ofs.precision(prec);
      }
    }
  }

  ofs.close();

  return std::make_tuple(min, max, sum);
}

template <typename MomentumTypeTag,
          class Problem,
          class GridVariables,
          class SolutionVector>
void writePressuresOnInterfaceToFile(const std::string             &meshName,
                                     const std::string             &filename,
                                     const std::shared_ptr<Problem> problem,
                                     const GridVariables           &gridVars,
                                     const SolutionVector          &sol)
{
  const auto &gridGeometry      = problem->gridGeometry();
  auto        fvGeometry        = localView(gridGeometry);
  auto        elemVolVars       = localView(gridVars.curGridVolVars());
  auto        elemFluxVarsCache = localView(gridVars.gridFluxVarsCache());

  const auto &couplingParticipant =
      Dumux::Precice::CouplingAdapter::getInstance();

  std::ofstream ofs(filename + ".csv",
                    std::ofstream::out | std::ofstream::trunc);
  ofs << "x,y,";
  if (couplingParticipant.getMeshDimensions(meshName) == 3)
    ofs << "z,";
  ofs << "pressure"
      << "\n";
  for (const auto &element : elements(gridGeometry.gridView())) {
    fvGeometry.bind(element);
    elemVolVars.bind(element, fvGeometry, sol);
    elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

    for (const auto &scvf : scvfs(fvGeometry)) {
      if (couplingParticipant.isCoupledEntity(scvf.index())) {
        const auto &pos = scvf.center();
        for (int i = 0;
             i < couplingParticipant.getMeshDimensions(meshName); ++i) {
          ofs << pos[i] << ",";
        }
        const double p = pressureAtInterface<MomentumTypeTag>(
            problem.get(), element, scvf, fvGeometry, elemVolVars,
            elemFluxVarsCache);
        ofs << p << "\n";
      }
    }
  }

  ofs.close();
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

  // Define the sub problem type tags
  using MomentumTypeTag = Properties::TTag::FreeFlowSubMomentum;
  using MassTypeTag     = Properties::TTag::FreeFlowSubMass;

  // try to create a grid (from the given grid file or the input file)
  using FreeFlowGridManager =
      Dumux::GridManager<GetPropType<MomentumTypeTag, Properties::Grid>>;
  FreeFlowGridManager freeFlowGridManager;
  freeFlowGridManager.init("FreeFlow"); // pass parameter group

  // we compute on the leaf grid view
  const auto &freeFlowGridView = freeFlowGridManager.grid().leafGridView();

  // create the finite volume grid geometry
  using MomentumGridGeometry =
      GetPropType<MomentumTypeTag, Properties::GridGeometry>;
  auto momentumGridGeometry =
      std::make_shared<MomentumGridGeometry>(freeFlowGridView);
  using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
  auto massGridGeometry =
      std::make_shared<MassGridGeometry>(freeFlowGridView);

  // create the coupling manager to couple the two subproblems of the freeflow participant
  using CouplingManager =
      GetPropType<MomentumTypeTag, Properties::CouplingManager>;
  auto           couplingManager = std::make_shared<CouplingManager>();
  constexpr auto momentumIdx     = CouplingManager::freeFlowMomentumIndex;
  constexpr auto massIdx         = CouplingManager::freeFlowMassIndex;

  // the problem (initial and boundary conditions)
  using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
  auto momentumProblem  = std::make_shared<MomentumProblem>(
      momentumGridGeometry, couplingManager);
  using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
  auto massProblem =
      std::make_shared<MassProblem>(massGridGeometry, couplingManager);

  // the solution vector
  using Traits         = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
  using SolutionVector = Traits::SolutionVector;
  SolutionVector sol;

  // Initialize preCICE.Tell preCICE about:
  std::string preciceConfigFilename = "../precice-config.xml";

  auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();
  couplingParticipant.announceConfig(mpiHelper.rank(), mpiHelper.size());

  const std::string meshName = couplingParticipant.getMeshNames()[0]; // mesh name
  const int         dim      = couplingParticipant.getMeshDimensions(meshName);
  std::cout << dim << "  " << int(MassGridGeometry::GridView::dimension)
            << std::endl;
  if (dim != int(MassGridGeometry::GridView::dimension))
    DUNE_THROW(Dune::InvalidStateException, "Dimensions do not match");

  // GET mesh corodinates
  const double xMin =
      getParamFromGroup<std::vector<double>>("Darcy", "Grid.LowerLeft")[0];
  const double xMax =
      getParamFromGroup<std::vector<double>>("Darcy", "Grid.UpperRight")[0];
  std::vector<double> coords; //( dim * vertexSize );
  std::vector<int>    coupledScvfIndices;

  for (const auto &element : elements(freeFlowGridView)) {
    auto fvGeometry = localView(*momentumGridGeometry);
    fvGeometry.bindElement(element);

    for (const auto &scvf : scvfs(fvGeometry)) {
      if (!scvf.isFrontal())
        continue;
      static constexpr auto eps = 1e-7;
      const auto           &pos = scvf.center();
      if (pos[1] < momentumGridGeometry->bBoxMin()[1] + eps) {
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

  const std::string dataNameV = couplingParticipant.getReadDataNamesOnMesh(meshName)[0];
  const std::string dataNameP = couplingParticipant.getWriteDataNamesOnMesh(meshName)[0];

  // apply initial solution for instationary problems
  momentumProblem->applyInitialSolution(sol[momentumIdx]);
  massProblem->applyInitialSolution(sol[massIdx]);

  // the grid variables
  using MomentumGridVariables =
      GetPropType<MomentumTypeTag, Properties::GridVariables>;
  auto momentumGridVariables = std::make_shared<MomentumGridVariables>(
      momentumProblem, momentumGridGeometry);
  using MassGridVariables =
      GetPropType<MassTypeTag, Properties::GridVariables>;
  auto massGridVariables =
      std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

  // initialize the coupling manager and the grid variables of the subproblems
  couplingManager->init(
      momentumProblem, massProblem,
      std::make_tuple(momentumGridVariables, massGridVariables), sol);
  momentumGridVariables->init(sol[momentumIdx]);
  massGridVariables->init(sol[massIdx]);

  bool writeInterfaceDataToFile =
      getParamFromGroup<bool>("Output", "EnableCSVWriter", false);

  // intialize the vtk output module
  using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
  VtkOutputModule freeFlowVtkWriter(*massGridVariables, sol[massIdx],
                                    massProblem->name());
  IOFields::initOutputModule(freeFlowVtkWriter);
  freeFlowVtkWriter.addVelocityOutput(
      std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
  freeFlowVtkWriter.addField(massProblem->getAnalyticalVelocityX(),
                             "analyticalV_x");
  freeFlowVtkWriter.write(0.0);

  if (couplingParticipant.requiresToWriteInitialData()) {
    setInterfacePressures<MomentumTypeTag>(
        momentumProblem, *momentumGridVariables, sol[momentumIdx], meshName,
        dataNameP);
    couplingParticipant.writeQuantityToOtherSolver(meshName, dataNameP);
  }
  couplingParticipant.initialize();
  couplingParticipant.initializeCheckpoint(sol[momentumIdx],
                                           *momentumGridVariables);
  couplingParticipant.initializeCheckpoint(sol[massIdx], *massGridVariables);

  // the assembler for a stationary problem
  using Assembler =
      MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
  auto assembler = std::make_shared<Assembler>(
      std::make_tuple(momentumProblem, massProblem),
      std::make_tuple(momentumGridGeometry, massGridGeometry),
      std::make_tuple(momentumGridVariables, massGridVariables),
      couplingManager);

  // the linear solver
  using LinearSolver =
      UMFPackIstlSolver<SeqLinearSolverTraits,
                        LinearAlgebraTraitsFromAssembler<Assembler>>;
  auto linearSolver = std::make_shared<LinearSolver>();

  // the non-linear solver
  using NewtonSolver =
      MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
  NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

  double preciceDt = couplingParticipant.getMaxTimeStepSize();
  auto   dt        = preciceDt;

  double vtkTime = 1.0;

  while (couplingParticipant.isCouplingOngoing()) {
    couplingParticipant.writeCheckpointIfRequired();

    couplingParticipant.readQuantityFromOtherSolver(meshName, dataNameV,
                                                    dt);
    // solve the non-linear system
    nonLinearSolver.solve(sol);

    if (writeInterfaceDataToFile) {
      writeVelocitiesOnInterfaceToFile(
          meshName, Dumux::Fmt::format("ff_interface_velocities_{}", vtkTime),
          momentumProblem, *momentumGridVariables, sol[momentumIdx]);
      writePressuresOnInterfaceToFile<MomentumTypeTag>(
          meshName, Dumux::Fmt::format("ff_interface_pressures_{}", vtkTime),
          momentumProblem, *momentumGridVariables, sol[momentumIdx]);
    }

    setInterfacePressures<MomentumTypeTag>(
        momentumProblem, *momentumGridVariables, sol[momentumIdx], meshName,
        dataNameP);
    couplingParticipant.writeQuantityToOtherSolver(meshName, dataNameP);
    couplingParticipant.advance(dt);
    preciceDt = couplingParticipant.getMaxTimeStepSize();
    dt        = std::min(preciceDt, dt);

    if (!couplingParticipant.readCheckpointIfRequired()) {
      vtkTime += 1.;
      freeFlowVtkWriter.write(vtkTime);
    }
  }
  ////////////////////////////////////////////////////////////
  // finalize, print dumux message to say goodbye
  ////////////////////////////////////////////////////////////

  couplingParticipant.finalize();

  // print dumux end message
  if (mpiHelper.rank() == 0) {
    Parameters::print();
    DumuxMessage::print(/*firstCall=*/false);
  }

  return 0;
}
