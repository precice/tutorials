#include <pybind11/numpy.h> // numpy arrays
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // std::vector conversion

#include <config.h>

#include <ctime>
#include <fstream>
#include <iostream>
#include <numbers>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>

#include <dumux/linear/pdesolver.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/method.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include "properties_allencahn.hh"
#include "properties_cellproblem.hh"

namespace py = pybind11;
namespace Dumux {
class MicroSimulation {
  using AllenCahnTypeTag   = Dumux::Properties::TTag::AllenCahn;
  using CellProblemTypeTag = Dumux::Properties::TTag::CellProblem;
  using ACSolutionVector   = Dumux::GetPropType<AllenCahnTypeTag, Dumux::Properties::SolutionVector>;
  using CPSolutionVector   = Dumux::GetPropType<CellProblemTypeTag, Dumux::Properties::SolutionVector>;
  using ACProblem          = Dumux::GetPropType<AllenCahnTypeTag, Dumux::Properties::Problem>;
  using CPProblem          = Dumux::GetPropType<CellProblemTypeTag, Dumux::Properties::Problem>;
  using ACGridVariables    = Dumux::GetPropType<AllenCahnTypeTag, Dumux::Properties::GridVariables>;
  using CPGridVariables    = Dumux::GetPropType<CellProblemTypeTag, Dumux::Properties::GridVariables>;
  using ACAssembler        = Dumux::FVAssembler<AllenCahnTypeTag, Dumux::DiffMethod::numeric>;
  using CPAssembler        = Dumux::FVAssembler<CellProblemTypeTag, Dumux::DiffMethod::numeric>;
  using LinearSolver =
      Dumux::UMFPackIstlSolver<SeqLinearSolverTraits,
                               LinearAlgebraTraitsFromAssembler<ACAssembler>>;
  using CPLinearSolver    = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits,
                                                  LinearAlgebraTraitsFromAssembler<CPAssembler>>;
  using CPLinearPDESolver = Dumux::LinearPDESolver<CPAssembler, CPLinearSolver>;
  using ACNewtonSolver    = Dumux::NewtonSolver<ACAssembler, LinearSolver>;
  using GridGeometry      = Dumux::GetPropType<AllenCahnTypeTag, Dumux::Properties::GridGeometry>;
  using Scalar            = Dumux::GetPropType<AllenCahnTypeTag, Dumux::Properties::Scalar>;
  using GridManager       = Dumux::GridManager<Dumux::GetPropType<AllenCahnTypeTag, Dumux::Properties::Grid>>;
  using JacobianMatrix    = Dumux::GetPropType<CellProblemTypeTag, Dumux::Properties::JacobianMatrix>;
  using SolutionVector    = Dumux::GetPropType<CellProblemTypeTag, Dumux::Properties::SolutionVector>;

public:
  MicroSimulation(int simulationID);
  py::dict initialize();

  // solve takes python dict for macro_write data, dt, and returns python dict for macro_read data
  py::dict solve(py::dict macro_write_data, double dt);

  // void save_checkpoint();
  // void reload_checkpoint();

  void      setState(py::tuple phi);
  py::tuple getState() const;

private:
  const double pi_ = 3.14159265358979323846;
  double       _k_00;
  double       _k_01;
  double       _k_10;
  double       _k_11;
  double       _porosity;

  ACSolutionVector _phi;    // Solution of Allen Cahn Problem
  ACSolutionVector _phiOld; // for checkpointing
  CPSolutionVector _psi;    // Solutions(s) of Cell Problem

  // shared pointers are necessary due to partitioned nature of micro simulation
  std::shared_ptr<ACNewtonSolver>                    _acNonLinearSolver;
  std::shared_ptr<LinearSolver>                      _acLinearSolver;
  std::shared_ptr<CPLinearSolver>                    _cpLinearSolver;
  std::shared_ptr<CPLinearPDESolver>                 _cpLinearPDESolver;
  std::shared_ptr<ACAssembler>                       _acAssembler;
  std::shared_ptr<CPAssembler>                       _cpAssembler;
  std::shared_ptr<Dumux::CheckPointTimeLoop<double>> _timeLoop;
  std::shared_ptr<ACProblem>                         _acProblem;
  std::shared_ptr<CPProblem>                         _cpProblem;
  std::shared_ptr<CPGridVariables>                   _cpGridVariables;
  std::shared_ptr<ACGridVariables>                   _acGridVariables;
  std::shared_ptr<GridGeometry>                      _gridGeometry;
  GridManager                                        _gridManager;
};

// Constructor
MicroSimulation::MicroSimulation(int simulationID)
{
  using namespace Dumux;

  std::cout << "Initialize micro problem \n";

  // parse the input file
  Parameters::init("params.input");

  // try to create a grid (from the given grid file or the input file)
  _gridManager.init();

  // we compute on the leaf grid view
  const auto &leafGridView = _gridManager.grid().leafGridView();

  // create the finite volume grid geometry
  _gridGeometry = std::make_shared<GridGeometry>(leafGridView);
  _gridGeometry->update(leafGridView);

  // get some time loop parameters
  const auto tEnd  = getParam<Scalar>("TimeLoop.TEnd");
  const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
  auto       dt    = getParam<Scalar>("TimeLoop.DtInitial");

  // instantiate time loop
  _timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
  _timeLoop->setMaxTimeStepSize(maxDt);

  ////////////////////////////////////
  // Set up the Allen-Cahn Problem  //
  ////////////////////////////////////

  // the AC problem
  _acProblem = std::make_shared<ACProblem>(_gridGeometry);

  // the AC solution vector
  _acProblem->applyInitialSolution(_phi);
  _phiOld = _phi;

  // the AC grid variables
  _acGridVariables = std::make_shared<ACGridVariables>(_acProblem, _gridGeometry);
  _acGridVariables->init(_phi);

  // the AC assembler with time loop for the instationary problem
  _acAssembler = std::make_shared<ACAssembler>(_acProblem, _gridGeometry, _acGridVariables, _timeLoop, _phiOld);

  // the non-linear (Newton) solver based on the linear solver for the AC problem
  _acLinearSolver    = std::make_shared<LinearSolver>();
  _acNonLinearSolver = std::make_shared<ACNewtonSolver>(_acAssembler, _acLinearSolver);

  ////////////////////////////////////////////////////////////
  // Set up the Cell Problem
  ////////////////////////////////////////////////////////////

  // the cell problem
  _cpProblem = std::make_shared<CPProblem>(_gridGeometry);

  // the CP solution vector
  CPSolutionVector psi(_gridGeometry->numDofs());
  _psi = psi;

  // the CP grid variables
  _cpGridVariables = std::make_shared<CPGridVariables>(_cpProblem, _gridGeometry);
  _cpGridVariables->init(_psi);

  // the CP assembler for the stationary problem
  _cpAssembler = std::make_shared<CPAssembler>(_cpProblem, _gridGeometry, _cpGridVariables);

  // the CP linear solver for the CP problem
  _cpLinearSolver    = std::make_shared<CPLinearSolver>();
  _cpLinearPDESolver = std::make_shared<CPLinearPDESolver>(_cpAssembler, _cpLinearSolver);

  // start tracking time
  _timeLoop->start();
};

// Initialize micro-data to be used in initial adaptivity
py::dict MicroSimulation::initialize()
{
  // update Phi in the cell problem
  _cpProblem->spatialParams().updatePhi(_phi);

  // solve the cell problems
  _cpLinearPDESolver->solve(_psi);

  // calculate porosity
  _porosity = _acProblem->calculatePorosity(_phi);

  // compute the psi derivatives (required for conductivity tensor)
  _cpProblem->computePsiDerivatives(*_cpProblem, *_cpAssembler, *_cpGridVariables, _psi);

  // calculate the conductivity tensor
  _k_00 = _cpProblem->calculateConductivityTensorComponent(0, 0);
  _k_11 = _cpProblem->calculateConductivityTensorComponent(1, 1);

  // create python dict for micro_write_data
  py::dict micro_write_data;

  // add micro_scalar_data and micro_vector_data to micro_write_data
  micro_write_data["k_00"]     = _k_00;
  micro_write_data["k_11"]     = _k_11;
  micro_write_data["porosity"] = _porosity;

  return micro_write_data;
}

// Solve
py::dict MicroSimulation::solve(py::dict macro_write_data, double dt)
{
  // call leafgridView and point gridGeometry to it
  const auto &leafGridView = _gridManager.grid().leafGridView();
  _gridGeometry->update(leafGridView);

  std::cout << "Solve timestep of micro problem \n";

  // assert(dt != 0);
  if (dt == 0) {
    std::cout << "dt is zero\n";
    exit(1);
  }

  _timeLoop->setTimeStepSize(dt);

  // read concentration from preCICE
  double conc = macro_write_data["concentration"].cast<double>();

  // input macro concentration into allen-cahn problem
  _acProblem->updateConcentration(conc);

  // linearize & solve the allen cahn problem
  _acNonLinearSolver->solve(_phi, *_timeLoop);

  // u pdate Phi in the cell problem
  _cpProblem->spatialParams().updatePhi(_phi);

  // solve the cell problems
  _cpLinearPDESolver->solve(_psi);

  std::cout << "Compute upscaled quantities \n";

  // calculate porosity
  _porosity = _acProblem->calculatePorosity(_phi);

  // compute the psi derivatives (required for conductivity tensor)
  _cpProblem->computePsiDerivatives(*_cpProblem, *_cpAssembler, *_cpGridVariables, _psi);

  // calculate the conductivity tensor
  _k_00 = _cpProblem->calculateConductivityTensorComponent(0, 0);
  _k_10 = _cpProblem->calculateConductivityTensorComponent(1, 0);
  _k_01 = _cpProblem->calculateConductivityTensorComponent(0, 1);
  _k_11 = _cpProblem->calculateConductivityTensorComponent(1, 1);

  // create python dict for micro_write_data
  py::dict micro_write_data;

  // add micro_scalar_data and micro_vector_data to micro_write_data
  micro_write_data["k_00"]       = _k_00;
  micro_write_data["k_10"]       = _k_10;
  micro_write_data["k_01"]       = _k_01;
  micro_write_data["k_11"]       = _k_11;
  micro_write_data["porosity"]   = _porosity;
  micro_write_data["grain_size"] = std::sqrt((1 - _porosity) / pi_);

  // write current primary variables to previous primary variables
  _acGridVariables->advanceTimeStep();

  // return micro_write_data
  return micro_write_data;
}

// This function needs to set the complete state of a micro simulation
void MicroSimulation::setState(py::tuple phi)
{
  py::list phi_py    = phi[0];
  py::list phiOld_py = phi[1];

  for (int i = 0; i < py::len(phi_py); ++i) {
    _phi[i]    = phi_py[i].cast<double>();
    _phiOld[i] = phiOld_py[i].cast<double>();
  }
  _acGridVariables->update(_phiOld);
  _acGridVariables->advanceTimeStep();
  _acGridVariables->update(_phi);
}

// This function needs to return variables which can fully define the state of a micro simulation
py::tuple MicroSimulation::getState() const
{
  py::list phi_py;
  py::list phiOld_py;

  for (const auto &x : this->_phi) {
    phi_py.append(x[0]);
  }

  for (const auto &xOld : this->_phiOld) {
    phiOld_py.append(xOld[0]);
  }

  return py::make_tuple(phi_py, phiOld_py);
}

PYBIND11_MODULE(micro_sim, m)
{
  m.doc() = "pybind11 example plugin"; // optional module docstring

  py::class_<MicroSimulation>(m, "MicroSimulation")
      .def(py::init<int>())
      .def("initialize", &MicroSimulation::initialize)
      .def("solve", &MicroSimulation::solve)
      //.def("save_checkpoint", &MicroSimulation::save_checkpoint)
      //.def("reload_checkpoint", &MicroSimulation::reload_checkpoint)
      .def("get_state", &MicroSimulation::getState)
      .def("set_state", &MicroSimulation::setState)
      .def(py::pickle(
          [](const MicroSimulation &ms) { // __getstate__
            /* Return a tuple that fully encodes the state of the object */
            return ms.getState();
          },
          [](py::tuple t) { // __setstate__
            if (t.size() != 2)
              throw std::runtime_error("Invalid state!");

            /* Create a new C++ instance */
            MicroSimulation ms(0);
            ms.initialize();
            ms.setState(t);

            return ms;
          }));
}
} //  end namespace Dumux
