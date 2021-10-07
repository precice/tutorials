// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <stdio.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dune/istl/matrixmarket.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bdmatrix.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/elastodynamics/assemblers/operatorassembler.hh>
#include <dune/elastodynamics/assemblers/stiffnessassembler.hh>
#include <dune/elastodynamics/assemblers/hrzlumpedmassassembler.hh>

#include <dune/elastodynamics/timesteppers/coefficients.hh>
#include <dune/elastodynamics/timesteppers/timestepcontroller.hh>
#include <dune/elastodynamics/timesteppers/rungekuttanystroem.hh>

#include <dune/precice/couplinginterface.hh>

using namespace Dune;

struct couplingParameters {

  const std::string config_file      = "../precice-config.xml";
  const std::string participant_name = "Solid";
  const std::string mesh_name        = "Solid-Mesh";
  const std::string read_data_name   = "Force";
  const std::string write_data_name  = "Displacement";

};

int main(int argc, char** argv) {

  // set up MPI
  const MPIHelper& mpiHelper = MPIHelper::instance(argc, argv);
  int mpiRank = mpiHelper.rank();
  int mpiSize = mpiHelper.size();
  
  const int dim = 2;
  const int p = 2;

  // load and generate grid
  using Grid = UGGrid<dim>;
  using GridView = Grid::LeafGridView;

  GridFactory<Grid> factory;
  GmshReader<Grid>::read(factory, "Solid.msh");
  std::shared_ptr<Grid> grid(factory.createGrid());    

  // partition grid on different processors
  GridView gridView = grid->leafGridView();
  std::cout << "\nThere are " << grid->leafGridView().comm().size() << " processes active!" << std::endl;
 
  // generate Basis
  using namespace Functions::BasisFactory;
  auto basis = makeBasis(gridView, power<dim>(lagrange<p>()));
  using Basis = decltype(basis);

  // define operators needed
  using operatorType = BCRSMatrix<FieldMatrix<double, dim, dim>>;
  using diagonalType = BDMatrix<FieldMatrix<double, dim, dim>>;
  using blockVector  = BlockVector<FieldVector<double, dim>>;

  // assemble problem
  Elastodynamics::OperatorAssembler<Basis> operatorAssembler(basis);
  double E = 4000000.0, nu = 0.3;
  operatorType stiffnessMatrix;
  operatorAssembler.initialize(stiffnessMatrix);
  Elastodynamics::StiffnessAssembler stiffnessAssembler(E, nu);
  operatorAssembler.assemble(stiffnessAssembler, stiffnessMatrix, false);
  
  double rho = 3000.0;
  diagonalType massMatrix(basis.size());
  Elastodynamics::HRZLumpedMassAssembler massAssembler(rho);
  operatorAssembler.assemble(massAssembler, massMatrix, true);
  massMatrix.invert();

  blockVector loadVector(basis.size()), displacementVector(basis.size()), velocityVector(basis.size()), accelerationVector(basis.size());
  loadVector = 0.0, displacementVector = 0.0, velocityVector = 0.0, accelerationVector = 0.0;

  // state quantatiy vector
  std::vector<blockVector*> stateQuantaties;
  stateQuantaties = {&displacementVector, &velocityVector, &accelerationVector};

  // set dirichlet boundary
  BitSetVector<dim> dirichletDofs(basis.size(), false);
  auto dirichletPredicate = [](auto position) {return position[1] < 0.00001;};  
  Functions::interpolate(basis, dirichletDofs, dirichletPredicate);
 
  DiagonalMatrix<double, dim> I(1.0);
  FieldMatrix<double, dim> OM(0.0);
  
  // set fixed Dirichlet entries in stiffnessmatrix
  for( int i=0; i<stiffnessMatrix.N(); i++) {
    if(dirichletDofs[i][0]) {
      auto cIt = stiffnessMatrix[i].begin();
      auto cEndIt = stiffnessMatrix[i].end();  
      for( ; cIt!=cEndIt; ++cIt) {
        *cIt = (i==cIt.index()) ? I : OM;
      }
    }
  }
  
  // set fixed Dirichlet entries in massmatrix
  for( int i=0; i<massMatrix.N(); i++) {
    if(dirichletDofs[i][0]) {
      auto cIt = massMatrix[i].begin();
      auto cEndIt = massMatrix[i].end();  
      for( ; cIt!=cEndIt; ++cIt) {
        *cIt = (i==cIt.index()) ? I : OM;
      }
    }
  }

  FieldVector<double, dim> OV(0.0);

  // set fixed Dirichlet entries in load vector
  for( int i=0; i<loadVector.N(); i++) {
    if(dirichletDofs[i][0]) {
      loadVector[i] = OV;
    }
  }

  // set neumann boundary
  BlockVector<FieldVector<bool, dim>> neumannDofs(basis.size());
  neumannDofs = false;
  auto neumannPredicate = [](auto position) {return position[1] > -3.000001 || position[0] < 2.955555 || position[0] > 3.055555;};  
  
  Functions::interpolate(basis, neumannDofs, neumannPredicate);
  
  // generate global displacement function
  using displacementRange = Dune::FieldVector<double, dim>;
  auto displacementFunction = Functions::makeDiscreteGlobalBasisFunction<displacementRange> (basis, displacementVector);
  
  // generate output writer
  double dt_out = 0.1;
  auto vtkWriter = std::make_shared<SubsamplingVTKWriter<GridView>> (gridView, refinementLevels(0));
  VTKSequenceWriter<GridView> vtkSequenceWriter(vtkWriter, "solid");
  vtkWriter->addVertexData(displacementFunction, VTK::FieldInfo("displacement", VTK::FieldInfo::Type::vector, dim));
  vtkSequenceWriter.write(0.0);

  // setup looping
  double t = 0.0;
  double dt = 0.0001;
  int iteration_count = 0;

  // set up RKN method
  FixedStepController fixed(t, dt);
  RKNCoefficients coefficients = RKN5();
  RungeKuttaNystroem<operatorType, blockVector> rkn(massMatrix, stiffnessMatrix, coefficients, fixed);
  rkn.initialize(loadVector);

  // set up coupling interface
  couplingParameters parameters; 
  preCICE::CouplingInterface<dim, blockVector, couplingParameters> couplingInterface(parameters, neumannDofs, mpiRank, mpiSize); 
  couplingInterface.initialize(basis);
  
  while(couplingInterface.is_coupling_ongoing()) {
  
    if(couplingInterface.is_save_required()) {
        couplingInterface.save_current_state(stateQuantaties, t, iteration_count);
        couplingInterface.mark_save_fullfilled();
    }
    
    // get neumann values from fluid solver
    couplingInterface.read_blockvector_data(loadVector);
 
    // setup fixed dirichlet values
    for( int i=0; i<loadVector.N(); i++) {
      if(dirichletDofs[i][0]) {
        loadVector[i] = OV;
      }
    }

    rkn.step(displacementVector, velocityVector, accelerationVector, loadVector);
            
    dt = std::min(couplingInterface.precice_timestep_length, dt);
                          
    couplingInterface.send_blockvector_data(displacementVector);
    couplingInterface.advance(dt);
                     
    if(couplingInterface.is_load_required()) {
        couplingInterface.reload_old_state(stateQuantaties, t, iteration_count);
        couplingInterface.mark_load_fullfilled();
    }
    else {
      t += dt;
      ++iteration_count;
      std::cout << "=== iteration: " << iteration_count << " time: " << t << std::endl;
    }
    
    if(couplingInterface.is_time_window_complete()) {
      double tol = 10e-5;
      if(std::abs(std::fmod((t+tol), dt_out)) < 2.0*tol) {
        std::cout << "output vtk for time = " << t << "\n" << std::endl;
        vtkSequenceWriter.write(t);
      }
    }
  }
  
  couplingInterface.finalize();

}
