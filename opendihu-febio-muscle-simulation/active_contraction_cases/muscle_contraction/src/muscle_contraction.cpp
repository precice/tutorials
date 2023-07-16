#include <Python.h>
#include <iostream>
#include <cstdlib>
#include <iostream>
#include "easylogging++.h"
#include "opendihu.h"

int main(int argc, char *argv[])
{
  // initialize OpenDihu 
  DihuContext settings(argc, argv);
  
  // define problem with FastMonodomainSolver and MuscleContractionSolver
  Control::Coupling<					// Couple fibers (FastMonodomain) with solid mechanics (MuscleContraction)
    FastMonodomainSolver<                
      Control::MultipleInstances<                       // subdomains in xy-plane
        OperatorSplitting::Strang<
          Control::MultipleInstances<			// fiber reaction term
            TimeSteppingScheme::Heun<                   
              CellmlAdapter<
                9, 19,  					// nStates, nAlgebraics
                FunctionSpace::FunctionSpace<
                  Mesh::StructuredDeformableOfDimension<1>,
                  BasisFunction::LagrangeOfOrder<1>
                >
              >
            >
          >,
          Control::MultipleInstances<			// fiber diffusion
            TimeSteppingScheme::ImplicitEuler<          // note that implicit euler gives lower error in this case than crank nicolson
              SpatialDiscretization::FiniteElementMethod<
                Mesh::StructuredDeformableOfDimension<1>,
                BasisFunction::LagrangeOfOrder<1>,
                Quadrature::Gauss<2>,
                Equation::Dynamic::IsotropicDiffusion
              >
            >
          >
        >
      >
    >,
    MuscleContractionSolver< 				// 3D solid mechanics equation
      Mesh::StructuredDeformableOfDimension<3>,
      Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressibleActive3D
    >
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}




