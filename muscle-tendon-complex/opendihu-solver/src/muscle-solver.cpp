#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

int main(int argc, char *argv[]) {
 
  DihuContext settings(argc, argv);

  Control::PreciceAdapter<Control::Coupling<
      FastMonodomainSolver<           
          Control::MultipleInstances< 
              OperatorSplitting::Strang<
                  Control::MultipleInstances<
                      TimeSteppingScheme::Heun< 
                          CellmlAdapter<
                              44, 19, 
                              FunctionSpace::FunctionSpace<
                                  Mesh::StructuredDeformableOfDimension<1>,
                                  BasisFunction::LagrangeOfOrder<1>>>>>,
                  Control::MultipleInstances<
                      TimeSteppingScheme::CrankNicolson< 
                          SpatialDiscretization::FiniteElementMethod<
                              Mesh::StructuredDeformableOfDimension<1>,
                              BasisFunction::LagrangeOfOrder<1>,
                              Quadrature::Gauss<2>,
                              Equation::Dynamic::IsotropicDiffusion>>>>>>,
      MuscleContractionSolver<>>>
      problem(settings);

  problem.run();

  return EXIT_SUCCESS;
}
