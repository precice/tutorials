#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

int main(int argc, char *argv[]) {
  // multiple fibers in arbitrary partitioning, coupled to dynamic nonlinear
  // elasticity

  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  // define problem
  Control::PreciceAdapter<Control::Coupling<
      FastMonodomainSolver<           // a wrapper that improves performance of
                                      // multidomain
          Control::MultipleInstances< // fibers
              OperatorSplitting::Strang<
                  Control::MultipleInstances<
                      TimeSteppingScheme::Heun< // fiber
                                                // reaction
                                                // term
                          CellmlAdapter<
                              44, 19, // nStates,nAlgebraics: 57,71 = Shorten,
                                      // 4,9 = Hodgkin Huxley
                              FunctionSpace::FunctionSpace<
                                  Mesh::StructuredDeformableOfDimension<1>,
                                  BasisFunction::LagrangeOfOrder<1>>>>>,
                  Control::MultipleInstances<
                      TimeSteppingScheme::CrankNicolson< // fiber diffusion
                          SpatialDiscretization::FiniteElementMethod<
                              Mesh::StructuredDeformableOfDimension<1>,
                              BasisFunction::LagrangeOfOrder<1>,
                              Quadrature::Gauss<2>,
                              Equation::Dynamic::IsotropicDiffusion>>>>>>,
      MuscleContractionSolver<>>>
      problem(settings);

  // run problem
  problem.run();

  return EXIT_SUCCESS;
}
