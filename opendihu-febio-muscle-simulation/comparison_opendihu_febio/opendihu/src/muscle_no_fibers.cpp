#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[]) {
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  TimeSteppingScheme::DynamicHyperelasticitySolver<> solver(settings);

  solver.run();

  return EXIT_SUCCESS;
}

/*
int main(int argc, char *argv[])
{
  // linear elasticity
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<3>,
    Equation::Static::LinearElasticity
  > equationDiscretized(settings);
  
  equationDiscretized.run();
  
  return EXIT_SUCCESS;
}
*/