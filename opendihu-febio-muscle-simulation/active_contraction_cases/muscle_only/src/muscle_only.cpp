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
  
  TimeSteppingScheme::RepeatedCall<
    MuscleContractionSolver< 				// 3D solid mechanics equation
      Mesh::StructuredDeformableOfDimension<3>,
      Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressibleActive3D
    >
  > problem(settings);
  
  problem.run();
  
  return EXIT_SUCCESS;
}




