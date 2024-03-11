#include <Python.h>
#include <cstdlib>
#include <iostream>
#include <easylogging++.h>
#include <opendihu.h>

int main(int argc, char *argv[])
{

  DihuContext settings(argc, argv);

  Control::PreciceAdapter<TimeSteppingScheme::DynamicHyperelasticitySolver<
      Equation::SolidMechanics::HyperelasticTendon>>
      problem(settings);

  problem.run();

  return EXIT_SUCCESS;
}
