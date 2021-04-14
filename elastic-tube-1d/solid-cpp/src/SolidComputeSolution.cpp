#include <cmath>
#include "SolidSolver.h"

void SolidComputeSolution(int chunkLength, double const * const pressure, double *crossSectionLength)
{
  /*
   * Update displacement of membrane based on pressure data from the fluid solver
   */
  const double PI    = 3.14159265359;
  const double E = 10000;
  const double c_mk2 = E / 2 * std::sqrt(PI);
  const double pressure0 = 0.0;
  for (int i = 0; i < chunkLength; i++) {
    //crossSectionLength[i] = 4.0 / ((2.0 - pressure[i]) * (2.0 - pressure[i]));
    crossSectionLength[i] = std::pow((pressure0 - 2.0 * c_mk2) / (pressure[i] - 2.0 * c_mk2), 2);
  }
}
