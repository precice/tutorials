#include <cmath>
#include "SolidSolver.h"

void SolidComputeSolution(int chunkLength, double const *const pressure, double *crossSectionLength)
{
  /*
   * Update displacement of membrane based on pressure data from the fluid solver
   */
  const double PI        = 3.14159265359;
  const double E         = 10000;
  const double r0        = 1.0 / std::sqrt(PI);
  const double c_mk      = std::sqrt(E / 2 / r0);
  const double c_mk2     = c_mk * c_mk;
  const double pressure0 = 0.0;

  for (int i = 0; i < chunkLength; i++) {
    crossSectionLength[i] = std::pow((pressure0 - 2.0 * c_mk2), 2.0) / std::pow((pressure[i] - 2.0 * c_mk2), 2.0);
  }
}
