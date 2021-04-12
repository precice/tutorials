#include <cmath>
#include "SolidSolver.h"

void SolidComputeSolution(int rank, int size, int chunkLength, double *pressure, double *crossSectionLength)
{
  /*
   * Update displacement of membrane based on pressure data from the fluid solver
   */
  const double PI    = 3.14159265359;
  double       c_mk2 = 10000 / 2 * sqrt(PI);
  for (int i = 0; i < chunkLength; i++) {
    //crossSectionLength[i] = 4.0 / ((2.0 - pressure[i]) * (2.0 - pressure[i]));
    crossSectionLength[i] = (0.0 - 2.0 * c_mk2) * (0.0 - 2.0 * c_mk2) / (pressure[i] - 2.0 * c_mk2) / (pressure[i] - 2.0 * c_mk2);
    ;
  }
}
