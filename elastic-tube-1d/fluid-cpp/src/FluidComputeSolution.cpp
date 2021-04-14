#include "FluidComputeSolution.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <numeric>

using std::sin;
using std::sqrt;

// Simplifies strided access on a 1d buffer
template<class T>
class StridedAccess {
  public:
    StridedAccess(T*first, int stride) : first(first), stride(stride) {};
    // This accessor allows for row-major access
    T& operator()(int i, int j) { return first[j*stride+i]; }
  private:
    T*first;
    int stride;
};

extern "C" {
void dgesv_(
    int *   n,
    int *   nrhs,
    double *A,
    int *   lda,
    int *   ipiv,
    double *b,
    int *   ldb,
    int *   info);
}

/* Function for fluid_nl i.e. non-linear */
int fluidComputeSolutionSerial(
    double const * const velocity_old,
    double const * const pressure_old,
    double const * const crossSectionLength_old,
    double const * const crossSectionLength,
    double  t,
    int     N,
    double  kappa,
    double  tau,
    double * velocity,
    double * pressure)
{
  const double PI = 3.141592653589793;

  /* fluid_nl Variables */
  const double E = 10000;
  // c_mk^2
  const double c_mk2 = E / 2 * sqrt(PI);

  const int chunkLength = N + 1;
  std::copy(velocity_old, velocity_old + chunkLength, velocity);
  std::copy(pressure_old, pressure_old + chunkLength, pressure);

  // Used as Ax = b
  // i.e. LHS*x = Res
  std::vector<double> Res(2 * N + 2);
  std::vector<double> LHS_buffer(std::pow(2*N + 2, 2));
  StridedAccess<double> LHS(LHS_buffer.data(), 2* N + 2);

  /* LAPACK Variables */
  int nlhs = (2 * N + 2);
  int nrhs = 1;
  std::vector<int> ipiv(nlhs);

  /* Stabilization Intensity */
  const double alpha = 0.0; //(N * kappa * tau) / (N * tau + 1);
  const double L = 10.0;
  const double dx = L / kappa; // 1.0 / (N * kappa);

  // k is the iteration counter
  for(int k = 1; ; ++k) {
    std::fill(Res.begin(), Res.end(), 0.0);

    for (int i = 1; i < N; i++) {
      /* Momentum */ 
      Res[i] = (velocity_old[i] * crossSectionLength_old[i] - velocity[i] * crossSectionLength[i]) * dx / tau;

      Res[i] += 0.25 * (- crossSectionLength[i + 1] * velocity[i] * velocity[i + 1]
                        - crossSectionLength[i] * velocity[i] * velocity[i + 1]);

      Res[i] += 0.25 * (- crossSectionLength[i + 1] * velocity[i] * velocity[i]
                        - crossSectionLength[i] * velocity[i] * velocity[i]
                        + crossSectionLength[i] * velocity[i - 1] * velocity[i] 
                        + crossSectionLength[i - 1] * velocity[i - 1] * velocity[i]);

      Res[i] += 0.25 * (+ crossSectionLength[i - 1] * velocity[i - 1] * velocity[i - 1]
                        + crossSectionLength[i] * velocity[i - 1] * velocity[i - 1]);

      Res[i] += 0.25*(+ crossSectionLength[i - 1] * pressure[i - 1]
                      + crossSectionLength[i]     * pressure[i - 1]
                      - crossSectionLength[i - 1] * pressure[i]
                      + crossSectionLength[i + 1] * pressure[i]
                      - crossSectionLength[i]     * pressure[i + 1]
                      - crossSectionLength[i + 1] * pressure[i + 1]);

      /* Continuity */
      Res[i + N + 1] = (crossSectionLength_old[i] - crossSectionLength[i]) * dx / tau;
      Res[i + N + 1] += 0.25*(+ crossSectionLength[i - 1] * velocity[i - 1]
                              + crossSectionLength[i]     * velocity[i - 1]
                              + crossSectionLength[i - 1] * velocity[i]
                              - crossSectionLength[i + 1] * velocity[i]
                              - crossSectionLength[i]     * velocity[i + 1]
                              - crossSectionLength[i + 1] * velocity[i + 1]);

      Res[i + N + 1] += alpha * (pressure[i - 1] - 2 * pressure[i] + pressure[i + 1]);
    }

    /* Boundary */

    /* Velocity Inlet is prescribed */
    const double u0 = 10.0;
    const double ampl = 3.0;
    const double frequency = 10.0;
    const double t_shift = 0.0;
    const double velocity_in = u0 + ampl * sin(frequency * (t + t_shift) * PI);
    Res[0] = velocity_in - velocity[0];

    /* Pressure Inlet is lineary interpolated */
    Res[N + 1] = -pressure[0] + 2 * pressure[1] - pressure[2];

    /* Velocity Outlet is lineary interpolated */
    Res[N] = -velocity[N] + 2 * velocity[N - 1] - velocity[N - 2];

    /* Pressure Outlet is "non-reflecting" */
    const double tmp2 = sqrt(c_mk2 - pressure_old[N] / 2) - (velocity[N] - velocity_old[N]) / 4;
    Res[2 * N + 1] = -pressure[N] + 2 * (c_mk2 - std::pow(tmp2, 2));

    // compute norm of residual
    const double norm_1 = std::sqrt(
        std::inner_product(Res.begin(), Res.end(), Res.begin(), 0.0)
        );

    const double norm_2 = std::sqrt(
        std::inner_product(pressure, pressure + chunkLength, pressure, 0.0) +
        std::inner_product(velocity, velocity + chunkLength, velocity, 0.0)
        );
    const double norm = norm_1 / norm_2;

    // NOTE tolerance is 1e-10 and max iterations is 1000 in python
    if ((norm < 1e-15 && k > 1) || k > 50) {
      std::cout << "Nonlinear Solver break, iterations: " << k << ", residual norm: " << norm << '\n';
      break;
    }

    /* Initilizing the the LHS i.e. Left Hand Side */
    std::fill(LHS_buffer.begin(), LHS_buffer.end(), 0.0);

    for (int i = 1; i < N; i++) {
      // Momentum, Velocity
      LHS(i, i - 1) +=0.25*(-2 * crossSectionLength[i - 1] * velocity[i - 1]
                            -2 * crossSectionLength[i]     * velocity[i - 1]
                            -    crossSectionLength[i]     * velocity[i]
                            -    crossSectionLength[i - 1] * velocity[i]);

      LHS(i, i)     += crossSectionLength[i] * dx / tau;
      LHS(i, i)     +=  0.25 * (+    crossSectionLength[i + 1] * velocity[i + 1]
                                +    crossSectionLength[i]     * velocity[i + 1]
                                +2 * crossSectionLength[i + 1] * velocity[i] 
                                +2 * crossSectionLength[i]     * velocity[i]
                                -    crossSectionLength[i]     * velocity[i - 1]
                                -    crossSectionLength[i - 1] * velocity[i - 1]);
      LHS(i, i + 1) +=  0.25 * (+ crossSectionLength[i + 1] * velocity[i]    
                                + crossSectionLength[i]     * velocity[i]);

      // Momentum, Pressure
      LHS(i, N + 1 + i - 1) += -0.25 * crossSectionLength[i - 1] - 0.25 * crossSectionLength[i];
      LHS(i, N + 1 + i)     +=  0.25 * crossSectionLength[i - 1] - 0.25 * crossSectionLength[i + 1];
      LHS(i, N + 1 + i + 1) +=  0.25 * crossSectionLength[i]     + 0.25 * crossSectionLength[i + 1];

      // Continuity, Velocity
      LHS(i + N + 1, i - 1) += -0.25 * crossSectionLength[i - 1] - 0.25 * crossSectionLength[i];
      LHS(i + N + 1, i)     += -0.25 * crossSectionLength[i - 1] + 0.25 * crossSectionLength[i + 1];
      LHS(i + N + 1, i + 1) +=  0.25 * crossSectionLength[i]     + 0.25 * crossSectionLength[i + 1];

      // Continuity, Pressure
      LHS(i + N + 1, N + 1 + i - 1) -= alpha;
      LHS(i + N + 1, N + 1 + i)     += 2 * alpha;
      LHS(i + N + 1, N + 1 + i + 1) -= alpha;
    }

    /* Boundary */

    // Velocity Inlet is prescribed
    LHS(0, 0) = 1;
    // Pressure Inlet is lineary interpolated
    LHS(N + 1, N + 1) = 1;
    LHS(N + 1, N + 2) = -2;
    LHS(N + 1, N + 3) = 1;
    // Velocity Outlet is lineary interpolated
    LHS(N, N)     = 1;
    LHS(N, N - 1) = -2;
    LHS(N, N - 2) = 1;
    // Pressure Outlet is Non-Reflecting
    LHS(2 * N + 1, 2 * N + 1) = 1;
    LHS(2 * N + 1, N)         = -(sqrt(c_mk2 - pressure_old[N] / 2.0) - (velocity[N] - velocity_old[N]) / 4.0);

    /* LAPACK Function call to solve the linear system */
    int info{0};
    dgesv_(&nlhs, &nrhs, LHS_buffer.data(), &nlhs, ipiv.data(), Res.data(), &nlhs, &info);

    if (info != 0) {
      std::cerr << "Linear Solver not converged!, Info: " << info << '\n';
    }

    for (int i = 0; i <= N; i++) {
      velocity[i] += Res[i];
      pressure[i] += Res[i + N + 1];
    }
  }
  return 0;
}
