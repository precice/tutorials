#include "FluidComputeSolution.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <vector>

using std::sin;
using std::sqrt;

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

void fluidComputeSolutionParallel(
    int     rank,
    int     size,
    int     domainSize,
    int     chunkLength,
    double  kappa,
    double  tau,
    double  gamma,
    double  scaled_t,
    double *pressure,
    double *crossSectionLength,
    double *velocity)
{
  /*
   * Step 1: Recieve the complete dataset in process 0.
   */
  if (rank != 0) {
    int    N = domainSize;
    int    i = rank;
    double chunkLength_temp;
    if ((N + 1) % size == 0) {
      chunkLength_temp = (N + 1) / size;
    } else if (i < (N + 1) % size) {
      chunkLength_temp = (N + 1) / size + 1;
    } else {
      chunkLength_temp = (N + 1) / size;
    }
    chunkLength   = chunkLength_temp;
    double *velo  = new double[chunkLength];
    double *press = new double[chunkLength];
    double *diam  = new double[chunkLength];

    int tagStart = 7 * rank; // dont needed
    MPI_Send(pressure, chunkLength, MPI_DOUBLE, 0, tagStart + 0, MPI_COMM_WORLD);
    MPI_Send(crossSectionLength, chunkLength, MPI_DOUBLE, 0, tagStart + 3, MPI_COMM_WORLD);
    MPI_Send(velocity, chunkLength, MPI_DOUBLE, 0, tagStart + 5, MPI_COMM_WORLD);

    MPI_Status status;
    MPI_Recv(press, chunkLength, MPI_DOUBLE, 0, tagStart + 0, MPI_COMM_WORLD, &status);
    MPI_Recv(diam, chunkLength, MPI_DOUBLE, 0, tagStart + 3, MPI_COMM_WORLD, &status);
    MPI_Recv(velo, chunkLength, MPI_DOUBLE, 0, tagStart + 5, MPI_COMM_WORLD, &status);

    for (int k = 0; k < chunkLength; k++) {
      velocity[k]           = velo[k];
      pressure[k]           = press[k];
      crossSectionLength[k] = diam[k];
    }

  } else {
    double *pressure_NLS;
    double *crossSectionLength_NLS;
    double *velocity_NLS;

    int N = domainSize;

    pressure_NLS = new double[N + 1];

    crossSectionLength_NLS = new double[N + 1];

    velocity_NLS = new double[N + 1];

    for (int i = 0; i < chunkLength; i++) {
      pressure_NLS[i]           = pressure[i];
      crossSectionLength_NLS[i] = crossSectionLength[i];
      velocity_NLS[i]           = velocity[i];
    }

    for (int i = 1; i < size; i++) {
      int tagStart = 7 * i;
      int chunkLength_temp;
      int gridOffset;

      if ((N + 1) % size == 0) {
        chunkLength_temp = (N + 1) / size;
        gridOffset       = i * chunkLength_temp;
      } else if (i < (N + 1) % size) {
        chunkLength_temp = (N + 1) / size + 1;
        gridOffset       = i * chunkLength_temp;
      } else {
        chunkLength_temp = (N + 1) / size;
        gridOffset       = ((N + 1) % size) * ((N + 1) / size + 1) + (i - ((N + 1) % size)) * (N + 1) / size;
      }

      MPI_Recv(pressure_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(crossSectionLength_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(velocity_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // LAPACK Variables here
    double *Res, **LHS, *A, alpha, dx, tmp, tmp2, temp_sum, norm_1, norm_2, norm = 1.0;
    int *   ipiv, nlhs = 2 * N + 2, nrhs = 1, info, ampl;

    Res = new double[2 * N + 2];
    LHS = new double *[2 * N + 2];
    for (int i = 0; i < 2 * N + 2; i++) {
      LHS[i] = new double[2 * N + 2];
    }
    A    = new double[(2 * N + 2) * (2 * N + 2)];
    ipiv = new int[2 * N + 2];

    // Stabilization intensity
    alpha = (N * kappa * tau) / (N * tau + 1);
    dx    = 1.0 / (N * kappa); //tau
    ampl  = 100;

    int whileLoopCounter = 0;
    while (1) { // Add stopping criterion
      for (int i = 0; i < 2 * N + 2; i++) {
        Res[i] = 0.0;
      }

      for (int i = 1; i < N; i++) { // chnaged to 1
        /* Momentum */
        Res[i] = velocity_NLS[i] * crossSectionLength_NLS[i] * dx;
        Res[i] = Res[i] - 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i] * velocity_NLS[i + 1] - 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i] * velocity_NLS[i + 1];
        Res[i] = Res[i] - crossSectionLength_NLS[i] * dx * velocity_NLS[i] - 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i] * velocity_NLS[i] - 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i] * velocity_NLS[i] + 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i - 1] * velocity_NLS[i] + 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i - 1] * velocity_NLS[i];
        Res[i] = Res[i] + 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i - 1] * velocity_NLS[i - 1] + 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i - 1] * velocity_NLS[i - 1];
        Res[i] = Res[i] + 0.25 * crossSectionLength_NLS[i - 1] * pressure_NLS[i - 1] + 0.25 * crossSectionLength_NLS[i] * pressure_NLS[i - 1] - 0.25 * crossSectionLength_NLS[i - 1] * pressure_NLS[i] + 0.25 * crossSectionLength_NLS[i + 1] * pressure_NLS[i] - 0.25 * crossSectionLength_NLS[i] * pressure_NLS[i + 1] - 0.25 * crossSectionLength_NLS[i + 1] * pressure_NLS[i + 1];

        /* Continuity */
        Res[i + N + 1] = -(crossSectionLength_NLS[i] - crossSectionLength_NLS[i]) * dx;
        Res[i + N + 1] = Res[i + N + 1] + 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i - 1] + 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i - 1] + 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i] - 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i] - 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i + 1] - 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i + 1];
        Res[i + N + 1] = Res[i + N + 1] + alpha * pressure_NLS[i - 1] - 2 * alpha * pressure_NLS[i] + alpha * pressure_NLS[i + 1];
      }

      // Boundary

      // Velocity
      tmp    = sin(PI * scaled_t);
      Res[0] = (1.0 / kappa) + (1.0 / (kappa * ampl)) * tmp * tmp - velocity_NLS[0];

      // Pressure Inlet is lineary interpolated
      Res[N + 1] = -pressure_NLS[0] + 2 * pressure_NLS[1] - pressure_NLS[2];

      // Velocity Outlet is lineary interpolated
      Res[N] = -velocity_NLS[N] + 2 * velocity_NLS[N - 1] - velocity_NLS[N - 2];

      // Pressure Outlet is "non-reflecting"
      tmp2           = sqrt(1 - pressure_NLS[N] / 2) - (velocity_NLS[N] - velocity_NLS[N]) / 4;
      Res[2 * N + 1] = -pressure_NLS[N] + 2 * (1 - tmp2 * tmp2);

      // Stopping Criteria
      whileLoopCounter += 1; // Iteration Count

      temp_sum = 0;
      for (int i = 0; i < (2 * N + 2); i++) {
        temp_sum += Res[i] * Res[i];
      }
      norm_1 = sqrt(temp_sum);

      temp_sum = 0;
      for (int i = 0; i < (N + 1); i++) {
        temp_sum += (pressure_NLS[i] * pressure_NLS[i]) + (velocity_NLS[i] * velocity_NLS[i]);
      }
      norm_2 = sqrt(temp_sum);

      norm = norm_1 / norm_2; // Norm

      if ((norm < 1e-15 && whileLoopCounter > 1) || whileLoopCounter > 50) {
        std::cout << "Nonlinear Solver break, Its: " << whileLoopCounter << ", norm: " << norm << std::endl;
        break;
      }

      // Initilizing the the LHS i.e. Left Hand Side
      for (int i = 0; i <= (2 * N + 1); i++) {
        for (int j = 0; j <= (2 * N + 1); j++) {
          LHS[i][j] = 0.0;
        }
      }

      for (int i = 1; i < N; i++) {
        // Momentum, Velocity
        LHS[i][i - 1] = LHS[i][i - 1] - 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i - 1] * 2 - 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i - 1] * 2 - 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i] - 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i];
        LHS[i][i]     = LHS[i][i] + 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i + 1] + 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i + 1] + crossSectionLength_NLS[i] * dx + 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i] * 2 + 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i] * 2 - 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i - 1] - 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i - 1];
        LHS[i][i + 1] = LHS[i][i + 1] + 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i] + 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i];

        // Momentum, Pressure
        LHS[i][N + 1 + i - 1] = LHS[i][N + 1 + i - 1] - 0.25 * crossSectionLength_NLS[i - 1] - 0.25 * crossSectionLength_NLS[i];
        LHS[i][N + 1 + i]     = LHS[i][N + 1 + i] + 0.25 * crossSectionLength_NLS[i - 1] - 0.25 * crossSectionLength_NLS[i + 1];
        LHS[i][N + 1 + i + 1] = LHS[i][N + 1 + i + 1] + 0.25 * crossSectionLength_NLS[i] + 0.25 * crossSectionLength_NLS[i + 1];

        // Continuity, Velocity
        LHS[i + N + 1][i - 1] = LHS[i + N + 1][i - 1] - 0.25 * crossSectionLength_NLS[i - 1] - 0.25 * crossSectionLength_NLS[i];
        LHS[i + N + 1][i]     = LHS[i + N + 1][i] - 0.25 * crossSectionLength_NLS[i - 1] + 0.25 * crossSectionLength_NLS[i + 1];
        LHS[i + N + 1][i + 1] = LHS[i + N + 1][i + 1] + 0.25 * crossSectionLength_NLS[i] + 0.25 * crossSectionLength_NLS[i + 1];

        // Continuity, Pressure
        LHS[i + N + 1][N + 1 + i - 1] = LHS[i + N + 1][N + 1 + i - 1] - alpha;
        LHS[i + N + 1][N + 1 + i]     = LHS[i + N + 1][N + 1 + i] + 2 * alpha;
        LHS[i + N + 1][N + 1 + i + 1] = LHS[i + N + 1][N + 1 + i + 1] - alpha;
      }

      // Boundary
      // Velocity Inlet is prescribed
      LHS[0][0] = 1;
      // Pressure Inlet is lineary interpolated
      LHS[N + 1][N + 1] = 1;
      LHS[N + 1][N + 2] = -2;
      LHS[N + 1][N + 3] = 1;
      // Velocity Outlet is lineary interpolated
      LHS[N][N]     = 1;
      LHS[N][N - 1] = -2;
      LHS[N][N - 2] = 1;
      // Pressure Outlet is Non-Reflecting
      LHS[2 * N + 1][2 * N + 1] = 1;
      LHS[2 * N + 1][N]         = -(sqrt(1 - pressure_NLS[N] / 2.0) - (velocity_NLS[N] - velocity_NLS[N]) / 4.0);

      // Solve Linear System using LAPACK

      /* LAPACK requires a 1D array
         i.e. Linearizing 2D
      */

      int counter = 0;
      for (int i = 0; i <= (2 * N + 1); i++) {
        for (int j = 0; j <= (2 * N + 1); j++) {
          A[counter] = LHS[j][i];
          counter++;
        }
      }

      /* LAPACK Function call to solve the linear system */
      dgesv_(&nlhs, &nrhs, A, &nlhs, ipiv, Res, &nlhs, &info);

      if (info != 0) {
        std::cout << "Linear Solver not converged!, Info: " << info << std::endl;
      }

      for (int i = 0; i < N; i++) {
        velocity_NLS[i] = velocity_NLS[i] + Res[i];
        pressure_NLS[i] = pressure_NLS[i] + Res[i + N + 1];
      }
    } // End of while loop

    for (int i = 0; i < chunkLength; i++) {
      pressure[i]           = pressure_NLS[i];
      crossSectionLength[i] = crossSectionLength_NLS[i];
      velocity[i]           = velocity_NLS[i];
    }

    for (int i = 1; i < size; i++) {
      int tagStart = 7 * i;
      int chunkLength_temp;
      int gridOffset;

      if ((N + 1) % size == 0) {
        chunkLength_temp = (N + 1) / size;
        gridOffset       = i * chunkLength_temp;
      } else if (i < (N + 1) % size) {
        chunkLength_temp = (N + 1) / size + 1;
        gridOffset       = i * chunkLength_temp;
      } else {
        chunkLength_temp = (N + 1) / size;
        gridOffset       = ((N + 1) % size) * ((N + 1) / size + 1) + (i - ((N + 1) % size)) * (N + 1) / size;
      }

      MPI_Send(pressure_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 0, MPI_COMM_WORLD);
      MPI_Send(crossSectionLength_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 3, MPI_COMM_WORLD);
      MPI_Send(velocity_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 5, MPI_COMM_WORLD);
    }

    delete[] pressure_NLS;
    delete[] crossSectionLength_NLS;
    delete[] velocity_NLS;
    delete[] Res;
    delete[] LHS;
    delete[] A;
    delete[] ipiv;
  }
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
  /* fluid_nl Variables */
  int      i, j = 0;
  double   alpha, dx;
  double   tmp1, tmp2;
  double * Res;
  double **LHS;
  double   temp_sum;
  double   norm_1, norm_2;
  double   norm = 1.0;
  double   c_mk2; //c_mk**2

  c_mk2 = 10000 / 2 * sqrt(PI);

  const int chunkLength = N + 1;
  std::copy(velocity_old, velocity_old + chunkLength, velocity);
  std::copy(pressure_old, pressure_old + chunkLength, pressure);


  // Used as Ax = b
  // i.e. LHS*x = Res
  Res = (double *) calloc((2 * N + 2), sizeof(double));
  LHS = (double **) calloc((2 * N + 2), sizeof(double *));
  for (i = 0; i < (2 * N + 2); ++i) {
    LHS[i] = (double *) calloc((2 * N + 2), sizeof(double));
  }
  /* LAPACK Variables */
  double *A = (double *) calloc((2 * N + 2) * (2 * N + 2), sizeof(double));
  ;
  int *ipiv = (int *) calloc((2 * N + 2), sizeof(int));
  ;
  int nlhs = (2 * N + 2);
  int nrhs = 1;
  int info;

  /* Stabilization Intensity */
  alpha = (N * kappa * tau) / (N * tau + 1);
  dx    = 1.0 / (N * kappa);
  //dx = 0.1;
  //alpha = 0;

  // k is the iteration counter
  for(int k = 1; ; ++k) {
    for (i = 0; i < (2 * N + 2); i++)
      Res[i] = 0.0;

    for (i = 1; i < N; i++) {
      /* Momentum */ //theta = 1
      Res[i] = (velocity_old[i] * crossSectionLength_old[i] - velocity[i] * crossSectionLength[i]) * dx / tau;
      //Res[i] = velocity_old[i] *crossSectionLength[i]* dx/tau;

      Res[i] = Res[i] - 0.25 * crossSectionLength[i + 1] * velocity[i] * velocity[i + 1] - 0.25 * crossSectionLength[i] * velocity[i] * velocity[i + 1];

      Res[i] = Res[i] - 0.25 * crossSectionLength[i + 1] * velocity[i] * velocity[i] - 0.25 * crossSectionLength[i] * velocity[i] * velocity[i] + 0.25 * crossSectionLength[i] * velocity[i - 1] * velocity[i] + 0.25 * crossSectionLength[i - 1] * velocity[i - 1] * velocity[i];

      Res[i] = Res[i] + 0.25 * crossSectionLength[i - 1] * velocity[i - 1] * velocity[i - 1] + 0.25 * crossSectionLength[i] * velocity[i - 1] * velocity[i - 1];
      Res[i] = Res[i] + 0.25 * crossSectionLength[i - 1] * pressure[i - 1] + 0.25 * crossSectionLength[i] * pressure[i - 1] - 0.25 * crossSectionLength[i - 1] * pressure[i] + 0.25 * crossSectionLength[i + 1] * pressure[i] - 0.25 * crossSectionLength[i] * pressure[i + 1] - 0.25 * crossSectionLength[i + 1] * pressure[i + 1];

      /* Continuity */
      Res[i + N + 1] = (crossSectionLength_old[i] - crossSectionLength[i]) * dx / tau;
      Res[i + N + 1] = Res[i + N + 1] + 0.25 * crossSectionLength[i - 1] * velocity[i - 1] + 0.25 * crossSectionLength[i] * velocity[i - 1] + 0.25 * crossSectionLength[i - 1] * velocity[i] - 0.25 * crossSectionLength[i + 1] * velocity[i] - 0.25 * crossSectionLength[i] * velocity[i + 1] - 0.25 * crossSectionLength[i + 1] * velocity[i + 1];
      //Res[i + N + 1] = Res[i + N + 1] + alpha * pressure[i - 1] - 2 * alpha * pressure[i] + alpha * pressure[i + 1];
      Res[i + N + 1] = Res[i + N + 1] + alpha * (pressure[i - 1] - 2 * pressure[i] + pressure[i + 1]);
    }

    /* Boundary */

    /* Velocity Inlet is prescribed */
    tmp1   = 10 + 3 * sin(10 * PI * (t)); //inlet velocity
    Res[0] = tmp1 - velocity[0];

    /* Pressure Inlet is lineary interpolated */
    Res[N + 1] = -pressure[0] + 2 * pressure[1] - pressure[2];

    /* Velocity Outlet is lineary interpolated */
    Res[N] = -velocity[N] + 2 * velocity[N - 1] - velocity[N - 2];

    /* Pressure Outlet is "non-reflecting" */
    tmp2           = sqrt(c_mk2 - pressure_old[N] / 2) - (velocity[N] - velocity_old[N]) / 4;
    Res[2 * N + 1] = -pressure[N] + 2 * (c_mk2 - tmp2 * tmp2);

    // compute norm of residual
    temp_sum = 0;
    for (i = 0; i < (2 * N + 2); i++) {
      temp_sum += Res[i] * Res[i];
    }
    norm_1   = sqrt(temp_sum);
    temp_sum = 0;
    for (i = 0; i < (N + 1); i++) {
      temp_sum += (pressure[i] * pressure[i]) + (velocity[i] * velocity[i]);
    }
    norm_2 = sqrt(temp_sum);
    norm   = norm_1 / norm_2;

    if ((norm < 1e-15 && k > 1) || k > 50) {
      printf("Nonlinear Solver break, iterations: %i, residual norm: %e\n", k, norm);
      break;
    }

    /* Initilizing the the LHS i.e. Left Hand Side */
    for (i = 0; i <= (2 * N + 1); i++)
      for (j = 0; j <= (2 * N + 1); j++)
        LHS[i][j] = 0.0;

    for (i = 1; i < N; i++) {
      // Momentum, Velocity
      LHS[i][i - 1] = LHS[i][i - 1] - 0.25 * crossSectionLength[i - 1] * velocity[i - 1] * 2 - 0.25 * crossSectionLength[i] * velocity[i - 1] * 2 - 0.25 * crossSectionLength[i] * velocity[i] - 0.25 * crossSectionLength[i - 1] * velocity[i];
      LHS[i][i]     = LHS[i][i] + 0.25 * crossSectionLength[i + 1] * velocity[i + 1] + 0.25 * crossSectionLength[i] * velocity[i + 1] + crossSectionLength[i] * dx / tau + 0.25 * crossSectionLength[i + 1] * velocity[i] * 2 + 0.25 * crossSectionLength[i] * velocity[i] * 2 - 0.25 * crossSectionLength[i] * velocity[i - 1] - 0.25 * crossSectionLength[i - 1] * velocity[i - 1];
      LHS[i][i + 1] = LHS[i][i + 1] + 0.25 * crossSectionLength[i + 1] * velocity[i] + 0.25 * crossSectionLength[i] * velocity[i];

      // Momentum, Pressure
      LHS[i][N + 1 + i - 1] = LHS[i][N + 1 + i - 1] - 0.25 * crossSectionLength[i - 1] - 0.25 * crossSectionLength[i];
      LHS[i][N + 1 + i]     = LHS[i][N + 1 + i] + 0.25 * crossSectionLength[i - 1] - 0.25 * crossSectionLength[i + 1];
      LHS[i][N + 1 + i + 1] = LHS[i][N + 1 + i + 1] + 0.25 * crossSectionLength[i] + 0.25 * crossSectionLength[i + 1];

      // Continuity, Velocity
      LHS[i + N + 1][i - 1] = LHS[i + N + 1][i - 1] - 0.25 * crossSectionLength[i - 1] - 0.25 * crossSectionLength[i];
      LHS[i + N + 1][i]     = LHS[i + N + 1][i] - 0.25 * crossSectionLength[i - 1] + 0.25 * crossSectionLength[i + 1];
      LHS[i + N + 1][i + 1] = LHS[i + N + 1][i + 1] + 0.25 * crossSectionLength[i] + 0.25 * crossSectionLength[i + 1];

      // Continuity, Pressure
      LHS[i + N + 1][N + 1 + i - 1] = LHS[i + N + 1][N + 1 + i - 1] - alpha;
      LHS[i + N + 1][N + 1 + i]     = LHS[i + N + 1][N + 1 + i] + 2 * alpha;
      LHS[i + N + 1][N + 1 + i + 1] = LHS[i + N + 1][N + 1 + i + 1] - alpha;
    }

    /* Boundary */

    // Velocity Inlet is prescribed
    LHS[0][0] = 1;
    // Pressure Inlet is lineary interpolated
    LHS[N + 1][N + 1] = 1;
    LHS[N + 1][N + 2] = -2;
    LHS[N + 1][N + 3] = 1;
    // Velocity Outlet is lineary interpolated
    LHS[N][N]     = 1;
    LHS[N][N - 1] = -2;
    LHS[N][N - 2] = 1;
    // Pressure Outlet is Non-Reflecting
    LHS[2 * N + 1][2 * N + 1] = 1;
    LHS[2 * N + 1][N]         = -(sqrt(c_mk2 - pressure_old[N] / 2.0) - (velocity[N] - velocity_old[N]) / 4.0);

    /* LAPACK requires a 1D array 
       i.e. Linearizing 2D 
    */
    int counter = 0;
    for (i = 0; i <= (2 * N + 1); i++) {
      for (j = 0; j <= (2 * N + 1); j++) {
        A[counter] = LHS[j][i];
        counter++;
      }
    }

    /* LAPACK Function call to solve the linear system */
    dgesv_(&nlhs, &nrhs, A, &nlhs, ipiv, Res, &nlhs, &info);

    if (info != 0) {
      printf("Linear Solver not converged!, Info: %i\n", info);
    }

    for (i = 0; i <= N; i++) {
      velocity[i] = velocity[i] + Res[i];
      pressure[i] = pressure[i] + Res[i + N + 1];
    }
  }
  return 0;
}
