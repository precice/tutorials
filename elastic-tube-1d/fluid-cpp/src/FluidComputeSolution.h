#pragma once

const double PI = 3.14159265359;

void fluidComputeSolutionParallel(
    int     rank,
    int     size,
    int     domainSize,
    int     chunkLength,
    double  kappa,
    double  tau,
    double  gamma,
    double  t,
    double *pressure,
    double *crossSectionLength,
    double *velocity);

int fluidComputeSolutionSerial(double *crossSectionLength,
                               double *crossSectionLength_old,
                               double *velocity,
                               double *pressure,
                               double  t,
                               int     N,
                               double  kappa,
                               double  tau);
