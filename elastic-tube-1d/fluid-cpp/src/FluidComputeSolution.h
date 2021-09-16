#pragma once

const double PI = 3.14159265359;

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
    double * pressure);
