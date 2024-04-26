#ifndef FLUID_NL_H_
#define FLUID_NL_H_

void write_vtk(double      t,
               int         iteration,
               const char *filename_prefix,
               int         N_slices,
               double     *grid,
               double     *velocity,
               double     *pressure,
               double     *diameter);

#endif
