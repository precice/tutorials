#include "utilities.h"
#include <fstream>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

/* 
   Function for solving the linear system 
   LAPACK is used DGESV computes the solution to a real system of linear equations
   A * x = b,
   where A is an N-by-N matrix and x and b are N-by-NRHS matrices.
*/
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

void initializeWriting(std::ofstream &filestream)
{
  filestream.setf(std::ios::showpoint);
  filestream.setf(std::ios::scientific);
  filestream << std::setprecision(16);
}

void writeHeader(std::ostream &outFile)
{
  outFile << "# vtk DataFile Version 2.0" << std::endl
          << std::endl
          << "ASCII" << std::endl
          << std::endl
          << "DATASET UNSTRUCTURED_GRID" << std::endl
          << std::endl;
}

void exportMesh(std::ofstream &outFile, int N_slices, double *grid)
{
  // Plot vertices
  outFile << "POINTS " << N_slices << " float " << std::endl
          << std::endl;

  for (int i = 0; i < N_slices; i++) {
    // read x,y from grid. Set z = 0
    // Values are stored in grid in the following way: [x_0,y_0,x_1,y_1,...x_n-1,y_n-1]
    double x = grid[2 * i + 0];
    double y = grid[2 * i + 1];
    double z = 0.0;
    outFile << x << "  " << y << "  " << z << std::endl;
  }
  outFile << std::endl;
}

void exportVectorData(std::ofstream &outFile, int N_slices, double *data, const char *dataname)
{
  outFile << "VECTORS " << dataname << " float" << std::endl;

  for (int i = 0; i < N_slices; i++) {
    // Plot vertex data
    // read x vector component from dataset. Set y,z = 0
    // Values are stored in dataset in the following way: [vx_0,vx_1,...vx_n-1]
    double vx = data[i];
    double vy = 0.0;
    double vz = 0.0;
    outFile << vx << "  " << vy << "  " << vz << std::endl;
  }

  outFile << std::endl;
}

void exportScalarData(std::ofstream &outFile, int N_slices, double *data, std::string dataname)
{
  outFile << "SCALARS " << dataname << " float" << std::endl;
  outFile << "LOOKUP_TABLE default" << std::endl;

  for (int i = 0; i < N_slices; i++) {
    // Plot vertex data
    outFile << data[i] << std::endl;
  }

  outFile << std::endl;
}

void write_vtk(double t, int iteration, const char *filename_prefix, int N_slices, double *grid, double *velocity, double *pressure, double *diameter)
{
  std::stringstream filename_stream;
  filename_stream << filename_prefix << "_" << iteration << ".vtk";
  std::string filename = filename_stream.str();
  printf("writing timestep at t=%f to %s\n", t, filename.c_str());

  std::ofstream outstream(filename);

  initializeWriting(outstream);
  writeHeader(outstream);
  exportMesh(outstream, N_slices, grid);

  outstream << "POINT_DATA " << N_slices << std::endl;
  outstream << std::endl;

  exportVectorData(outstream, N_slices, velocity, "velocity");
  exportScalarData(outstream, N_slices, pressure, "pressure");
  exportScalarData(outstream, N_slices, diameter, "diameter");

  outstream.close();
}
