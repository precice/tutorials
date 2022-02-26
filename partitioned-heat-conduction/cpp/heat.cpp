#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <precice/SolverInterface.hpp>

enum struct Problem { Dirichlet, Neumann };

[[noreturn]] void invalidUse() {
  std::cerr << "Usage: heat SOLVER TOL [NUM]\n";
  std::cerr << "       SOLVER - N for Neumann or D for Dirichlet\n";
  std::cerr << "       TOL    - the tolerance for the Gauss-Seidel solver\n";
  std::cerr << "       NUM    - the number of grid points in each direction\n";
  std::exit(1);
}

double u(double x, double y, double t) {
  double alpha = 3;
  double beta = 1.6;
  return 1 + x * x + alpha * y * y + beta * t;
}

enum struct Boundary { Top, Left, Bottom, Right };

struct Domain {
  double top, left, bottom, right;

  double width() const { return right - left; }
  double height() const { return top - bottom; }
};

Domain problemDomain(Problem problem) {
  Domain d;
  d.top = 1.;
  d.bottom = 0.;
  d.left = (problem == Problem::Dirichlet) ? 0. : 1.;
  d.right = (problem == Problem::Dirichlet) ? 1. : 2.;
  return d;
}

std::string to_padded_string(int n) {
  std::string num = std::to_string(n);
  const int pad = 4 - num.length();
  return (pad <= 0) ? num : num.insert(0, pad, '0');
}

struct Grid {
  Grid(int width, int height) : width(width), height(height) {
    data.resize(width * height);
  }

  template <typename Func>
  Grid(int width, int height, Func f) : Grid(width, height) {
    evaluate(std::move(f));
  }

  Grid(Grid const&) = default;

  Grid& operator=(Grid const& other) {
    assert(this != &other);
    assert(width == other.width);
    assert(height == other.height);
    std::copy(other.data.begin(), other.data.end(), data.begin());
    return *this;
  }

  int idx(int x, int y) const { return y * width + x; }

  double &at(int x, int y) {
    assert(x >= 0 && x < width);
    assert(y >= 0 && y < height);
    return data[idx(x, y)];
  }

  double at(int x, int y) const {
    assert(x >= 0 && x < width);
    assert(y >= 0 && y < height);
    return data[idx(x, y)];
  }

  template <typename Func> void evaluate(Func f) {
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        at(x, y) = f(x, y);
      }
    }
  }

  template <typename Func> void evaluate(Func f, Boundary b) {
    switch (b) {
    case Boundary::Top:
      for (int x = 0; x < width; ++x) {
        at(x, 0) = f(x, 0);
      }
      break;
    case Boundary::Left:
      for (int y = 0; y < height; ++y) {
        at(0, y) = f(0, y);
      }
      break;
    case Boundary::Bottom:
      for (int x = 0; x < width; ++x) {
        at(x, height - 1) = f(x, height - 1);
      }
      break;
    case Boundary::Right:
      for (int y = 0; y < height; ++y) {
        at(width - 1, y) = f(width - 1, y);
      }
      break;
    }
  }

  int size() const { return width * height; }

  void writeVTK(std::string const &name, Domain d, int number) const {

    std::ofstream file{name + "-" + to_padded_string(number) + ".vtk"};

    auto dx = d.width() / (width - 1);
    auto dy = d.height() / (height - 1);

    file << "# vtk DataFile Version 2.0\n";
    file << "Export number " << number << '\n';
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << width << " " << height << " 1\n";
    file << "ORIGIN " << d.left << " " << d.bottom << " 0\n";
    file << "SPACING " << dx << " " << dy << " 1\n\n";
    file << "POINT_DATA " << size() << '\n';
    file << "SCALARS T DOUBLE 1\n";
    file << "LOOKUP_TABLE default\n";

    file << std::fixed;
    for (double s : data) {
      file << s << '\n';
    }
  }

  template <typename Iter> void writeInnerColTo(int col, Iter iter) {
    assert(col >= 0 && col < width);
    for (int y = 1; y < height - 1; ++y) {
      *iter = at(col, y);
      ++iter;
    }
  }

  template <typename Iter> void readInnerColFrom(int col, Iter iter) {
    assert(col >= 0 && col < width);
    for (int y = 1; y < height - 1; ++y) {
      at(col, y) = *iter;
      ++iter;
    }
  }

  int width, height;
  std::vector<double> data;
};

std::ostream &operator<<(std::ostream &out, const Grid &g) {
  out << std::fixed;
  for (int y = 0; y < g.height; ++y) {
    for (int x = 0; x < g.width; ++x) {
      out << g.at(x, y) << " ";
    }
    out << '\n';
  }
  out << std::defaultfloat;
  return out;
}

template <typename Iter> void computeFlux(Grid const &T, double dx, Iter iter) {
  const double dx2 = dx * dx;
  std::array<double, 3> c{1.0 / dx2, -2.0 / dx2, 1.0 / dx2};
  const int x = T.width - 1;
  for (int y = 1; y < T.height - 1; ++y) {
    // backward finite difference 2nd order
    *iter = T.at(x - 2, y) * c[0] + T.at(x - 1, y) * c[1] + T.at(x, y) * c[2];
    ++iter;
  }
}

struct DoNothing {
  template <typename... T> void operator()(T...) const {};
};

// Uses the given flux to compute the left boundary
struct ApplyFlux {
  double *fluxBegin;
  double dx;

  void operator()(Grid &T) const {
    auto flux = fluxBegin;
    for (int y = 1; y < T.height - 1; ++y, ++flux) {
      T.at(0, y) = T.at(1, y) - dx * *flux;
    }
  }
};

/// Compute an implicit timestep using GaussSeidel
/// The preTimestep function will be applied to the result each iteration.
template <class Func = DoNothing>
void timestep(Grid &T, Grid const &T_last, double dx, double dy, double dt,
              double tol, Func preTimestep = DoNothing{}) {
  T = T_last;

  assert(T.size() == T_last.size());
  double dx2 = dx * dx;
  double dy2 = dy * dy;

  double xfac = -dt / dx2;
  double yfac = -dt / dy2;
  double dfac = 1 - 2 * xfac - 2 * yfac;

  constexpr int max_iter{100000};
  int iteration = 1;
  double residual = 0.0;
  while (true) {
    if (iteration > max_iter) {
      std::cerr << "FATAL ERROR: Reached maxium amount of iterations ("
                << max_iter << ")\n"
                << "Last residual: " << residual << " (tol: " << tol << ")\n";
      std::exit(1);
    }

    preTimestep(T);

    for (int y = 1; y < T.height - 1; ++y) {
      for (int x = 1; x < T.width - 1; ++x) {
        T.at(x, y) =
            (T_last.at(x, y) - xfac * (T.at(x - 1, y) + T.at(x + 1, y)) -
             yfac * (T.at(x, y - 1) + T.at(x, y + 1))) /
            dfac;
      }
    }

    residual = 0.0;
    for (int y = 1; y < T.height - 1; ++y) {
      for (int x = 1; x < T.width - 1; ++x) {
        residual += std::pow(T_last.at(x, y) - dfac * T.at(x, y) -
                                 xfac * (T.at(x - 1, y) + T.at(x + 1, y)) -
                                 yfac * (T.at(x, y - 1) + T.at(x, y + 1)),
                             2);
      }
    }
    residual = std::sqrt(residual / ((T.width - 2) * (T.height - 2)));

    if (residual < tol) {
      std::cout << "Converged after " << iteration << " iterations\n";
      return;
    }

    ++iteration;
  }
}

auto defineInterfaceMesh(precice::SolverInterface &interface, int meshID,
                         int ny) {
  std::array<double, 2> coords{1, 0};
  std::vector<int> vids;
  vids.reserve(ny - 2);
  double dy = 1.0 / ny;
  for (int y = 1; y < ny - 1; ++y) {
    coords[1] = y * dy;
    vids.push_back(interface.setMeshVertex(meshID, coords.data()));
  }
  return vids;
}

struct Checkpoint {
  Grid T;
  double t;
};

int main(int argc, char *argv[]) {
  if (!(argc == 3 || argc == 4)) {
    invalidUse();
  }
  Problem problem;
  if (*argv[1] == 'D') {
    problem = Problem::Dirichlet;
  } else if (*argv[1] == 'N') {
    problem = Problem::Neumann;
  } else {
    invalidUse();
  }
  Domain d = problemDomain(problem);

  const double tol = std::stod(argv[2]);

  const int defaultNum = 8;
  const int num = (argc == 4) ? std::stoi(argv[3]) : defaultNum;

  const int nx = num;
  const int ny = num;

  const double dx = 1.0 / (nx);
  const double dy = 1.0 / (ny);

  const std::string name =
      problem == Problem::Dirichlet ? "Dirichlet" : "Neumann";
  precice::SolverInterface interface{name, "../precice-config.xml", 0, 1};
  const auto meshID = interface.getMeshID(name + "-Mesh");
  const auto tempID = interface.getDataID("Temperature", meshID);
  const auto fluxID = interface.getDataID("Heat-Flux", meshID);

  // defines the "inner mesh"
  auto vids = defineInterfaceMesh(interface, meshID, ny);
  std::vector<double> interfaceTemp(ny - 2);
  std::vector<double> interfaceFlux(ny - 2);
  assert(vids.size() == interfaceFlux.size() &&
         interfaceFlux.size() == interfaceTemp.size());

  auto u_grid = [&](int x, int y, double t) {
    return u(d.left + x * dx, d.bottom + y * dy, t);
  };

  Grid T{nx, ny, [u_grid](int x, int y) { return u_grid(x, y, 0); }};
  Grid T_last = T;
  Checkpoint checkpoint{T, 0};

  double dt = interface.initialize();

  if (problem == Problem::Neumann &&
      interface.isActionRequired(
          precice::constants::actionWriteInitialData())) {
    std::cout << "Initializing Temperature\n";
    T.writeInnerColTo(0, interfaceTemp.begin());
    interface.writeBlockScalarData(tempID, ny - 2, vids.data(),
                                   interfaceTemp.data());
    interface.markActionFulfilled(precice::constants::actionWriteInitialData());
  }

  interface.initializeData();

#ifndef NDEBUG
  std::cout << "Initial state\n" << T << "\n\n";
#endif

  int filecnt = 0;
  double t = 0;
  while (interface.isCouplingOngoing()) {

    // write checkpoint
    if (interface.isActionRequired(
            precice::constants::actionWriteIterationCheckpoint())) {
      std::cout << "Writing checkpoint\n";
      checkpoint.T = T_last;
      checkpoint.t = t;
      interface.markActionFulfilled(
          precice::constants::actionWriteIterationCheckpoint());
    }

    // update boundaryies
    auto u_grid_t = [t, u_grid](int x, int y) { return u_grid(x, y, t); };
    if (filecnt > 0) {
      T_last.evaluate(u_grid_t, Boundary::Top);
      T_last.evaluate(u_grid_t, Boundary::Bottom);
      if (problem == Problem::Dirichlet) {
        T_last.evaluate(u_grid_t, Boundary::Left);
      }
      if (problem == Problem::Neumann) {
        T_last.evaluate(u_grid_t, Boundary::Right);
      }
    }

    // Read data
    if (problem == Problem::Dirichlet) {
      interface.readBlockScalarData(tempID, ny - 2, vids.data(),
                                    interfaceTemp.data());
      T_last.readInnerColFrom(nx - 1, interfaceTemp.begin());
#ifndef NDEBUG
      std::cout << "Updated grid\n" << T << "\n\n";
#endif
    } else {
      interface.readBlockScalarData(fluxID, ny - 2, vids.data(),
                                    interfaceFlux.data());
      std::cout << "Received flux\n";
#ifndef NDEBUG
      for (double d : interfaceFlux)
        std::cout << d << '\n';
      std::cout << "\n";
#endif
    }

    std::cout << "Computing timestep: t " << t << " + " << dt << '\n';
    if (problem == Problem::Dirichlet) {
      timestep(T, T_last, dx, dy, dt, tol);
    } else {
      timestep(T, T_last, dx, dy, dt, tol, ApplyFlux{interfaceFlux.data(), dx});
    }
#ifndef NDEBUG
    std::cout << T << "\n\n";
#endif

    // Write data
    if (problem == Problem::Dirichlet) {
      std::cout << "Computing flux\n";
      computeFlux(T, dx, interfaceFlux.begin());
#ifndef NDEBUG
      for (double d : interfaceFlux)
        std::cout << d << "\n";
      std::cout << "\n";
#endif
      interface.writeBlockScalarData(fluxID, ny - 2, vids.data(),
                                     interfaceFlux.data());
    } else {
      T.writeInnerColTo(0, interfaceTemp.begin());
      interface.writeBlockScalarData(tempID, ny - 2, vids.data(),
                                     interfaceTemp.data());
    }

    double next_dt = interface.advance(dt);

    // write checkpoint
    if (interface.isActionRequired(
            precice::constants::actionReadIterationCheckpoint())) {
      interface.markActionFulfilled(
          precice::constants::actionReadIterationCheckpoint());
      std::cout << "Restoring checkpoint\n";
      T_last = checkpoint.T;
      t = checkpoint.t;
    } else {
      // timestep successfull
      T.writeVTK(name, d, filecnt);
      ++filecnt;
      T_last = T;
      t += dt;
      dt = next_dt;
    }
  }

  return 0;
}
