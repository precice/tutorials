#include <iostream>
#include <math.h>
#include <precice/SolverInterface.hpp>

using Vector = std::vector<double>;

struct DataContainer {
  void store_data(const Vector &vertices,
                  const double &theta,
                  const double &theta_dot)
  {
    old_vertices  = vertices;
    old_theta     = theta;
    old_theta_dot = theta_dot;
  }

  void reload_data(Vector &vertices,
                   double &theta,
                   double &theta_dot) const
  {
    vertices  = old_vertices;
    theta     = old_theta;
    theta_dot = old_theta_dot;
  }

  Vector old_vertices;
  double old_theta;
  double old_theta_dot;
};

class Solver {
public:
  Solver(const double moment_of_inertia)
      : moment_of_inertia(moment_of_inertia)
  {
  }

  void
  solve(const Vector &forces,
        const Vector &initial_vertices,
        Vector &      vertices,
        double &      theta,
        double &      theta_dot,
        const double  delta_t)
  {
    // Compute total moment m = x^{n} x f^{n+1}
    double moment = 0;
    for (uint i = 0; i < forces.size() / 2; ++i)
      moment += vertices[2 * i] * forces[2 * i + 1] - vertices[2 * i + 1] * forces[2 * i];

    // Store rigid body angle theta^{n}
    const double theta_old = theta;
    // Update angle to theta^{n+1} according to forward Euler method (simplified moment
    // computation, which does not depend on the updated configuration)
    theta = moment * std::pow(delta_t, 2) / moment_of_inertia + delta_t * theta_dot + theta;

    // Update angular velocity
    theta_dot = (theta - theta_old) / delta_t;

    // Update vertices according to rigid body rotation using an out-of-plane (z-axis) rotation matrix
    for (uint i = 0; i < vertices.size() / 2; ++i) {
      const double x_coord = initial_vertices[2 * i];
      vertices[2 * i]      = x_coord * std::cos(theta) + initial_vertices[2 * i + 1] * std::sin(theta);
      vertices[2 * i + 1]  = -x_coord * std::sin(theta) + initial_vertices[2 * i + 1] * std::cos(theta);
    }
    std::cout << "Theta: " << theta << " Theta dot: " << theta_dot << " Moment: " << moment << std::endl;
  }

private:
  const double moment_of_inertia;
};

int main()
{
  std::cout << "Rigid body: starting... \n";

  // Configuration settings
  const std::string config_file_name("precice-config.xml");
  const std::string solver_name("Solid");
  const std::string mesh_name("Solid-mesh");
  const std::string data_write_name("Displacement");
  const std::string data_read_name("Force");

  // Mesh configuration
  constexpr int vertical_refinement   = 3;
  constexpr int horizontal_refinement = 6;
  // Rotation centre is at (0,0)
  constexpr double length = 0.25;
  constexpr double height = 0.02;

  constexpr double density = 100;

  //*******************************************************************************************//

  // Derived quantities
  constexpr int n_vertical_nodes   = vertical_refinement * 2 + 1;
  constexpr int n_horizontal_nodes = horizontal_refinement * 2 + 1;
  // Substract shared nodes at each rigid body corner
  constexpr int    n_nodes        = (n_vertical_nodes + n_horizontal_nodes - 2) * 2;
  constexpr double mass           = length * height * density;
  constexpr double inertia_moment = (1. / 12) * mass * (4 * std::pow(length, 2) + std::pow(height, 2));

  // Create Solverinterface
  precice::SolverInterface precice(solver_name,
                                   config_file_name,
                                   /*comm_rank*/ 0,
                                   /*comm_size*/ 1);

  const int mesh_id  = precice.getMeshID(mesh_name);
  const int dim      = precice.getDimensions();
  const int write_id = precice.getDataID(data_write_name, mesh_id);
  const int read_id  = precice.getDataID(data_read_name, mesh_id);

  // Set up data structures
  Vector           forces(dim * n_nodes);
  Vector           vertices(dim * n_nodes);
  Vector           displacement(dim * n_nodes);
  std::vector<int> vertex_ids(n_nodes);
  double           theta_dot  = 0.0;
  double           theta      = 0.0;

  {
    // Define a boundary mesh
    std::cout << "Rigid body: defining mesh: ";
    std::cout << n_nodes << " nodes " << std::endl;
    const double delta_y = height / (n_vertical_nodes - 1);
    const double delta_x = length / (n_horizontal_nodes - 1);

    // x planes
    for (int i = 0; i < n_vertical_nodes; ++i) {
      // negative x
      vertices[dim * i]     = 0.0;                         // fixed x
      vertices[dim * i + 1] = -height * 0.5 + delta_y * i; // increasing y
      // positive x
      vertices[(2 * n_vertical_nodes) + dim * i]     = length;                // fixed x
      vertices[(2 * n_vertical_nodes) + dim * i + 1] = vertices[dim * i + 1]; // increasing y
    }

    // y planes
    const unsigned int of = 2 * dim * n_vertical_nodes; // static offset
    // Lower and upper bounds are already included due to positive/negative x-planes
    const unsigned int n_remaining_nodes = n_horizontal_nodes - 2;
    for (int i = 0; i < n_remaining_nodes; ++i) {
      // negative y
      vertices[of + dim * i]     = delta_x * (i + 1); // increasing x
      vertices[of + dim * i + 1] = -height * 0.5;     // fixed y
      // positive y
      vertices[of + (2 * n_remaining_nodes) + dim * i]     = vertices[of + dim * i];      // increasing x
      vertices[of + (2 * n_remaining_nodes) + dim * i + 1] = -vertices[of + dim * i + 1]; // fixed y
    }
  }
  // Store the initial configuration
  const Vector initial_vertices = vertices;

  // Pass the vertices to preCICE
  precice.setMeshVertices(mesh_id,
                          n_nodes,
                          vertices.data(),
                          vertex_ids.data());

  // initialize the Solverinterface
  double dt = precice.initialize();

  // Compute the absolute displacement between the current vertices and the
  // initial configuration. Here, this is mostly done for consistency reasons.
  for (uint i = 0; i < displacement.size(); ++i)
    displacement[i] = vertices[i] - initial_vertices[i];

  if (precice.isActionRequired(precice::constants::actionWriteInitialData())) {
    precice.writeBlockVectorData(write_id,
                                 n_nodes,
                                 vertex_ids.data(),
                                 displacement.data());

    precice.markActionFulfilled(
        precice::constants::actionWriteInitialData());

    precice.initializeData();
  }

  std::cout << "Rigid body: reading initial data \n";
  if (precice.isReadDataAvailable())
    precice.readBlockVectorData(read_id,
                                n_nodes,
                                vertex_ids.data(),
                                forces.data());

  // Set up a struct in order to store time dependent values
  DataContainer data_container;
  // Set up an object which handles the time integration
  Solver solver(inertia_moment);

  // Start time loop
  double time = 0;
  while (precice.isCouplingOngoing()) {

    std::cout << "Rigid body: t = " << time << "s \n";

    // Store time dependent values
    if (precice.isActionRequired(
            precice::constants::actionWriteIterationCheckpoint())) {
      data_container.store_data(vertices, theta, theta_dot);

      precice.markActionFulfilled(
          precice::constants::actionWriteIterationCheckpoint());
    }

    // Solve system
    solver.solve(forces, initial_vertices, vertices, theta, theta_dot, dt);

    // Advance coupled system
    {
      // Compute absolute displacement with respect to the initial configuration
      for (uint i = 0; i < displacement.size(); ++i)
        displacement[i] = vertices[i] - initial_vertices[i];

      std::cout << "Rigid body: writing coupling data \n";
      if (precice.isWriteDataRequired(dt))
        precice.writeBlockVectorData(write_id,
                                     n_nodes,
                                     vertex_ids.data(),
                                     displacement.data());

      std::cout << "Rigid body: advancing in time\n";
      dt = precice.advance(dt);

      std::cout << "Rigid body: reading coupling data \n";
      if (precice.isReadDataAvailable())
        precice.readBlockVectorData(read_id,
                                    n_nodes,
                                    vertex_ids.data(),
                                    forces.data());
    }

    // Reload time dependent values
    if (precice.isActionRequired(
            precice::constants::actionReadIterationCheckpoint())) {
      data_container.reload_data(vertices, theta, theta_dot);

      precice.markActionFulfilled(
          precice::constants::actionReadIterationCheckpoint());
    }

    // Increment time in case the time window has been completed
    if (precice.isTimeWindowComplete())
      time += dt;
  }

  std::cout << "Rigid body: closing...\n";

  return 0;
}
