#include <cmath>
#include <iostream>
#include <precice/SolverInterface.hpp>
#include <string>

using Vector = std::vector<double>;

struct DataContainer {
  void save_old_state(const Vector &vertices,
                      const double &theta,
                      const double &theta_dot)
  {
    old_vertices  = vertices;
    old_theta     = theta;
    old_theta_dot = theta_dot;
  }

  void reload_old_state(Vector &vertices,
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
        const double  spring_constant,
        const double  delta_t) const
  {
    // Compute total moment M = x^{n} x f^{n+1}
    double moment = 0;
    for (unsigned int i = 0; i < forces.size() / 2; ++i)
      moment += vertices[2 * i] * forces[2 * i + 1] - vertices[2 * i + 1] * forces[2 * i];

    // Store rigid body angle at the previous time level theta^{n}
    const double theta_old = theta;

    // Update angle to theta^{n+1} according to forward Euler method (simplified moment
    // computation, which does not depend on the updated configuration). The contribution
    // of the spring with strength k is discretized implicitly.
    // theta^{n+1} = 1/ (1 - (k/I)*dt^2) * dt^2 * M / I + dt * \dot{theta}^{n} + theta^{n}
    theta = (1. / (1 - (spring_constant / moment_of_inertia) * std::pow(delta_t, 2))) *
            (std::pow(delta_t, 2) * moment / moment_of_inertia + delta_t * theta_dot + theta);

    // Update angular velocity
    // \dot{theta}^{n+1} = (theta^{n+1} - \dot{theta}^{n}) / dt
    theta_dot = (theta - theta_old) / delta_t;

    // Update vertices according to rigid body rotation using an out-of-plane (z-axis) rotation matrix
    // x^{n+1} = x^{0} * cos(theta^{n+1}) + y^{0} * sin(theta^{n+1})
    // y^{n+1} = -x^{0} * sin(theta^{n+1}) + y^{0} * cos(theta^{n+1})
    for (uint i = 0; i < vertices.size() / 2; ++i) {
      const double x_coord = initial_vertices[2 * i];
      vertices[2 * i]      = x_coord * std::cos(theta) + initial_vertices[2 * i + 1] * std::sin(theta);
      vertices[2 * i + 1]  = -x_coord * std::sin(theta) + initial_vertices[2 * i + 1] * std::cos(theta);
    }
    std::cout << "Theta: " << theta << " Theta dot: " << theta_dot << " Moment: " << moment << " Spring force: " << spring_constant * theta << std::endl;
  }

private:
  const double moment_of_inertia;
};

int main()
{
  std::cout << "Rigid body: starting... \n";

  // Configuration settings
  const std::string config_file_name("../precice-config.xml");
  const std::string solver_name("Solid");
  const std::string mesh_name("Solid-Mesh");
  const std::string data_write_name("Displacement");
  const std::string data_read_name("Force");

  // Mesh configuration
  constexpr int vertical_refinement   = 3;
  constexpr int horizontal_refinement = 6;
  // Rotation centre is at (0,0)
  constexpr double length = 0.2;
  constexpr double height = 0.02;

  constexpr double density         = 10000;
  constexpr double spring_constant = -25;

  // Time, where spring is stiffened
  constexpr double switch_time       = 1.5;
  constexpr double stiffening_factor = 8;
  //*******************************************************************************************//

  // Derived quantities
  constexpr int n_vertical_nodes   = vertical_refinement * 2 + 1;
  constexpr int n_horizontal_nodes = horizontal_refinement * 2 + 1;
  // Substract shared nodes at each rigid body corner
  constexpr int    n_nodes = (n_vertical_nodes + n_horizontal_nodes - 2) * 2;
  constexpr double mass    = length * height * density;
  // The moment of inertia is computed according to the rigid body configuration:
  // a thin rectangular plate of height h, length l and mass m with axis of rotation
  // at the end of the plate: I = (1/12)*m*(4*l^2+h^2)
  constexpr double inertia_moment = (1. / 12) * mass * (4 * std::pow(length, 2) + std::pow(height, 2));
  constexpr double delta_y        = height / (n_vertical_nodes - 1);
  constexpr double delta_x        = length / (n_horizontal_nodes - 1);

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
  double           theta_dot = 0.0;
  double           theta     = 0.0;

  {
    // Define a boundary mesh
    std::cout << "Rigid body: defining mesh: ";
    std::cout << n_nodes << " nodes " << std::endl;

    //                  upper y
    //          o---o---o---o---o---o---o
    //          |.......................|
    // lower x  o.......................o   upper x
    //          |.......................|
    //          o---o---o---o---o---o---o
    //                  lower y

    // x planes
    for (int i = 0; i < n_vertical_nodes; ++i) {
      // lower x plane
      vertices[dim * i]     = 0.0;                         // fixed x
      vertices[dim * i + 1] = -height * 0.5 + delta_y * i; // increasing y
      // upper x plane
      vertices[(2 * n_vertical_nodes) + dim * i]     = length;                // fixed x
      vertices[(2 * n_vertical_nodes) + dim * i + 1] = vertices[dim * i + 1]; // increasing y
    }

    // y planes
    const unsigned int of = 2 * dim * n_vertical_nodes; // static offset
    // Lower and upper bounds are already included due to positive/negative x-planes
    const unsigned int n_remaining_nodes = n_horizontal_nodes - 2;
    for (unsigned int i = 0; i < n_remaining_nodes; ++i) {
      // lower y plane
      vertices[of + dim * i]     = delta_x * (i + 1); // increasing x
      vertices[of + dim * i + 1] = -height * 0.5;     // fixed y
      // upper y plane
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

  // Set up a struct in order to store time dependent values
  DataContainer data_container;
  // Set up an object which handles the time integration
  const Solver solver(inertia_moment);

  // Start time loop
  double time = 0;
  while (precice.isCouplingOngoing()) {

    std::cout << "Rigid body: t = " << time << "s \n";

    std::cout << "Rigid body: reading initial data \n";
    if (precice.isReadDataAvailable())
      precice.readBlockVectorData(read_id,
                                  n_nodes,
                                  vertex_ids.data(),
                                  forces.data());

    // Store time dependent values
    if (precice.isActionRequired(
            precice::constants::actionWriteIterationCheckpoint())) {
      data_container.save_old_state(vertices, theta, theta_dot);

      precice.markActionFulfilled(
          precice::constants::actionWriteIterationCheckpoint());
    }

    const double current_spring = time > switch_time ? spring_constant * stiffening_factor : spring_constant;
    // Solve system
    solver.solve(forces, initial_vertices, vertices, theta, theta_dot, current_spring, dt);

    // Advance coupled system
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

    // Reload time dependent values
    if (precice.isActionRequired(
            precice::constants::actionReadIterationCheckpoint())) {
      data_container.reload_old_state(vertices, theta, theta_dot);

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
