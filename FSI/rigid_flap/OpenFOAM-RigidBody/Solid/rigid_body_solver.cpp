#include <iostream>
#include <math.h>
#include <precice/SolverInterface.hpp>

using Vector = std::vector<double>;

struct DataContainer {
  void store_data(const Vector &vertices, const double &old_moment,
                  const double &theta, const double &theta_dot)
  {
    old_vertices   = vertices;
    old_old_moment = old_moment;
    old_theta      = theta;
    old_theta_dot  = theta_dot;
  }

  void reload_data(Vector &vertices, double &old_moment,
                   double &theta, double &theta_dot) const
  {
    vertices   = old_vertices;
    old_moment = old_old_moment;
    theta      = old_theta;
    theta_dot  = old_theta_dot;
  }

  Vector old_vertices;
  double old_old_moment;
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
        double &,
        double &     theta,
        double &     theta_dot,
        const double delta_t)
  {
    // Compute total moment
    double moment = 0;
    for (uint i = 0; i < forces.size() / 2; ++i)
      moment += vertices[2 * i] * forces[2 * i + 1] - vertices[2 * i + 1] * forces[2 * i];
    // Store old angle
    const double theta_old = theta;
    // Update angle
    theta = moment * std::pow(delta_t, 2) / moment_of_inertia + delta_t * theta_dot + theta;
//    theta = -theta;
    // Update velocity
    theta_dot = (theta - theta_old) / delta_t;
    // Update vertices according to rigid body rotation
    for (uint i = 0; i < vertices.size() / 2; ++i) {
      //compute length
      const double x_coord = initial_vertices[2 * i];
      vertices[2 * i]      = x_coord * std::cos(theta) + initial_vertices[2 * i + 1] * std::sin(theta);
      vertices[2 * i + 1]  = -x_coord * std::sin(theta) + initial_vertices[2 * i + 1] * std::cos(theta);
    }
    // Compute recent moment on updated coordinates

    // Update recent angular velocity
    //    theta_dot = theta_dot + 0.5 * (moment + old_moment) * delta_t;
    //    theta_dot = -theta_dot;
    std::cout << "Theta: " << theta << " Theta dot: " << theta_dot << " Moment: " << moment << std::endl;
    // Update moment
    //    old_moment = moment;
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
  constexpr int    vertical_refinement   = 3;
  constexpr int    horizontal_refinement = 6;
  constexpr double length                = 0.25;
  constexpr double height                = 0.02;
  // Rotation centre is at (0,0)

  constexpr double density = 100;
  //***********************************************************************************//

  // Derived quantities
  // Substract common nodes at each corner
  constexpr int    n_vertical_nodes   = vertical_refinement * 2 + 1;
  constexpr int    n_horizontal_nodes = horizontal_refinement * 2 + 1;
  constexpr int    n_nodes            = (n_vertical_nodes + n_horizontal_nodes - 2) * 2;
  constexpr double mass               = length * height * density;
  constexpr double inertia_moment     = (1. / 12) * mass * (4 * std::pow(length, 2) + std::pow(height, 2));

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
  double           old_moment = 0.0;
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
    // Lower and upper bounds are already included
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
  const Vector initial_vertices = vertices;
  // Pass the vertices to preCICE
  precice.setMeshVertices(mesh_id,
                          n_nodes,
                          vertices.data(),
                          vertex_ids.data());

  // initialize the Solverinterface
  double dt = precice.initialize();

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

  DataContainer data_container;
  Solver        solver(inertia_moment);

  // Start time loop
  double time = 0;
  while (precice.isCouplingOngoing()) {

    // Store
    if (precice.isActionRequired(
            precice::constants::actionWriteIterationCheckpoint())) {
      data_container.store_data(vertices, old_moment, theta, theta_dot);

      precice.markActionFulfilled(
          precice::constants::actionWriteIterationCheckpoint());
    }
    // Solve system
    solver.solve(forces, initial_vertices, vertices, old_moment, theta, theta_dot, dt);
    // Advance
    {
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

    // Reload
    if (precice.isActionRequired(
            precice::constants::actionReadIterationCheckpoint())) {
      data_container.reload_data(vertices, old_moment, theta, theta_dot);

      precice.markActionFulfilled(
          precice::constants::actionReadIterationCheckpoint());
    }
    // Increment
    if (precice.isTimeWindowComplete())
      time += dt;
  }

  std::cout << "Rigid body: closing...\n";

  return 0;
}
