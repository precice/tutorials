/*
 *                flow-over-heated-plate/solid-mfem
 *
 * This file implements the solid heat equation solver for the preCICE
 * flow-over-heated-plate tutorial using MFEM in parallel. Only
 * AssemblyLevel::LEGACY is supported. It is recommended to review Example 16p
 * in the MFEM repository prior to this.
 *
 * The semi-discrete form can be written as:
 *
 *                       M (du/dt) = -Ku + N,
 *
 * where M is a mass matrix, K is a stiffness matrix, N is the Neumann
 * term, and u is the temperature solution vector.
 *
 * For explicit time integration in TimeDependentOperator::Mult, the
 * equation to solve for (du/dt) for time t_n is simply
 *
 *                (du/dt) = M^{-1}(-Ku_n + N).
 *
 * For implicit time integration in TimeDependentOperator::ImplicitSolve, the
 * equation to solve for k=du/dt can be written as
 *
 *                (M + dt K)(du/dt) = -Ku_n + N.
 *
 * Writing A = M + dt K, the equation is then
 *
 *                (du/dt) = A^{-1}(-Ku_n + N).
 *
 * In this case, M and K are constant, A depends on dt, and N depends on the
 * coupling data. For essential BCs, it is enforced that (du/dt) = 0.
 *
 * For the coupled heat flux, we set the interface mesh at the solution nodes,
 * and then GridFunctionCoefficient such that the received heat flux is then
 * interpolated and quadrature is performed using the interpolated heat flux
 * at quadrature points.
 *
 * In regards to parallelization, we use unique DOFs across all ranks.
 */

#include <mfem/mfem.hpp>
#include <precice/precice.hpp>

using namespace mfem;
using namespace precice;

class LinearHeatOperator : public TimeDependentOperator {
protected:
  const Array<int> &ess_tdof_list_;

  CGSolver      M_solver_;
  CGSolver      A_solver_;
  HypreSmoother M_prec_;
  HypreSmoother A_prec_;

  ParBilinearForm                 M_;
  std::unique_ptr<HypreParMatrix> M_mat_;

  ParBilinearForm                 K_;
  std::unique_ptr<HypreParMatrix> K_mat_;

  ParLinearForm b_;
  Vector        b_tvec_;

  std::unique_ptr<HypreParMatrix> A_mat_;

  mutable Vector rhs_;

  void EvalRHS(const Vector &u, Vector &rhs) const;

public:
  LinearHeatOperator(ParFiniteElementSpace &f,
                     const Array<int>      &ess_tdof_list,
                     Coefficient           &alpha,
                     Coefficient           &q_bdr);

  void Mult(const Vector &u, Vector &du_dt) const override;

  void ImplicitSolve(const real_t dt, const Vector &u, Vector &k) override;

  void UpdateNeumann();
};

int main(int argc, char *argv[])
{
  mfem::Mpi::Init(argc, argv);
  int size = mfem::Mpi::WorldSize();
  int rank = mfem::Mpi::WorldRank();
  mfem::Hypre::Init();

  int         order           = 1;
  int         ode_solver_type = 23; // SDIRK33Solver
  real_t      dt              = 0.01;
  int         pvdc_freq       = 20;
  std::string precice_config  = "../precice-config.xml";
  bool        visualization   = true;
  int         visport         = 19916;
  char        vishost[]       = "localhost";

  OptionsParser args(argc, argv);
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                 ODESolver::Types.c_str());
  args.AddOption(&order, "-dt", "--time-step",
                 "Time step.");
  args.AddOption(&pvdc_freq, "-pv", "--paraview-freq",
                 "ParaView data collection write frequency. 0 to disable.");
  args.AddOption(&precice_config, "-c", "--config-file",
                 "preCICE configuration file.");
  args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                 "--no-visualization",
                 "Enable or disable GLVis visualization.");
  args.AddOption(&visport, "-p", "--send-port", "Socket for GLVis.");
  args.Parse();
  if (!args.Good()) {
    if (rank == 0) {
      args.PrintUsage(std::cout);
    }
    return 1;
  }
  if (rank == 0) {
    args.PrintOptions(std::cout);
  }

  // Create mesh.
  // Bdr attributes are:
  // 1 = bottom, 2 = right, 3 = top, 4 = left.
  const int           nx        = 20;
  const int           ny        = 5;
  const Element::Type e_type    = Element::QUADRILATERAL;
  const bool          gen_edges = true;
  const real_t        Lx        = 1.0;
  const real_t        Ly        = 0.25;
  Mesh                mesh      = Mesh::MakeCartesian2D(nx, ny, e_type, gen_edges, Lx, Ly);
  mesh.Transform(
      [&](const Vector &x_old, Vector &x_new) {
        // Ensure top of plate at y = 0.
        x_new[1] = x_old[1] - Ly;
      });

  // Partition grid.
  ParMesh pmesh(MPI_COMM_WORLD, mesh);

  // Use 2D H1-conforming finite elements.
  const int       dim = 2;
  H1_FECollection fecoll(order, dim);

  Array<int>          interface_dofs;   // rank-unique ldofs (ldofs of the tdofs)
  std::vector<double> interface_coords; // preCICE requires double
  // Get the coordinates + indices of interface nodes on this rank.
  {
    // Create an FE space for the obtaining the nodal coordinates (2D).
    // Use Ordering::byVDIM to ensure XYZ XYZ ...
    ParFiniteElementSpace fes_nodes(&pmesh, &fecoll, 2, Ordering::byVDIM);

    // Get the interface coordinates.
    // (If using p > 1, SetNodalFESpace() is necessary).
    pmesh.SetNodalFESpace(&fes_nodes);
    const GridFunction &coords = *pmesh.GetNodes();

    // Get the interface DOFs/indices
    Array<int> interface_bdr_elems, interface_tdofs;
    fes_nodes.GetBoundaryElementsByAttribute(3, interface_bdr_elems);
    fes_nodes.GetBoundaryLoopEdgeDofs(interface_bdr_elems, interface_tdofs,
                                      interface_dofs);

    // Convert ldofs to vdofs (include indices for y-coords).
    // Indices are appended to end.
    Array<int> interface_vdofs = interface_dofs;
    fes_nodes.DofsToVDofs(interface_vdofs);

    // Sort indices so that x,y indices are adjacent.
    interface_vdofs.Sort();

    // Get the coordinates, correctly organized.
    Vector vec_coords;
    coords.GetSubVector(interface_vdofs, vec_coords);
    interface_coords.resize(vec_coords.Size());
    for (int i = 0; i < vec_coords.Size(); i++) {
      interface_coords[i] = vec_coords[i];
    }
  }

  // Initialize FE space for temperature solution.
  // Use same ordering to ensure indices `interface_dofs` match.
  ParFiniteElementSpace fespace(&pmesh, &fecoll, 1, Ordering::byVDIM);

  // Initialize temperature with IC
  ParGridFunction u_gf(&fespace);
  u_gf = 300;

  // Lower bdr is essential, fixed at T=300 K.
  // Get tdofs.
  Array<int> bdr_marker_arr(pmesh.bdr_attributes.Size());
  bdr_marker_arr                         = 0;
  bdr_marker_arr[bdr_marker_arr.Find(1)] = 1;
  Array<int> ess_tdof_list;
  fespace.GetEssentialTrueDofs(bdr_marker_arr, ess_tdof_list);

  // Initialize coefficient for alpha.
  ConstantCoefficient alpha(1.0);

  // Initialize interface boundary heat flux coefficient.
  ParGridFunction         q_interface(&fespace);
  GridFunctionCoefficient q_interface_coeff(&q_interface);

  // Initialize Neumann term coefficient.
  // Only add coefficient for interface.
  PWCoefficient q_bdr;
  q_bdr.UpdateCoefficient(3, q_interface_coeff);

  // Initialize the operator.
  LinearHeatOperator oper(fespace, ess_tdof_list, alpha, q_bdr);

  // Initialize preCICE
  Participant      participant("Solid", precice_config, rank, size);
  std::vector<int> mesh_vertices(interface_coords.size() / dim);

  participant.setMeshVertices("Solid", interface_coords, mesh_vertices);

  return 0;
}

LinearHeatOperator::LinearHeatOperator(
    ParFiniteElementSpace &f,
    const Array<int>      &ess_tdof_list,
    Coefficient           &alpha,
    Coefficient           &q_bdr)
    : ess_tdof_list_(ess_tdof_list),
      M_solver_(f.GetComm()),
      A_solver_(f.GetComm()),
      M_(&f),
      K_(&f),
      b_(&f),
      b_tvec_(f.GetTrueVSize()),
      rhs_(f.GetTrueVSize())
{

#if defined(MFEM_USE_DOUBLE)
  const real_t rel_tol = 1e-12;
#elif defined(MFEM_USE_SINGLE)
  const real_t rel_tol = 1e-6;
#else
#error "Only single and double are supported!"
  const real_t rel_tol = 1e-12;
#endif
  const int max_iter = 100;

  M_.AddDomainIntegrator(new MassIntegrator);
  M_.Assemble(0);
  M_.Finalize(0);
  M_mat_ = std::unique_ptr<HypreParMatrix>(M_.ParallelAssemble());
  M_mat_->EliminateBC(ess_tdof_list, Operator::DiagonalPolicy::DIAG_ONE);

  K_.AddDomainIntegrator(new DiffusionIntegrator(alpha));
  K_.Assemble(0);
  K_.Finalize(0);
  K_mat_ = std::unique_ptr<HypreParMatrix>(K_.ParallelAssemble());

  b_.AddBoundaryIntegrator(new BoundaryLFIntegrator(q_bdr));
  UpdateNeumann();

  M_solver_.iterative_mode = false;
  M_solver_.SetRelTol(rel_tol);
  M_solver_.SetAbsTol(0.0);
  M_solver_.SetMaxIter(max_iter);
  M_solver_.SetPrintLevel(IterativeSolver::PrintLevel().Warnings().Errors());
  M_solver_.SetPreconditioner(M_prec_);
  M_solver_.SetOperator(*M_mat_);

  A_solver_.iterative_mode = false;
  A_solver_.SetRelTol(rel_tol);
  A_solver_.SetAbsTol(0.0);
  A_solver_.SetMaxIter(max_iter);
  A_solver_.SetPrintLevel(IterativeSolver::PrintLevel().Warnings().Errors());
  A_solver_.SetPreconditioner(A_prec_);
}

void LinearHeatOperator::EvalRHS(const Vector &u, Vector &rhs) const
{
  rhs = 0.0;
  K_mat_->Mult(u, rhs);
  rhs.Neg();
  rhs += b_tvec_;
}

void LinearHeatOperator::Mult(const Vector &u, Vector &du_dt) const
{
  EvalRHS(u, rhs_);

  // Application of essential BC for M already done.
  rhs_.SetSubVector(ess_tdof_list_, 0.0);

  M_solver_.Mult(rhs_, du_dt);
}

void LinearHeatOperator::ImplicitSolve(const real_t dt, const Vector &u,
                                       Vector &k)
{
  A_mat_.reset();
  A_mat_ = std::unique_ptr<HypreParMatrix>(Add(1.0, *M_mat_, dt, *K_mat_));

  EvalRHS(u, rhs_);

  A_mat_->EliminateBC(ess_tdof_list_, Operator::DiagonalPolicy::DIAG_ONE);
  rhs_.SetSubVector(ess_tdof_list_, 0.0);

  A_solver_.SetOperator(*A_mat_);
  A_solver_.Mult(rhs_, k);
}

void LinearHeatOperator::UpdateNeumann()
{
  b_.Assemble();
  b_.ParallelAssemble(b_tvec_);
}