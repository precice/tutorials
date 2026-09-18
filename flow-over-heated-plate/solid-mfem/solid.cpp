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
 * term, and u is the temperature solution vector. In this case, N=0 (natural).
 *
 * For explicit time integration in TimeDependentOperator::Mult, the
 * equation to solve for (du/dt) for time t_n is simply
 *
 *                (du/dt) = M^{-1}(-Ku_n).
 *
 * For implicit time integration in TimeDependentOperator::ImplicitSolve, the
 * equation to solve for k=du/dt can be written as
 *
 *                (M + dt K)(du/dt) = -Ku_n.
 *
 * Writing A = M + dt K, the equation is then
 *
 *                (du/dt) = A^{-1}(-Ku_n).
 *
 * In this case, M and K are constant, and A depends on dt. For essential BCs,
 * it is enforced that (du/dt) = 0. Note that this is an approximation for the
 * interface, as it indeed does vary in time.
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

#ifndef MFEM_USE_DOUBLE
#error "Tutorial requires MFEM built with double-precision for consistency \
        with preCICE."
#endif

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

  std::unique_ptr<HypreParMatrix> A_mat_;

  mutable Vector rhs_;

  void EvalRHS(const Vector &u, Vector &rhs) const;

public:
  LinearHeatOperator(ParFiniteElementSpace &f,
                     const Array<int>      &ess_tdof_list,
                     Coefficient           &alpha);

  void Mult(const Vector &u, Vector &du_dt) const override;

  void ImplicitSolve(const real_t dt, const Vector &u, Vector &k) override;
};

// Class for extracting wall normal fluxes on a given boundary.
class WallNormalFluxExtractor {
private:
  const FiniteElementSpace &fes_;

  // Boundary element indices.
  Array<int> bdr_elems_;

  // Boundary dofs -> edofs for associated boundary element in bdr_elems_.
  Array<int> bdr_edofs_;

public:
  WallNormalFluxExtractor(const FiniteElementSpace &fes,
                          int                       bdr_attr,
                          const Array<int>         &bdr_dofs);

  // fluxes of size bdr_dofs
  void GetWallNormalFluxes(const GridFunction &gf, real_t *fluxes) const;
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

  Array<int>          interface_dofs; // rank-unique ldofs (ldofs of the tdofs)
  std::vector<double> interface_coords;
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
    interface_coords.resize(interface_vdofs.Size());
    coords.GetSubVector(interface_vdofs, interface_coords.data());
  }

  // Initialize FE space for temperature solution.
  // Use same ordering to ensure indices `interface_dofs` match.
  ParFiniteElementSpace fespace(&pmesh, &fecoll, 1, Ordering::byVDIM);

  // Initialize temperature with IC
  ParGridFunction u_gf(&fespace);
  u_gf = 300;

  // Lower + upper boundaries are essential --> Get tdofs.
  Array<int> bdr_marker_arr(pmesh.bdr_attributes.Size());
  bdr_marker_arr                         = 0;
  bdr_marker_arr[bdr_marker_arr.Find(1)] = 1; // lower
  bdr_marker_arr[bdr_marker_arr.Find(3)] = 1; // upper
  Array<int> ess_tdof_list;
  fespace.GetEssentialTrueDofs(bdr_marker_arr, ess_tdof_list);

  // Initialize coefficient for alpha.
  ConstantCoefficient alpha(1.0);

  // Initialize the operator.
  LinearHeatOperator oper(fespace, ess_tdof_list, alpha);

  // Initialize ODE Solver
  std::unique_ptr<ODESolver> ode_solver = ODESolver::Select(ode_solver_type);
  ode_solver->Init(oper);
  real_t t  = 0.0;
  int    ti = 0;

  // Initialize ParaViewDataCollection for output.
  std::unique_ptr<ParaViewDataCollection> pvdc;
  if (pvdc_freq > 0) {
    pvdc = std::make_unique<ParaViewDataCollection>("Solid", &pmesh);
    pvdc->SetLevelsOfDetail(order);
    pvdc->SetDataFormat(VTKFormat::BINARY);
    pvdc->SetHighOrderOutput(true);
    pvdc->SetCycle(0);
    pvdc->SetTime(0.0);
    pvdc->RegisterField("Temperature", &u_gf);
    pvdc->Save();
  }

  // Initialize preCICE
  const std::string mesh_name = "Solid-Mesh";
  Participant       participant("Solid", precice_config, rank, size);
  std::vector<int>  mesh_vertices(interface_coords.size() / dim);
  participant.setMeshVertices(mesh_name, interface_coords, mesh_vertices);
  participant.initialize();

  // Initialize preCICE-related variables.
  real_t                        precice_dt;
  real_t                        t_save;
  ParGridFunction               u_gf_save(&fespace);
  std::vector<double>           u_receive(interface_dofs.Size());
  std::vector<double>           qwall_write(interface_dofs.Size());
  const WallNormalFluxExtractor qwall_getter(fespace, 3, interface_dofs);

  // Main solver loop.
  while (participant.isCouplingOngoing()) {
    if (participant.requiresWritingCheckpoint()) {
      t_save    = t;
      u_gf_save = u_gf;
    }

    precice_dt = participant.getMaxTimeStepSize();
    dt         = std::min(dt, precice_dt);

    // Get coupling data + set
    participant.readData(mesh_name, "Temperature", mesh_vertices, dt, u_receive);
    u_gf.SetSubVector(interface_dofs, u_receive.data());

    // Step in time
    ode_solver->Step(u_gf.GetTrueVector(), t, dt);
    u_gf.SetFromTrueVector();

    // Get and write the heat fluxes.
    qwall_getter.GetWallNormalFluxes(u_gf, qwall_write.data());
    participant.writeData(mesh_name, "Heat-Flux", mesh_vertices, qwall_write);

    if (participant.requiresReadingCheckpoint()) {
      t    = t_save;
      u_gf = u_gf_save;
    } else {
      if (rank == 0) {
        std::cout << "step " << ti << ", t = " << t << std::endl;
      }
      ti++;
      if (pvdc && ti % pvdc_freq == 0) {
        pvdc->SetTime(t);
        pvdc->SetCycle(ti);
        pvdc->Save();
      }
    }
  }

  participant.finalize();

  return 0;
}

LinearHeatOperator::LinearHeatOperator(
    ParFiniteElementSpace &f,
    const Array<int>      &ess_tdof_list,
    Coefficient           &alpha)
    : ess_tdof_list_(ess_tdof_list),
      M_solver_(f.GetComm()),
      A_solver_(f.GetComm()),
      M_(&f),
      K_(&f),
      rhs_(f.GetTrueVSize())
{

  // Double-precision relative tolerance:
  const real_t rel_tol = 1e-12;

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

WallNormalFluxExtractor::WallNormalFluxExtractor(
    const FiniteElementSpace &fes,
    int                       bdr_attr,
    const Array<int>         &bdr_dofs)
    : fes_(fes),
      bdr_elems_(bdr_dofs.Size()),
      bdr_edofs_(bdr_dofs.Size())
{
  // Loop over all boundary elements
  // Save boundary element and associated edof for every bdr_dof
  for (int i = 0; i < fes.GetNBE(); i++) {
    // Check if this bdr element is associated with bdr_attr
    if (fes.GetBdrAttribute(i) != bdr_attr) {
      continue;
    }

    Array<int> edofs;
    fes.GetBdrElementDofs(i, edofs);

    for (int j = 0; j < edofs.Size(); j++) {
      int idx = bdr_dofs.Find(edofs[j]);

      // Skip if not in bdr_dofs
      if (idx < 0) {
        continue;
      }

      bdr_elems_[idx] = i;
      bdr_edofs_[idx] = j;
    }
  }
}

void WallNormalFluxExtractor::GetWallNormalFluxes(const GridFunction &gf,
                                                  real_t             *fluxes) const
{
  for (int i = 0; i < bdr_elems_.Size(); i++) {

    const FiniteElement   &fe = *fes_.GetBE(bdr_elems_[i]);
    ElementTransformation &tr =
        *fes_.GetBdrElementTransformation(bdr_elems_[i]);

    // Set point to compute normal at
    const IntegrationPoint &ip = fe.GetNodes().IntPoint(bdr_edofs_[i]);
    tr.SetIntPoint(&ip);

    // Compute the normal
    Vector normal(tr.Jacobian().Height());
    CalcOrtho(tr.Jacobian(), normal);

    // Get the gradient of the grid function at the point
    Vector grad_u(normal.Size());
    gf.GetGradient(tr, grad_u);

    // Compute grad_u dot normal
    fluxes[i] = grad_u * normal;
  }
}