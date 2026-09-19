/*
 *                flow-over-heated-plate/solid-mfem
 *
 * This file implements the solid heat equation solver for the preCICE
 * flow-over-heated-plate tutorial using MFEM in parallel, where the solid
 * recieves the temperature and sends the heat flux. In regards to
 * parallelization, we use unique DOFs across all ranks. Only
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
 */

#include <mfem/mfem.hpp>
#include <precice/precice.hpp>

#ifndef MFEM_USE_DOUBLE
#error "Tutorial requires MFEM built with double-precision for consistency \
        with preCICE."
#endif

#ifndef MFEM_USE_MPI
#error "Tutorial uses MFEM built with MPI."
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

// Adapted from miniapps/common/pfem-extras
void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    const ParGridFunction &gf, const char *title,
                    int x = 0, int y = 0, int w = 400, int h = 400,
                    const char *keys = NULL, bool vec = false,
                    const char *commands = NULL);

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
  const int           nx        = 100;
  const int           ny        = 25;
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

  // Initialize FE space for temperature solution.
  ParFiniteElementSpace fespace(&pmesh, &fecoll, 1, Ordering::byVDIM);

  // For setting temperature, need ldof (index) of each interface node.
  Array<int> interface_dofs;

  // For getting heat flux, need bdr element and element-local dof of each
  // interface node.
  Array<int> interface_elems;
  Array<int> interface_eldofs;

  // Loop over all boundary elements
  for (int i = 0; i < fespace.GetNBE(); i++) {
    // Skip if not interface boundary element
    if (fespace.GetBdrAttribute(i) != 3) {
      continue;
    }

    // Get ldofs from this element
    Array<int> elem_dofs;
    fespace.GetBdrElementDofs(i, elem_dofs);

    for (int j = 0; j < elem_dofs.Size(); j++) {
      int ldof = elem_dofs[j];

      // Avoid duplication at bdr element ends
      if (interface_dofs.Find(ldof) >= 0) {
        continue;
      }

      // Avoid duplication at MPI rank interfaces
      if (fespace.GetLocalTDofNumber(ldof) < 0) {
        continue;
      }
      interface_dofs.Append(ldof);
      interface_elems.Append(i);
      interface_eldofs.Append(j);
    }
  }

  // Initialize interface coordinates for preCICE
  std::vector<double> interface_coords;

  // Get the coordinates of each node.
  pmesh.SetCurvature(order); // ensure that nodes are defined for p>1
  const GridFunction &coords = *pmesh.GetNodes();

  // Convert ldofs to vdofs (include indices for y-coords).
  // Indices are appended to end.
  Array<int> interface_vdofs = interface_dofs;
  pmesh.GetNodalFESpace()->DofsToVDofs(interface_vdofs);

  // Sort indices so that x,y indices are adjacent.
  interface_vdofs.Sort();

  // Get the coordinates, correctly organized.
  interface_coords.resize(interface_vdofs.Size());
  coords.GetSubVector(interface_vdofs, interface_coords.data());

  pmesh.SetCurvature(1);

  // Initialize temperature with IC
  ParGridFunction u_gf(&fespace);
  u_gf = 310;

  // Lower + upper boundaries are essential --> Get tdofs.
  Array<int> bdr_marker_arr(pmesh.bdr_attributes.Size());
  bdr_marker_arr                               = 0;
  bdr_marker_arr[pmesh.bdr_attributes.Find(1)] = 1; // lower
  bdr_marker_arr[pmesh.bdr_attributes.Find(3)] = 1; // upper
  Array<int> ess_tdof_list;
  fespace.GetEssentialTrueDofs(bdr_marker_arr, ess_tdof_list);

  // Initialize coefficient for alpha.
  ConstantCoefficient alpha(1.0);

  // Initialize coefficient for thermal conductivity.
  ConstantCoefficient k(100.0);

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
    pvdc = std::make_unique<ParaViewDataCollection>("ParaView", &pmesh);
    pvdc->SetLevelsOfDetail(order);
    pvdc->SetDataFormat(VTKFormat::BINARY);
    pvdc->SetHighOrderOutput(true);
    pvdc->SetCycle(0);
    pvdc->SetTime(0.0);
    pvdc->RegisterField("Temperature", &u_gf);
    pvdc->Save();
  }

  // GLVis parameters.
  socketstream sout;
  const auto   WriteGLVis = [&]() {
    VisualizeField(sout, vishost, visport, u_gf,
                     "flow-over-heated-plate", 10, 10, 600, 300,
                     "ARj******]]]]]]]c", false,
                     "autoscale off\n valuerange 300 310\n");
  };

  if (visualization) {
    WriteGLVis();
  }

  // Initialize preCICE
  const std::string mesh_name = "Solid-Mesh";
  Participant       participant("Solid", precice_config, rank, size);
  std::vector<int>  mesh_vertices(interface_coords.size() / dim);
  participant.setMeshVertices(mesh_name, interface_coords, mesh_vertices);
  participant.initialize();

  // Initialize preCICE-related variables.
  real_t              precice_dt;
  real_t              t_save;
  ParGridFunction     u_gf_save(&fespace);
  std::vector<double> u_receive(interface_dofs.Size());
  std::vector<double> qwall_write(interface_dofs.Size());

  // Main solver loop.
  while (participant.isCouplingOngoing()) {
    if (participant.requiresWritingCheckpoint()) {
      t_save    = t;
      u_gf_save = u_gf;
    }

    precice_dt = participant.getMaxTimeStepSize();
    dt         = std::min(dt, precice_dt);

    // Get temperatures and set.
    participant.readData(mesh_name, "Temperature", mesh_vertices, dt,
                         u_receive);
    u_gf.SetSubVector(interface_dofs, u_receive.data());
    u_gf.SetTrueVector();

    // Step in time
    ode_solver->Step(u_gf.GetTrueVector(), t, dt);
    u_gf.SetFromTrueVector();

    // Get and write the heat fluxes.
    // See https://mfem.org/howto/outer_normals/
    Vector normal(dim), grad_u(dim);
    for (int i = 0; i < interface_elems.Size(); i++) {

      const int be_i  = interface_elems[i];
      const int eldof = interface_eldofs[i];

      const FiniteElement    &fe  = *fespace.GetBE(be_i);
      ElementTransformation  &trf = *fespace.GetBdrElementTransformation(be_i);
      const IntegrationPoint &ip  = fe.GetNodes().IntPoint(eldof);
      trf.SetIntPoint(&ip);

      CalcOrtho(trf.Jacobian(), normal);
      normal /= normal.Norml2();
      u_gf.GetGradient(trf, grad_u);
      qwall_write[i] = -k.Eval(trf, ip) * (grad_u * normal);
    }
    participant.writeData(mesh_name, "Heat-Flux", mesh_vertices, qwall_write);
    participant.advance(dt);

    if (participant.requiresReadingCheckpoint()) {
      t    = t_save;
      u_gf = u_gf_save;
    } else {
      if (visualization) {
        WriteGLVis();
      }

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
    : TimeDependentOperator(f.GetTrueVSize()),
      ess_tdof_list_(ess_tdof_list),
      M_solver_(f.GetComm()),
      A_solver_(f.GetComm()),
      M_(&f),
      K_(&f),
      rhs_(f.GetTrueVSize())
{
  // Double-precision relative tolerance:
  const real_t rel_tol = 1e-12;

  const int max_iter = 100;

  M_.AddDomainIntegrator(new MassIntegrator());
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

void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    const ParGridFunction &gf, const char *title,
                    int x, int y, int w, int h, const char *keys, bool vec,
                    const char *commands)
{
  ParMesh &pmesh = *gf.ParFESpace()->GetParMesh();
  MPI_Comm comm  = pmesh.GetComm();

  int num_procs, myid;
  MPI_Comm_size(comm, &num_procs);
  MPI_Comm_rank(comm, &myid);

  bool newly_opened = false;
  int  connection_failed;

  do {
    if (myid == 0) {
      if (!sock.is_open() || !sock) {
        sock.open(vishost, visport);
        sock.precision(8);
        newly_opened = true;
      }
      sock << "solution\n";
    }

    pmesh.PrintAsOne(sock);
    gf.SaveAsOne(sock);

    if (myid == 0 && newly_opened) {
      if (commands) {
        sock << commands << "\n";
      }
      sock << "window_title '" << title << "'\n"
           << "window_geometry "
           << x << " " << y << " " << w << " " << h << "\n";
      if (keys) {
        sock << "keys " << keys << "\n";
      } else {
        sock << "keys maaAc";
      }
      if (vec) {
        sock << "vvv";
      }
      sock << std::endl;
    }

    if (myid == 0) {
      connection_failed = !sock && !newly_opened;
    }
    MPI_Bcast(&connection_failed, 1, MPI_INT, 0, comm);
  } while (connection_failed);
}
