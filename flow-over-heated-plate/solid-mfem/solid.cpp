#include <mfem/mfem.hpp>
#include <precice/precice.hpp>

using namespace mfem;
using namespace precice;


class LinearDiffusionOperator : public TimeDependentOperator
{
protected:
   ParFiniteElementSpace &f_

   const Array<int> &ess_tdof_list_;

   // Mass bilinear form
   ParBilinearForm M_;

   // Mass matrix
   std::unique_ptr<HypreParMatrix> M_mat_;

   // Stiffness bilinear form.
   ParBilinearForm K_;

   // Stiffness matrix.
   std::unique_ptr<HypreParMatrix> K_mat_;

   // Neumann linear form.
   ParLinearForm b_;

public:
   LinearDiffusionOperator(ParFiniteElementSpace &f, 
                           const Array<int> &ess_tdof_list,
                           Coefficient &alpha,
                           Coefficient &q_bdr);
   
   void Mult(const Vector &u, Vector &du_dt) const override;

   void ImplicitSolve(const real_t dt, const Vector &u, Vector &k) override;

   void ReassembleNeumann(Coefficient &q_bdr_);
};

int main(int argc, char *argv[])
{
   mfem::Mpi::Init(argc, argv);
   int size = mfem::Mpi::WorldSize();
   int rank = mfem::Mpi::WorldRank();
   mfem::Hypre::Init();

   std::string precice_config = "../precice-config.xml";
   bool visualization = true;
   int visport = 19916;
   char vishost[] = "localhost";
   
   OptionsParser args(argc, argv);
   args.AddOption(&precice_config, "-c", "--config-file", 
                  "preCICE configuration file.")
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&visport, "-p", "--send-port", "Socket for GLVis.");
   args.Parse();
   if (!args.Good())
   {
      if (rank == 0)
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (rank == 0)
   {
      args.PrintOptions(cout);
   }

   // Create mesh. Boundary attribute 3 = interface, 1 = lower boundary.
   const int nx = 100;
   const int ny = 25;
   const Element e_type = Element::QUADRILATERAL;
   const bool gen_edges = true;
   const real_t Lx = 1.0;
   const real_t Ly = 0.25;
   Mesh mesh = Mesh::MakeCartesian2D(nx, ny, e_type, gen_edges, Lx, Ly);
   mesh.Transform(
   [](const Vector &x_old, Vector &x_new)
   {
      // Ensure top of plate at y = 0.
      x_new[1] = x_old[1] - Ly;
   });

   // Partition grid.
   ParMesh pmesh(MPI_COMM_WORLD, &mesh);

   // Use order p = 1, 2D H1-conforming finite elements.
   const int p = 1;
   const int dim = 2;
   H1_FECollection fecoll(p, dim);

   Vector interface_ldofs;
   Vector lcoords;
   // Get the coordinates of nodes on this rank.
   {
      // Create an FE space for the obtaining the nodal coordinates (2D).
      // Use Ordering::byVDIM to ensure XYZ XYZ ...
      ParFiniteElementSpace fes_nodes(&mesh, &fecoll, 2, Ordering::byVDIM);

      // Get the interface coordinates.
      // (If using p > 1, SetNodalFESpace() is necessary).
      pmesh.SetNodalFESpace(&fespace);
      const GridFunction &coords = *pmesh.GetNodes();

      // Get the interface DOFs.
      Array<int> interface_bdr_elems, interface_bdr_dofs;
      fes_nodes.GetBoundaryElementsByAttribute(3, interface_bdr_elems);
      fes_nodes.FiniteElementSpace::GetBoundaryLoopEdgeDofs(
                                       interface_bdr_elems,
                                       interface_bdr_dofs);
      // Note that this will yield duplicate mesh vertices across ranks.                                       
      // See ParFiniteElementSpace::GetBoundaryLoopEdgeDofs to get unique
      //    vertices across all ranks.                                     
      interface_ldofs = interface_bdr_dofs;

      // Convert ldofs to vdofs (include indices for y-coords).
      // Indices are appended to end.
      fes_nodes.DofsToVDofs(interface_bdr_dofs);

      // Sort indices so that x,y indices are adjacent.
      interface_bdr_dofs.Sort();

      // Get the coordinates, correctly organized.
      coords.GetSubVector(interface_bdr_dofs, lcoords);
   }

   // Initialize FE space for temperature solution.
   // Use same ordering to ensure indices `interface_ldofs` match.
   ParFiniteElementSpace fespace(&pmesh, &fecoll, ordering=Ordering::byVDIM);
   
   // Initialize temperature with IC
   ParGridFunction temperature(&fespace);
   temperature = 300;

   // Lower boundary will be essential, fixed at T=300 K.
   // Get tdofs.
   Array<int> bdr_marker_arr(pmesh.bdr_attributes.Size());
   bdr_marker_arr = 0;
   bdr_marker_arr[bdr_marker_arr.Find(1)] = 1;
   Array<int> ess_tdof_list;
   fespace.GetEssentialTrueDofs(bdr_marker_arr, ess_tdof_list);

   // Initialize coefficient for alpha.
   ConstantCoefficient alpha(1.0);

   // Initialize interface boundary heat flux coefficient.
   // We use a GridFunction to store q_wall.
   ParGridFunction q_interface(&fespace);
   GridFunctionCoefficient q_interface_coeff(&q_wall);

   // Initialize Neumann term coefficient.
   // Only add coefficient for interface.
   PWCoefficient q_bdr;
   q_bdr.UpdateCoefficient(3, q_interface_coeff);

   // Initialize the operator.
   const real_t alpha = 1.0;
   LinearDiffusionOperator oper(fespace, ess_tdof_list, alpha, q_bdr);

   // Initialize preCICE
   Participant participant("Solid", "../precice-config.xml", rank, size);
   Array<int> mesh_vertices(lcoords/dim);
   participant.setMeshVertices("Solid", lcoords, mesh_vertices);



   return 0;
}


LinearDiffusionOperator::LinearDiffusionOperator(
   ParFiniteElementSpace &f,
   const Array<int> &ess_tdof_list,
   Coefficient &alpha,
   Coefficient &q_bdr)
: f_(f),
  ess_tdof_list_(ess_tdof_list),
  M_(&f),
  K_(&f),
  b_(&f),
  alpha_(alpha)
{
   M_.AddDomainIntegrator(new MassIntegrator);
   M_.Assemble(0);
   M_.Finalize(0);
   M_mat_ = std::make_unique<HypreParMatrix>(M.ParallelAssemble());
   M_mat_->EliminateBC(ess_tdof_list, Operator::DiagonalPolicy::DIAG_ONE);

   K_.AddDomainIntegrator(new DiffusionIntegrator(alpha_));
   K_.Assemble(0);
   K_.Finalize(0);
   K_mat_ = std::make_unique<HypreParMatrix>(K_.ParallelAssemble());
   // Do not eliminate BCs of K, as it is on RHS and performed after evaluation.
}