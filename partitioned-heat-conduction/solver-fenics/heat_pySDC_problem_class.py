from fenics import *
import dolfin as df

from pySDC.core.Problem import ptype
from pySDC.implementations.datatype_classes.fenics_mesh import fenics_mesh, rhs_fenics_mesh


class fenics_heat_2d(ptype):
    dtype_u = fenics_mesh
    dtype_f = rhs_fenics_mesh
    
    def __init__(self, mesh, function_space, forcing_term_expr, solution_expr, coupling_boundary, remaining_boundary, coupling_expr, precice_ref, participant_name):
        # Set precice reference and coupling expression reference to update coupling boundary 
        # at every step within pySDC
        self.precice = precice_ref
        self.coupling_expression = coupling_expr
        self.t_start = 0.0
        self.participant_name = participant_name
        
        # set mesh and function space for future reference
        self.mesh = mesh
        self.V = function_space
        
        # invoke super init
        super(fenics_heat_2d, self).__init__(self.V)
        
        # Define Trial and Test function
        u = df.TrialFunction(self.V)
        v = df.TestFunction(self.V)
        
        # Mass term
        a_M = u * v * df.dx
        self.M = df.assemble(a_M)
        
        # Stiffness term (Laplace)
        a_K = -1.0 * df.inner(df.nabla_grad(u), df.nabla_grad(v)) * df.dx
        self.K = df.assemble(a_K)         
        
        # Forcing term
        self.forcing_term_expr = forcing_term_expr
        
        # Solution expression for error comparison and as boundary condition 
        # on the non-coupling boundary
        self.solution_expr = solution_expr
         
        # Currently only for Dirichlet boundary, has to be changed for Neumann boundary
        if self.participant_name == 'Dirichlet':
            self.couplingBC = df.DirichletBC(self.V, coupling_expr, coupling_boundary)
            
        self.remainingBC = df.DirichletBC(self.V, solution_expr, remaining_boundary)
        
        # Allow for fixing the boundary conditions for the residual computation
        # Necessary if imex-1st-order-mass sweeper is used
        self.fix_bc_for_residual = True
        
        # define the homogeneous Dirichlet boundary for residual correction
        def FullBoundary(x, on_boundary):
            return on_boundary
        self.homogenousBC = df.DirichletBC(self.V, df.Constant(0), FullBoundary)
        
        
    def solve_system(self, rhs, factor, u0, t):
        r"""
        Solves (M - factor \cdot A) \vec{u} = \vec{rhs}.
        """
        u = self.dtype_u(u0)
        T = self.M - factor * self.K
        b = self.dtype_u(rhs)

        self.solution_expr.t = t
        self.remainingBC.apply(T, b.values.vector())
         
        if self.participant_name == 'Dirichlet':     
            # Coupling BC needs to point to correct time
            dt = t - self.t_start
            read_data = self.precice.read_data(dt)
            self.precice.update_coupling_expression(self.coupling_expression, read_data)   
            
            self.couplingBC.apply(T, b.values.vector())

        df.solve(T, u.values.vector(), b.values.vector())

        return u

    def eval_f(self, u, t):
        """
            Derivative computation.
            Returns a dtype_f object (rhs_fenics_mesh in our case),
            with an explicit and implicit part.
        """
        f = self.dtype_f(self.V)

        # Implicit part: K*u
        self.K.mult(u.values.vector(), f.impl.values.vector())
        
        # Explicit part: M*g   (g = forcing term)
        self.forcing_term_expr.t = t
        f.expl = self.dtype_u(df.interpolate(self.forcing_term_expr, self.V))
        f.expl = self.apply_mass_matrix(f.expl)

        return f
    
    def fix_residual(self, res):
        """
            Applies homogeneous Dirichlet boundary conditions to the residual
        """
        self.homogenousBC.apply(res.values.vector())
        return None
    
    def apply_mass_matrix(self, u):
        r"""
            Returns M*u.
        """

        me = self.dtype_u(self.V)
        self.M.mult(u.values.vector(), me.values.vector())

        return me

    def u_exact(self, t):
        r"""
            Returns the exact solution at time t.
        """
        self.solution_expr.t = t
        
        me = self.dtype_u(interpolate(self.solution_expr, self.V), val=self.V)
        return me

    def set_t_start(self, t_start):
        self.t_start = t_start