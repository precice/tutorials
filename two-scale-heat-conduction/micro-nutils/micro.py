"""
Micro simulation
In this script we solve the Laplace equation with a grain depicted by a phase field on a square domain :math:`Ω`
with boundary :math:`Γ`, subject to periodic boundary conditions in both dimensions
"""
import math

from nutils import mesh, function, solver, cli
# import treelog, export
import numpy as np
from copy import deepcopy


class MicroSimulation:

    def __init__(self, sim_id):
        """
        Constructor of MicroSimulation class.

        Parameters
        ----------
        sim_id : int
            Unique global ID assigned to this micro simulation instance by the Micro Manager.
        """
        self._sim_id = sim_id

        # Initial parameters
        self._nelems = 6  # Elements in one direction

        self._ref_level = 3  # Number of levels of mesh refinement
        self._r_initial = 0.4  # Initial radius of the grain

        # Interpolation order for phi and u
        self._degree_phi = 2
        self._degree_u = 2

        # Set up mesh with periodicity in both X and Y directions
        self._topo, self._geom = mesh.rectilinear([np.linspace(-0.5, 0.5, self._nelems + 1)] * 2, periodic=(0, 1))
        self._topo_coarse = self._topo  # Save original coarse topology to use to re-refinement

        self._solu = None  # Solution of weights for which cell problem is solved for
        self._solphi = None  # Solution of phase field
        self._solphinm1 = None  # Solution of phase field at t_{n-1}
        self._solphi_checkpoint = None  # Save the current solution of the phase field as a checkpoint.
        self._topo_checkpoint = None  # Save the refined mesh as a checkpoint.
        self._ucons = None
        self._first_iter_done = False
        self._initial_condition_is_set = False
        self._k_nm1 = None  # Average effective conductivity of last time step

    def initialize(self):
        """
        Set up the initial phase field and solve the heat cell problem on a refined mesh.

        Computes the analytical initial grain geometry, refines the mesh around the
        grain interface, and solves the cell problem to obtain initial effective
        conductivity and porosity values.

        Returns
        -------
        dict
            Initial output data with keys ``K00``, ``K11``, and ``Porosity``.
        """
        # Define initial namespace
        self._ns = function.Namespace()
        self._ns.x = self._geom

        self._ns.phibasis = self._topo.basis('std', degree=self._degree_phi)
        self._ns.phi = 'phibasis_n ?solphi_n'  # Initial phase field
        self._ns.coarsephibasis = self._topo_coarse.basis('std', degree=self._degree_phi)
        self._ns.coarsephi = 'coarsephibasis_n ?coarsesolphi_n'  # Phase field on original coarse topology
        self._ns.lam = (4 / self._nelems) / (2 ** self._ref_level)
        self._ns.coarselam = 3 / self._nelems

        # Initialize phase field
        solphi = self._get_analytical_phasefield(
            self._topo,
            self._ns,
            self._degree_phi,
            self._ns.coarselam,
            self._r_initial)

        # Refine the mesh
        self._topo, self._solphi = self._refine_mesh(self._topo, solphi)
        self._reinitialize_namespace(self._topo)
        self._initial_condition_is_set = True

        # Initialize phase field once more on refined topology
        solphi = self._get_analytical_phasefield(self._topo, self._ns, self._degree_phi, self._ns.lam, self._r_initial)
        self._solphi = solphi  # Save solution of phi
        psi = self._get_avg_porosity(self._topo, solphi)
        self._psi_nm1 = psi  # Average porosity value of last time step

        # Solve the heat cell problem
        solu = self._solve_heat_cell_problem(self._topo, solphi)
        k = self._get_eff_conductivity(self._topo, solu, solphi)

        self._solu = solu  # Save solution for output

        output_data = dict()
        output_data["K00"] = k[0][0]
        output_data["K11"] = k[1][1]
        output_data["Porosity"] = psi

        return output_data

    def _reinitialize_namespace(self, topo):
        """
        Rebuild the nutils function namespace for a given (possibly refined) topology.

        Called whenever the mesh changes (after each refinement step) to keep all
        basis functions and symbolic expressions consistent with the current topology.

        Parameters
        ----------
        topo : nutils.topology.Topology
            The topology (coarse or hierarchically refined) to build the namespace for.
        """
        self._ns = None  # Clear old namespace
        self._ns = function.Namespace()
        self._ns.x = self._geom
        self._ns.ubasis = topo.basis('h-std', degree=self._degree_u).vector(topo.ndims)
        self._ns.phibasis = topo.basis('h-std', degree=self._degree_phi)
        self._ns.coarsephibasis = self._topo_coarse.basis('std', degree=self._degree_phi)

        self._ns.lam = (4 / self._nelems) / (2 ** self._ref_level)  # Diffuse interface width
        self._ns.gam = 0.05
        self._ns.kt = 1.0
        self._ns.eqconc = 0.5  # Equilibrium concentration
        self._ns.kg = 0.0  # Conductivity of grain material
        self._ns.ks = 1.0  # Conductivity of sand material
        self._ns.reacrate = 'kt (?conc / eqconc)^2 - 1'  # Constructed reaction rate based on macro temperature
        self._ns.u = 'ubasis_ni ?solu_n'  # Weights for which cell problem is solved for
        self._ns.du_ij = 'u_i,j'  # Gradient of weights field
        self._ns.phi = 'phibasis_n ?solphi_n'  # Phase field
        self._ns.coarsephi = 'coarsephibasis_n ?coarsesolphi_n'  # Phase field on original coarse topology
        self._ns.ddwpdphi = '16 phi (1 - phi) (1 - 2 phi)'  # gradient of double-well potential
        self._ns.dphidt = 'phibasis_n (?solphi_n - ?solphinm1_n) / ?dt'  # Implicit time evolution of phase field

        self._ucons = np.zeros(len(self._ns.ubasis), dtype=bool)
        self._ucons[-1] = True  # constrain u to zero at a point

    @staticmethod
    def _analytical_phasefield(x, y, r, lam):
        """
        Evaluate the analytical diffuse-interface phase field at coordinates (x, y).

        Uses a smooth step function centred on the circle of radius ``r``, with
        interface width controlled by ``lam``.

        Parameters
        ----------
        x : array_like
            X coordinates.
        y : array_like
            Y coordinates.
        r : float
            Radius of the grain.
        lam : float
            Diffuse interface width.

        Returns
        -------
        array_like
            Phase field value in [0, 1]: 1 inside the grain, 0 outside.
        """
        return 1. / (1. + np.exp(-4. / lam * (np.sqrt(x ** 2 + y ** 2) - r)))

    @staticmethod
    def _get_analytical_phasefield(topo, ns, degree_phi, lam, r):
        """
        Project the analytical phase field onto the FE basis of the given topology.

        Minimises the L2 error between the analytical expression and the discrete
        representation to obtain DOF coefficients.

        Parameters
        ----------
        topo : nutils.topology.Topology
            Topology on which the projection is performed.
        ns : nutils.function.Namespace
            Namespace containing the ``phi`` and ``x`` symbols.
        degree_phi : int
            Polynomial degree of the phase field basis.
        lam : float
            Diffuse interface width passed to the analytical formula.
        r : float
            Grain radius passed to the analytical formula.

        Returns
        -------
        numpy.ndarray
            DOF coefficient vector for the projected phase field.
        """
        phi_ini = MicroSimulation._analytical_phasefield(ns.x[0], ns.x[1], r, lam)
        sqrphi = topo.integral((ns.phi - phi_ini) ** 2, degree=degree_phi * 2)
        solphi = solver.optimize('solphi', sqrphi, droptol=1E-12)

        return solphi

    # def output(self):
    #  bezier = self._topo.sample('bezier', 2)
    #  x, u, phi = bezier.eval(['x_i', 'u_i', 'phi'] @ self._ns, solu=self._solu, solphi=self._solphi)
    #  with treelog.add(treelog.DataLog()):
    #      export.vtk("micro-heat-{}".format(self._sim_id), bezier.tri, x, T=u, phi=phi)

    def get_state(self):
        """
        Return the current simulation state for checkpointing or state transfer.

        Returns
        -------
        list
            ``[solphi, topo, psi_nm1, k_nm1]`` — phase field DOF vector (copy),
            the current refined topology (deep copy), average porosity, and
            effective conductivity tensor from the last time step.
        """
        return [self._solphi.copy(), deepcopy(self._topo), self._psi_nm1, self._k_nm1]

    def set_state(self, state):
        """
        Restore the simulation state from a checkpoint or a donor simulation.

        Also rebuilds the nutils namespace for the restored topology and marks
        the initial condition as set so that subsequent :meth:`solve` calls use
        the correct (projection-based) mesh-refinement branch.

        Parameters
        ----------
        state : list
            State list as returned by :meth:`get_state`:
            ``[solphi, topo, psi_nm1, k_nm1]``.
        """
        self._solphi = state[0]
        self._topo = state[1]
        self._psi_nm1 = state[2]
        self._k_nm1 = state[3]
        self._initial_condition_is_set = True  # state comes from a fully initialized sim, use the 'else' branch in _refine_mesh
        self._reinitialize_namespace(self._topo)  # The namespace also needs to reloaded to its earlier state

    def _refine_mesh(self, topo_nm1, solphi_nm1):
        """
        Adaptively refine the mesh around the grain interface.

        Two code paths exist:

        * **Initial refinement** (``_initial_condition_is_set = False``): starts
          from the coarse topology and uses the analytical phase field projected
          onto the coarse basis (``coarsesolphi``) to drive refinement.
        * **Subsequent refinement** (``_initial_condition_is_set = True``): starts
          from the coarse topology but guides refinement with the previous
          refined solution (``solphi``), then projects it onto the new refined mesh.

        Parameters
        ----------
        topo_nm1 : nutils.topology.Topology
            Refined topology from the previous time step (used in the projection step).
        solphi_nm1 : numpy.ndarray
            Phase field DOF vector on ``topo_nm1``.

        Returns
        -------
        topo : nutils.topology.HierarchicalTopology
            Newly refined topology.
        solphi : numpy.ndarray
            Phase field DOF vector expressed on the new ``topo``.
        """
        if not self._initial_condition_is_set:
            solphi = solphi_nm1
            # ----- Refine the coarse mesh according to the projected solution to get a predicted refined topology ----
            topo = self._topo_coarse
            for _ in range(self._ref_level):
                smpl = topo.sample('uniform', 5)
                ielem, criterion = smpl.eval([topo.f_index, abs(self._ns.coarsephi - .5) < .4], coarsesolphi=solphi)

                # Refine the elements for which at least one point tests true.
                topo = topo.refined_by(np.unique(ielem[criterion]))

                self._reinitialize_namespace(topo)

                phi_ini = MicroSimulation._analytical_phasefield(self._ns.x[0], self._ns.x[1], self._r_initial,
                                                                 self._ns.lam)
                sqrphi = topo.integral((self._ns.coarsephi - phi_ini) ** 2, degree=self._degree_phi * 2)
                solphi = solver.optimize('coarsesolphi', sqrphi, droptol=1E-12)
            # ----------------------------------------------------------------------------------------------------
        else:
            # ----- Refine the coarse mesh according to the projected solution to get a predicted refined topology ----
            topo = self._topo_coarse
            for _ in range(self._ref_level):
                topo_union1 = topo_nm1 & topo
                smpl = topo_union1.sample('uniform', 5)
                ielem, criterion = smpl.eval([topo.f_index, abs(self._ns.phi - .5) < .4], solphi=solphi_nm1)

                # Refine the elements for which at least one point tests true.
                topo = topo.refined_by(np.unique(ielem[criterion]))
            # ----------------------------------------------------------------------------------------------------

            # Create a new projection mesh which is the union of the previous refined mesh and the predicted mesh
            topo_union = topo_nm1 & topo

            # ----- Project the solution of the last time step on the projection mesh -----
            self._ns.projectedphi = function.dotarg('projectedsolphi', topo.basis('h-std', degree=self._degree_phi))
            sqrphi = topo_union.integral((self._ns.projectedphi - self._ns.phi) ** 2, degree=self._degree_phi * 2)
            solphi = solver.optimize('projectedsolphi', sqrphi, droptol=1E-12, arguments=dict(solphi=solphi_nm1))

        # Clip values of phase field to be in [0 1] to avoid overshoots and undershoots due to projection
        solphi = np.clip(solphi, 0, 1)

        return topo, solphi

    def _solve_allen_cahn(self, topo, phi_coeffs_nm1, concentration, dt):
        """
        Solve the Allen-Cahn phase field equation for one time step.

        Uses a fully implicit Newton iteration.  The reaction term is driven by
        the macro-scale concentration field.

        Parameters
        ----------
        topo : nutils.topology.Topology
            Current refined topology.
        phi_coeffs_nm1 : numpy.ndarray
            Phase field DOF vector at the previous time step (initial guess).
        concentration : float
            Macro-scale concentration value at this micro simulation location.
        dt : float
            Time step size.

        Returns
        -------
        numpy.ndarray
            Updated phase field DOF vector at the current time step.
        """
        self._first_iter_done = True
        resphi = topo.integral('(lam^2 phibasis_n dphidt + gam phibasis_n ddwpdphi + gam lam^2 phibasis_n,i phi_,i + '
                               '4 lam reacrate phibasis_n phi (1 - phi)) d:x' @ self._ns, degree=self._degree_phi * 2)

        args = dict(solphinm1=phi_coeffs_nm1, dt=dt, conc=concentration)
        phi_coeffs = solver.newton('solphi', resphi, lhs0=phi_coeffs_nm1, arguments=args).solve(tol=1E-12)

        return phi_coeffs

    def _get_avg_porosity(self, topo, phi_coeffs):
        """
        Compute the volume-averaged porosity of the micro domain.

        Porosity is defined as the integral of the phase field ``phi`` over the
        unit cell (phi = 1 inside the grain, 0 outside).

        Parameters
        ----------
        topo : nutils.topology.Topology
            Current topology.
        phi_coeffs : numpy.ndarray
            Phase field DOF vector.

        Returns
        -------
        float
            Average porosity value in [0, 1].
        """
        psi = topo.integral('phi d:x' @ self._ns, degree=self._degree_phi * 2).eval(solphi=phi_coeffs)

        return psi

    def _solve_heat_cell_problem(self, topo, phi_coeffs):
        """
        Solve the periodic cell problem for heat conduction homogenisation.

        Computes the corrector field ``u`` by solving the linear cell problem
        arising from first-order homogenisation of the heat equation on the
        micro domain with a phase-field-dependent conductivity.

        Parameters
        ----------
        topo : nutils.topology.Topology
            Current topology.
        phi_coeffs : numpy.ndarray
            Phase field DOF vector.

        Returns
        -------
        numpy.ndarray
            DOF vector for the corrector field ``u``.
        """
        res = topo.integral('((phi ks + (1 - phi) kg) u_i,j ubasis_ni,j - '
                            '(ks - kg) phi_,j $_ij ubasis_ni) d:x' @ self._ns, degree=self._degree_u * 2)

        args = dict(solphi=phi_coeffs)
        u_coeffs = solver.solve_linear('solu', res, constrain=self._ucons, arguments=args)

        return u_coeffs

    def _get_eff_conductivity(self, topo, u_coeffs, phi_coeffs):
        """
        Compute the effective (upscaled) conductivity tensor.

        Evaluates the homogenised conductivity by integrating the cell-problem
        solution ``u`` together with the phase-field-weighted local conductivity
        over the micro domain.

        Parameters
        ----------
        topo : nutils.topology.Topology
            Current topology.
        u_coeffs : numpy.ndarray
            Corrector field DOF vector from :meth:`_solve_heat_cell_problem`.
        phi_coeffs : numpy.ndarray
            Phase field DOF vector.

        Returns
        -------
        numpy.ndarray
            2x2 effective conductivity tensor as a dense array.
        """
        b = topo.integral(
            self._ns.eval_ij('(phi ks + (1 - phi) kg) ($_ij + du_ij) d:x'),
            degree=self._degree_u *
            2).eval(
            solu=u_coeffs,
            solphi=phi_coeffs)

        return b.export("dense")

    def solve(self, macro_data, dt):
        """
        Advance the micro simulation by one time step.

        Refines the mesh, evolves the phase field via the Allen-Cahn equation,
        and solves the cell problem to obtain updated effective conductivity and
        porosity.  If the porosity has reached 0.95 (grain fully dissolved), the
        simulation is frozen and the last computed values are returned.

        Parameters
        ----------
        macro_data : dict
            Dictionary containing macro-scale input data.  Must include the key
            ``"Concentration"`` with the macro concentration at this location.
        dt : float
            Time step size.

        Returns
        -------
        dict
            Output data with keys ``K00``, ``K11``, ``Porosity``, and
            ``Grain-Size``.
        """
        if self._psi_nm1 < 0.95:
            topo, solphi = self._refine_mesh(self._topo, self._solphi)
            self._reinitialize_namespace(topo)

            assert ((solphi >= 0.0) & (solphi <= 1.0)).all()

            solphi = self._solve_allen_cahn(topo, solphi, macro_data["Concentration"], dt)
            psi = self._get_avg_porosity(topo, solphi)

            solu = self._solve_heat_cell_problem(topo, solphi)
            k = self._get_eff_conductivity(topo, solu, solphi)

            # Save state variables
            self._topo = topo
            self._solphi = solphi
            self._solu = solu
            self._psi_nm1 = psi
            self._k_nm1 = k
        else:
            # Micro simulation has reached max porosity limit and hence is not solved
            k = self._k_nm1
            psi = self._psi_nm1

        output_data = dict()
        output_data["K00"] = k[0][0]
        output_data["K11"] = k[1][1]

        # Off-diagonal terms are not sent back as they are zero in this case
        # output_data["K01"] = k[0][1]
        # output_data["K10"] = k[1][0]

        output_data["Porosity"] = psi
        output_data["Grain-Size"] = math.sqrt((1 - psi) / math.pi)

        return output_data


def main():
    """
    Standalone driver for a single micro simulation (used for testing / debugging).

    Creates one :class:`MicroSimulation` instance, initialises it, and runs two
    solve steps with prescribed concentration values.
    """
    micro_problem = MicroSimulation(0)
    dt = 1e-3
    micro_problem.initialize()
    concentrations = [0.5, 0.4]
    t = 0.0
    n = 0
    concentration = dict()

    for conc in concentrations:
        concentration["Concentration"] = conc

        micro_sim_output = micro_problem.solve(concentration, dt)

        # micro_problem.output()
        t += dt
        n += 1
        print(micro_sim_output)


if __name__ == "__main__":
    cli.run(main)
