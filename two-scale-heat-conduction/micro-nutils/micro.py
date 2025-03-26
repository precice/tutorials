"""
Micro simulation
In this script we solve the Laplace equation with a grain depicted by a phase field on a square domain :math:`Ω`
with boundary :math:`Γ`, subject to periodic boundary conditions in both dimensions
"""
import math

from nutils import mesh, function, solver, export, cli
import treelog
import numpy as np
from copy import deepcopy


class MicroSimulation:

    def __init__(self, sim_id):
        """
        Constructor of MicroSimulation class.
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
        output_data["k_00"] = k[0][0]
        output_data["k_11"] = k[1][1]
        output_data["porosity"] = psi

        return output_data

    def _reinitialize_namespace(self, topo):
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
        return 1. / (1. + np.exp(-4. / lam * (np.sqrt(x ** 2 + y ** 2) - r)))

    @staticmethod
    def _get_analytical_phasefield(topo, ns, degree_phi, lam, r):
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
        return [self._solphi.copy(), deepcopy(self._topo)]

    def set_state(self, state):
        self._solphi = state[0]
        self._topo = state[1]
        self._reinitialize_namespace(self._topo)  # The namespace also needs to reloaded to its earlier state

    def _refine_mesh(self, topo_nm1, solphi_nm1):
        """
        At the time of the calling of this function a predicted solution exists in ns.phi
        """
        if not self._initial_condition_is_set:
            solphi = solphi_nm1
            # ----- Refine the coarse mesh according to the projected solution to get a predicted refined topology ----
            topo = self._topo_coarse
            for level in range(self._ref_level):
                # print("refinement level = {}".format(level))
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
            for level in range(self._ref_level):
                # print("refinement level = {}".format(level))
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
        Solving the Allen-Cahn equation using a Newton solver.
        Returns porosity of the micro domain.
        """
        self._first_iter_done = True
        resphi = topo.integral('(lam^2 phibasis_n dphidt + gam phibasis_n ddwpdphi + gam lam^2 phibasis_n,i phi_,i + '
                               '4 lam reacrate phibasis_n phi (1 - phi)) d:x' @ self._ns, degree=self._degree_phi * 2)

        args = dict(solphinm1=phi_coeffs_nm1, dt=dt, conc=concentration)
        phi_coeffs = solver.newton('solphi', resphi, lhs0=phi_coeffs_nm1, arguments=args).solve(tol=1E-12)

        return phi_coeffs

    def _get_avg_porosity(self, topo, phi_coeffs):
        psi = topo.integral('phi d:x' @ self._ns, degree=self._degree_phi * 2).eval(solphi=phi_coeffs)

        return psi

    def _solve_heat_cell_problem(self, topo, phi_coeffs):
        """
        Solving the P1 homogenized heat equation
        Returns upscaled conductivity for the micro domain
        """
        res = topo.integral('((phi ks + (1 - phi) kg) u_i,j ubasis_ni,j - '
                            '(ks - kg) phi_,j $_ij ubasis_ni) d:x' @ self._ns, degree=self._degree_u * 2)

        args = dict(solphi=phi_coeffs)
        u_coeffs = solver.solve_linear('solu', res, constrain=self._ucons, arguments=args)

        return u_coeffs

    def _get_eff_conductivity(self, topo, u_coeffs, phi_coeffs):
        b = topo.integral(
            self._ns.eval_ij('(phi ks + (1 - phi) kg) ($_ij + du_ij) d:x'),
            degree=self._degree_u *
            2).eval(
            solu=u_coeffs,
            solphi=phi_coeffs)

        return b.export("dense")

    def solve(self, macro_data, dt):
        if self._psi_nm1 < 0.95:
            topo, solphi = self._refine_mesh(self._topo, self._solphi)
            self._reinitialize_namespace(topo)

            assert ((solphi >= 0.0) & (solphi <= 1.0)).all()

            solphi = self._solve_allen_cahn(topo, solphi, macro_data["concentration"], dt)
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
        output_data["k_00"] = k[0][0]
        output_data["k_01"] = k[0][1]
        output_data["k_10"] = k[1][0]
        output_data["k_11"] = k[1][1]
        output_data["porosity"] = psi
        output_data["grain_size"] = math.sqrt((1 - psi) / math.pi)

        return output_data


def main():
    micro_problem = MicroSimulation(0)
    dt = 1e-3
    micro_problem.initialize()
    concentrations = [0.5, 0.4]
    t = 0.0
    n = 0
    concentration = dict()

    for conc in concentrations:
        concentration["concentration"] = conc

        micro_sim_output = micro_problem.solve(concentration, dt)

        # micro_problem.output()
        t += dt
        n += 1
        print(micro_sim_output)


if __name__ == "__main__":
    cli.run(main)
