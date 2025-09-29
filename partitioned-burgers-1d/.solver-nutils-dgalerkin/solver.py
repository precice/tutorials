# Original source: https://github.com/evalf/nutils/blob/d73749ff7d64c9ccafdbb88cd442f80b9448c118/examples/burgers.py

from nutils import mesh, function, export, testing
from nutils.solver import System
from nutils.expression_v2 import Namespace
import treelog as log
import numpy as np
import itertools
import precice
import json
import os
import argparse

def _generate_initial_condition(x_coords, ic_config, epoch):
    np.random.seed(epoch)
    ic_values = np.zeros(len(x_coords))
    if ic_config["type"] == "sinusoidal":
        num_modes = ic_config.get("num_modes", 1)
        superpositions = np.random.randint(2, num_modes + 1)
        for _ in range(superpositions):
            amp = np.random.uniform(0.1, 2)
            k = np.random.randint(ic_config["wavenumber_range"][0], ic_config["wavenumber_range"][1] + 1)
            phase_shift = np.random.uniform(0, 2 * np.pi)
            ic_values += amp * np.sin(2 * np.pi * k * x_coords + phase_shift)
    return ic_values

def project_initial_condition(domain_min, domain_max, nelems, ic_config, epoch):
    # 1. Generate a high-resolution "truth" on a fine grid
    fine_res = nelems * 10
    fine_x = np.linspace(domain_min[0], domain_max[0], fine_res, endpoint=False)
    fine_u = _generate_initial_condition(fine_x, ic_config, epoch)

    # 2. Average the high-resolution truth over each coarse cell
    u_projected = np.zeros(nelems)
    for i in range(nelems):
        cell_start = i * 10
        cell_end = (i + 1) * 10
        u_projected[i] = np.mean(fine_u[cell_start:cell_end])
        
    return u_projected

def main(dim: int,
         epoch: int = 0,
         btype: str = 'discont',
         degree: int = 1,
         newtontol: float = 1e-5,
         config_file: str = "precice-config.xml"):

    script_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(script_dir, "..", "python_participant", "config.json"), 'r') as f:
        config = json.load(f)["solver"]
    with open(os.path.join(script_dir, "ic_params.json"), 'r') as f:
        ic_config = json.load(f)["initial_conditions"]
    
    config_path = os.path.join(script_dir, "..", f"{dim}d", config_file)
    participant = precice.Participant("Solver", config_path, 0, 1)
    
    mesh_internal_name = f"Solver-Mesh-{dim}D-Internal"
    mesh_boundaries_name = f"Solver-Mesh-{dim}D-Boundaries"
    data_name = f"Data_{dim}D"
    res = config[f"{dim}d_resolution"]
    domain_min = config[f"{dim}d_domain_min"]
    domain_max = config[f"{dim}d_domain_max"]
    nelems = res[0]

    domain, geom = mesh.line(np.linspace(domain_min[0], domain_max[0], nelems + 1), periodic=True)
    
    # Define all nelems +1 nodes for evaluation
    eval_coords_x = np.linspace(domain_min[0], domain_max[0], nelems + 1)
    
    # Define the nelems vertices for saving (all but the last)
    trunc_coords_x = eval_coords_x[:-1]
    internal_coords = np.array([trunc_coords_x, np.full(len(trunc_coords_x), domain_min[1])]).T
    boundary_coords = np.array([[domain_min[0], domain_min[1]], [domain_max[0], domain_max[1]]])
    
    internal_vertex_ids = participant.set_mesh_vertices(mesh_internal_name, internal_coords)
    boundary_vertex_ids = participant.set_mesh_vertices(mesh_boundaries_name, boundary_coords)

    sample = domain.locate(geom, eval_coords_x, tol=1e-5)
    
    ns = Namespace()
    ns.x = geom
    ns.define_for('x', gradient='∇', normal='n', jacobians=('dV', 'dS'))
    ns.u = domain.field('u', btype=btype, degree=degree)
    ns.du = ns.u - function.replace_arguments(ns.u, 'u:u0')
    ns.v = domain.field('v', btype=btype, degree=degree)
    ns.t = function.field('t')
    ns.dt = ns.t - function.field('t0')
    ns.f = '.5 u^2'
    ns.C = 1

    res_pde = domain.integral('(v du / dt - ∇(v) f) dV' @ ns, degree=degree*2)
    res_pde -= domain.interfaces.integral('[v] n ({f} - .5 C [u] n) dS' @ ns, degree=degree*2)
    system = System(res_pde, trial='u', test='v')

    # Project the initial condition
    u_averaged = project_initial_condition(domain_min, domain_max, nelems, ic_config, epoch)
    ns.uic = domain.basis('discont', degree=0).dot(u_averaged)

    sqr = domain.integral('(u - uic)^2 dV' @ ns, degree=max(degree*2, 5))
    args = System(sqr, trial='u').solve()

    if participant.requires_initial_data():
        # Evaluate at all nodes
        all_data = sample.eval(ns.u, arguments=args)
        # Truncate last element
        trunc_data = all_data[:-1]
        boundary_data_values = np.array([trunc_data[0], trunc_data[-1]])
        
        participant.write_data(mesh_internal_name, data_name, internal_vertex_ids, trunc_data)
        participant.write_data(mesh_boundaries_name, data_name, boundary_vertex_ids, boundary_data_values)

    participant.initialize()

    args['t'] = 0.

    with log.iter.plain('timestep', itertools.count()) as steps:
        for _ in steps:
            if not participant.is_coupling_ongoing():
                break

            timestep = participant.get_max_time_step_size()

            args = system.step(timestep=timestep, arguments=args, timearg='t', suffix='0', tol=newtontol)

            all_data = sample.eval(ns.u, arguments=args)
            trunc_data = all_data[:-1]
            boundary_data_values = np.array([trunc_data[0], trunc_data[-1]])

            participant.write_data(mesh_internal_name, data_name, internal_vertex_ids, trunc_data)
            participant.write_data(mesh_boundaries_name, data_name, boundary_vertex_ids, boundary_data_values)

            participant.advance(timestep)

    participant.finalize()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("dim", type=int, choices=[1, 2], help="Dimension of the simulation")
    parser.add_argument('--config_file', type=str, default="precice-config.xml")
    parser.add_argument('--btype', type=str, default='discont')
    parser.add_argument('--degree', type=int, default=1)
    parser.add_argument('--newtontol', type=float, default=1e-5)
    parser.add_argument("--epoch", type=int, default=0, help="Current epoch number")

    args_cli = parser.parse_args()

    main(dim=args_cli.dim, epoch=args_cli.epoch, btype=args_cli.btype, degree=args_cli.degree, newtontol=args_cli.newtontol, config_file=args_cli.config_file)
