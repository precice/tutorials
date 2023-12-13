use nalgebra as na;
use std::env;
use std::process::ExitCode;

use precice;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path;

fn write_vtk(
    file: &str,
    dx: f64,
    velocity: &[f64],
    pressure: &[f64],
    cross_section_length: &[f64],
) -> std::io::Result<()> {
    let filepath = path::Path::new(file);
    if !filepath.parent().unwrap().exists() {
        fs::create_dir(filepath.parent().unwrap())?;
    }

    // let mut f = File::create(file).expect("Unable to open vtk file");
    let file = File::create(file).expect("Unable to open vtk file");
    let mut f = BufWriter::new(file);

    let n_points = velocity.len();

    write!(f, "# vtk DataFile Version 2.0\n\n")?;
    write!(f, "ASCII\n\n")?;
    write!(f, "DATASET UNSTRUCTURED_GRID\n\n")?;
    write!(f, "POINTS {n_points} float\n\n")?;
    for i in 0..n_points {
        write!(
            f,
            "{:.16e} 0.0000000000000000e+00 0.0000000000000000e+00\n",
            i as f64 * dx
        )?;
    }
    writeln!(f)?;

    write!(f, "POINT_DATA {n_points}\n\n")?;

    write!(f, "VECTORS velocity float\n")?;
    for v in velocity {
        write!(
            f,
            "{:.16e} 0.0000000000000000e+00 0.0000000000000000e+00\n",
            v
        )?;
    }
    writeln!(f)?;

    write!(f, "SCALARS pressure float\n")?;
    write!(f, "LOOKUP_TABLE default\n")?;
    for v in pressure {
        write!(f, "{:.16e}\n", v)?;
    }
    writeln!(f)?;

    write!(f, "SCALARS diameter float\n")?;
    write!(f, "LOOKUP_TABLE default\n")?;
    for v in cross_section_length {
        write!(f, "{:.16e}\n", v)?;
    }
    writeln!(f)?;
    Ok(())
}

fn fluid_compute_solution(
    velocity_old: &[f64],
    pressure_old: &[f64],
    cross_section_length_old: &[f64],
    cross_section_length: &[f64],
    t: f64,
    n: usize,
    kappa: f64,
    tau: f64,
    velocity: &mut [f64],
    pressure: &mut [f64],
) {
    // Initial guess
    pressure.copy_from_slice(&pressure_old[..]);

    const E: f64 = 10000_f64;
    //const C_MK2 : f64 = E / std::f64::consts::FRAC_2_SQRT_PI;
    let c_mk2: f64 = E / 2.0 * std::f64::consts::PI.sqrt();

    // lhs = 2*N+2

    const ALPHA: f64 = 0.0;
    const L: f64 = 10.0;
    let dx = L / kappa;

    let equations = 2 * n + 2;

    let mut k = 0;
    loop {
        k += 1;

        let res = {
            let mut res = na::DVector::from_element(equations, 0.0);
            for i in 1..n {
                let [p_l, p_c, p_r]: [f64; 3] = pressure[i - 1..i + 2].try_into().unwrap();
                let [v_l, v_c, v_r]: [f64; 3] = velocity[i - 1..i + 2].try_into().unwrap();
                let [c_l, c_c, c_r]: [f64; 3] =
                    cross_section_length[i - 1..i + 2].try_into().unwrap();

                /* Momentum */
                res[i] = (velocity_old[i] * cross_section_length_old[i] - v_c * c_c) * dx / tau
                    + 0.25 * (-c_r * v_c * v_r - c_c * v_c * v_r)
                    + 0.25
                        * (-c_r * v_c * v_c - c_c * v_c * v_c + c_c * v_l * v_c + c_l * v_l * v_c)
                    + 0.25 * (c_l * v_l * v_l + c_c * v_l * v_l)
                    + 0.25
                        * (c_l * p_l + c_c * p_l - c_l * p_c + c_r * p_c - c_c * p_r - c_r * p_r);

                /* Continuity */
                res[i + n + 1] = (cross_section_length_old[i] - c_c) * dx / tau
                    + 0.25
                        * (c_l * v_l + c_c * v_l + c_l * v_c - c_r * v_c - c_c * v_r - c_r * v_r)
                    + ALPHA * (p_l - 2.0 * p_c + p_r);
            }
            /* Boundary */

            /* Velocity Inlet is prescribed */
            let u0 = 10.0;
            let ampl = 3.0;
            let frequency = 10.0;
            let t_shift = 0.0;
            let velocity_in =
                u0 + ampl * (frequency * (t + tau + t_shift) * std::f64::consts::PI).sin();
            res[0] = velocity_in - velocity[0];

            /* Pressure Inlet is linearly interpolated */
            res[n + 1] = -pressure[0] + 2.0 * pressure[1] - pressure[2];

            /* Velocity Outlet is linearly interpolated */
            res[n] = -velocity[n] + 2.0 * velocity[n - 1] - velocity[n - 2];

            /* Pressure Outlet is "non-reflecting" */
            let tmp2 =
                (c_mk2 - pressure_old[n] / 2.0).sqrt() - (velocity[n] - velocity_old[n]) / 4.0;
            res[2 * n + 1] = -pressure[n] + 2.0 * (c_mk2 - tmp2.powi(2));
            res
        };

        // compute norm of residual
        let norm: f64 = {
            let norm1 = res.map(|x: f64| x * x).sum().sqrt();
            let norm2 = pressure
                .iter()
                .chain(velocity.iter())
                .map(|x| x * x)
                .sum::<f64>()
                .sqrt();
            norm1 / norm2
        };

        const TOL: f64 = 1e-10;
        const MAX_ITER: usize = 1000;
        if norm < TOL && k > 1 {
            //println!("Success at k={} res={} tol={}", k, norm, TOL);
            break;
        } else if k > MAX_ITER {
            panic!("Nonlinear Solver break, iterations: {k}, residual norm: {norm}");
        } else {
            //println!("Continuing k={} res={:.10} tol={:.10}", k, norm, TOL);
        }

        /* Initializing the the LHS i.e. Left Hand Side */
        let lhs = {
            let mut lhs = na::DMatrix::from_element(equations, equations, 0.0);
            for i in 1..n {
                let [v_l, v_c, v_r]: [f64; 3] = velocity[i - 1..i + 2].try_into().unwrap();
                let [c_l, c_c, c_r]: [f64; 3] =
                    cross_section_length[i - 1..i + 2].try_into().unwrap();

                // Momentum, Velocity
                lhs[(i, i - 1)] =
                    0.25 * (-2.0 * c_l * v_l - 2.0 * c_c * v_l - c_c * v_c - c_l * v_c);

                lhs[(i, i)] = c_c * dx / tau
                    + 0.25
                        * (c_r * v_r + c_c * v_r + 2.0 * c_r * v_c + 2.0 * c_c * v_c
                            - c_c * v_l
                            - c_l * v_l);
                lhs[(i, i + 1)] = 0.25 * (c_r * v_c + c_c * v_c);

                // Momentum, Pressure
                lhs[(i, n + 1 + i - 1)] = -0.25 * (c_l + c_c);
                lhs[(i, n + 1 + i)] = 0.25 * (c_l - c_r);
                lhs[(i, n + 1 + i + 1)] = 0.25 * (c_c + c_r);

                // Continuity, Velocity
                lhs[(i + n + 1, i - 1)] = -0.25 * (c_l + c_c);
                lhs[(i + n + 1, i)] = -0.25 * (c_l - c_r);
                lhs[(i + n + 1, i + 1)] = 0.25 * (c_c + c_r);

                // Continuity, Pressure
                lhs[(i + n + 1, n + 1 + i - 1)] = -ALPHA;
                lhs[(i + n + 1, n + 1 + i)] = 2.0 * ALPHA;
                lhs[(i + n + 1, n + 1 + i + 1)] = -ALPHA;
            }

            /* Boundary */

            // Velocity Inlet is prescribed
            lhs[(0, 0)] = 1.0;
            // Pressure Inlet is linearly interpolated
            lhs[(n + 1, n + 1)] = 1.0;
            lhs[(n + 1, n + 2)] = -2.0;
            lhs[(n + 1, n + 3)] = 1.0;
            // Velocity Outlet is linearly interpolated
            lhs[(n, n)] = 1.0;
            lhs[(n, n - 1)] = -2.0;
            lhs[(n, n - 2)] = 1.0;
            // Pressure Outlet is non-Reflecting
            lhs[(2 * n + 1, 2 * n + 1)] = 1.0;
            lhs[(2 * n + 1, n)] =
                -((c_mk2 - pressure_old[n] / 2.0).sqrt() - (velocity[n] - velocity_old[n]) / 4.0);
            lhs
        };

        let sol = lhs
            .lu()
            .solve(&res)
            .expect("Linear Solver did not converge!");

        velocity
            .iter_mut()
            .zip(sol.as_slice()[..n + 1].iter())
            .for_each(|p| *p.0 += p.1);
        pressure
            .iter_mut()
            .zip(sol.as_slice()[n + 1..].iter())
            .for_each(|p| *p.0 += p.1);
    }
}

fn main() -> ExitCode {
    println!("Starting Fluid Solver...");

    let args: Vec<_> = env::args().collect();

    if args.len() != 2 {
        println!("Fluid: Usage: {} <configurationFileName>", args[0]);
        return ExitCode::from(1);
    }

    let config = &args[1];

    const DOMAIN_SIZE: usize = 100;
    const CHUNK_SIZE: usize = DOMAIN_SIZE + 1;

    let mut participant = precice::Participant::new("Fluid", &config, 0, 1);

    println!("preCICE configured...");

    let mesh_name = "Fluid-Nodes-Mesh";
    let dimensions = participant.get_mesh_dimensions(mesh_name);
    assert!(participant.get_data_dimensions(mesh_name, "CrossSectionLength") == 1);
    assert!(participant.get_data_dimensions(mesh_name, "Pressure") == 1);

    const KAPPA: f64 = 100_f64;
    const L: f64 = 10_f64;
    const DX: f64 = L / KAPPA;

    const R0: f64 = 2.0 * std::f64::consts::FRAC_2_SQRT_PI;
    const A0: f64 = R0 * R0 * std::f64::consts::PI;
    const U0: f64 = 10.0;
    const AMPL: f64 = 3.0;
    const FREQUENCY: f64 = 10.0;
    const T_SHIFT: f64 = 0.0;
    const P0: f64 = 0.0;
    let vel_in0: f64 = U0 + AMPL * (FREQUENCY * T_SHIFT * std::f64::consts::PI).sin();

    let mut pressure: Vec<f64> = vec![P0; CHUNK_SIZE as usize];
    let mut cross_section_length: Vec<f64> = vec![A0; CHUNK_SIZE as usize];
    let mut velocity: Vec<f64> = vec![vel_in0; CHUNK_SIZE as usize];

    let mut pressure_old = pressure.clone();

    const CELLWIDTH: f64 = L / DOMAIN_SIZE as f64;
    let grid_size = CHUNK_SIZE * dimensions as usize;
    let grid: Vec<f64> = {
        let mut v: Vec<f64> = vec![0_f64; grid_size];
        for i in 0..CHUNK_SIZE {
            let idx = i * dimensions as usize;
            v[idx] = i as f64 * CELLWIDTH as f64;
        }
        v
    };

    let vertex_ids = {
        let mut ids = vec![-1; CHUNK_SIZE];
        participant.set_mesh_vertices(mesh_name, &grid[..], &mut ids[..]);
        ids
    };

    println!("Initializing preCICE...");

    if participant.requires_initial_data() {
        participant.write_data(mesh_name, "Pressure", &vertex_ids[..], &pressure[..]);
    }

    participant.initialize();

    participant.read_data(
        mesh_name,
        "CrossSectionLength",
        &vertex_ids[..],
        0.0,
        &mut cross_section_length[..],
    );

    let mut cross_section_length_old = cross_section_length.clone();

    // initialize such that mass conservation is fulfilled
    let mut velocity_old = {
        let csl0 = cross_section_length[0];
        cross_section_length
            .iter()
            .map(|csl| vel_in0 * csl0 / csl)
            .collect::<Vec<f64>>()
    };

    let mut out_counter = 0;

    let mut t = 0.0;

    while participant.is_coupling_ongoing() {
        if participant.requires_writing_checkpoint() {
            // do nothing
        }

        let dt = participant.get_max_time_step_size();

        fluid_compute_solution(
            &velocity_old[..],
            &pressure_old[..],
            &cross_section_length_old[..],
            &cross_section_length[..],
            t,
            DOMAIN_SIZE,
            KAPPA,
            dt,
            &mut velocity[..],
            &mut pressure[..],
        );

        participant.write_data(mesh_name, "Pressure", &vertex_ids[..], &pressure[..]);

        participant.advance(dt);

        participant.read_data(
            mesh_name,
            "CrossSectionLength",
            &vertex_ids[..],
            participant.get_max_time_step_size(),
            &mut cross_section_length[..],
        );

        if participant.requires_reading_checkpoint() {
            // i.e. fluid not yet converged

            //pressure.copy_from_slice(&pressure_old[..]);
            //velocity.copy_from_slice(&velocity_old[..]);
        } else {
            t += dt;
            pressure_old.copy_from_slice(&pressure[..]);
            velocity_old.copy_from_slice(&velocity[..]);
            cross_section_length_old.copy_from_slice(&cross_section_length[..]);

            let filename = format!("./output/out_fluid{out_counter}.vtk");
            println!("writing timestep at t={} to {}", t, filename);
            write_vtk(
                &filename,
                DX,
                &velocity_old[..],
                &pressure_old[..],
                &cross_section_length_old[..],
            )
            .expect("Unable to write the vtk file");
            out_counter += 1;
        }
    }

    println!("Exiting FluidSolver at t={}", t);
    participant.finalize();

    return ExitCode::SUCCESS;
}
