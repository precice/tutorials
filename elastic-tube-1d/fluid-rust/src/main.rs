use std::env;
use std::process::ExitCode;

mod solver;
mod utils;

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

    let mut participant = precice::Participant::new("Fluid", config, 0, 1);

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

    let mut pressure: Vec<f64> = vec![P0; CHUNK_SIZE];
    let mut cross_section_length: Vec<f64> = vec![A0; CHUNK_SIZE];
    let mut velocity: Vec<f64> = vec![vel_in0; CHUNK_SIZE];

    let mut pressure_old = pressure.clone();

    const CELLWIDTH: f64 = L / DOMAIN_SIZE as f64;
    let grid_size = CHUNK_SIZE * dimensions as usize;
    let grid: Vec<f64> = {
        let mut v: Vec<f64> = vec![0_f64; grid_size];
        for i in 0..CHUNK_SIZE {
            let idx = i * dimensions as usize;
            v[idx] = i as f64 * CELLWIDTH;
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

        solver::fluid_compute_solution(
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

            let filename = format!("./output/out_fluid_{out_counter}.vtk");
            println!("writing timestep at t={} to {}", t, filename);
            utils::write_vtk(
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

    ExitCode::SUCCESS
}
