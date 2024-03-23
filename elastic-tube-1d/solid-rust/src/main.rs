use std::env;
use std::process::ExitCode;

mod solver;

fn main() -> ExitCode {
    println!("Starting Solid Solver...");

    let args: Vec<_> = env::args().collect();

    if args.len() != 2 {
        println!("Fluid: Usage: {} <configurationFileName>", args[0]);
        return ExitCode::from(1);
    }

    let config = &args[1];

    const DOMAIN_SIZE: usize = 100;
    const CHUNK_SIZE: usize = DOMAIN_SIZE + 1;
    const TUBE_LENGTH: f64 = 10.0;

    let mut participant = precice::Participant::new("Solid", config, 0, 1);

    println!("preCICE configured...");

    let mesh_name = "Solid-Nodes-Mesh";
    let dimensions = participant.get_mesh_dimensions(mesh_name);
    assert!(participant.get_data_dimensions(mesh_name, "CrossSectionLength") == 1);
    assert!(participant.get_data_dimensions(mesh_name, "Pressure") == 1);

    let mut pressure: Vec<f64> = vec![0.0; CHUNK_SIZE];
    let mut cross_section_length: Vec<f64> = vec![1.0; CHUNK_SIZE];

    let grid_size = CHUNK_SIZE * dimensions as usize;
    let grid: Vec<f64> = {
        let mut v: Vec<f64> = vec![0_f64; grid_size];
        const DX: f64 = TUBE_LENGTH / DOMAIN_SIZE as f64;
        for i in 0..CHUNK_SIZE - 1 {
            v[i * dimensions as usize] = DX * i as f64;
        }
        v
    };

    let vertex_ids = {
        let mut ids = vec![-1; CHUNK_SIZE];
        participant.set_mesh_vertices(mesh_name, &grid[..], &mut ids[..]);
        ids
    };

    if participant.requires_initial_data() {
        participant.write_data(
            mesh_name,
            "CrossSectionLength",
            &vertex_ids[..],
            &cross_section_length[..],
        );
    }

    println!("Initializing preCICE...");

    participant.initialize();

    let mut t = 0.0;
    while participant.is_coupling_ongoing() {
        if participant.requires_writing_checkpoint() {
            // no nothing
        }

        let dt = participant.get_max_time_step_size();

        participant.read_data(
            mesh_name,
            "Pressure",
            &vertex_ids[..],
            dt,
            &mut pressure[..],
        );

        solver::solid_compute_solution(&pressure, &mut cross_section_length);

        participant.write_data(
            mesh_name,
            "CrossSectionLength",
            &vertex_ids[..],
            &cross_section_length[..],
        );

        participant.advance(dt);

        if participant.requires_reading_checkpoint() {
            // i.e. fluid not yet converged
            // do nothing
        } else {
            t += dt;
        }
    }

    println!("Exiting SolidSolver at t={}", t);
    participant.finalize();

    ExitCode::SUCCESS
}
