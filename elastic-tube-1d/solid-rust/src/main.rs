use std::env;
use std::process::ExitCode;

fn solid_compute_solution(pressure: &[f64], cross_section_length: &mut [f64]) {
    assert!(pressure.len() == cross_section_length.len());
    const E: f64 = 10000.0;
    const C_MK2: f64 = E / std::f64::consts::FRAC_2_SQRT_PI;
    const PRESSURE0: f64 = 0.0;
    let new: Vec<_> = pressure
        .iter()
        .map(|p| ((PRESSURE0 - 2.0 * C_MK2) / (p - 2.0 * C_MK2)).powi(2))
        .collect();
    cross_section_length.copy_from_slice(&new[..]);
}

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

    let mut participant = precice::new("Solid", &config, 0, 1);

    println!("preCICE configured...");

    let mesh_name = "Solid-Nodes-Mesh";
    let dimensions = participant.get_mesh_dimensions(mesh_name);
    assert!(participant.get_data_dimensions(mesh_name, "CrossSectionLength") == 1);
    assert!(participant.get_data_dimensions(mesh_name, "Pressure") == 1);

    let mut pressure: Vec<f64> = vec![0.0; CHUNK_SIZE as usize];
    let mut cross_section_length: Vec<f64> = vec![1.0; CHUNK_SIZE as usize];

    let grid_size = CHUNK_SIZE * dimensions as usize;
    let grid: Vec<f64> = {
        let mut v: Vec<f64> = vec![0_f64; grid_size];
        for i in 0..CHUNK_SIZE - 1 {
            for j in 0..(dimensions as usize) - 1 {
                let idx = i * dimensions as usize + j;
                v[idx] = (i * (1 - j)) as f64;
            }
        }
        v
    };

    let vertex_ids = {
        let mut ids = vec![-1; CHUNK_SIZE];
        participant
            .pin_mut()
            .set_mesh_vertices(mesh_name, &grid[..], &mut ids[..]);
        ids
    };

    if participant.pin_mut().requires_initial_data() {
        participant.pin_mut().write_data(
            mesh_name,
            "CrossSectionLength",
            &vertex_ids[..],
            &cross_section_length[..],
        );
    }


    println!("Initializing preCICE...");

    participant.pin_mut().initialize();

    let mut t = 0.0;
    while participant.is_coupling_ongoing() {
        if participant.pin_mut().requires_writing_checkpoint() {
            // no nothing
        }

        let dt = participant.get_max_time_step_size();

        participant.read_data(mesh_name, "Pressure", &vertex_ids[..], dt, &mut pressure[..]);

        solid_compute_solution(&pressure, &mut cross_section_length);

        participant.pin_mut().write_data(
            mesh_name,
            "CrossSectionLength",
            &vertex_ids[..],
            &cross_section_length[..],
        );

        participant.pin_mut().advance(dt);

        if participant.pin_mut().requires_reading_checkpoint() {
            // i.e. fluid not yet converged
            // do nothing
        } else {
            t += dt;
        }
    }

    println!("Exiting SolidSolver at t={}", t);
    participant.pin_mut().finalize();

    return ExitCode::SUCCESS;
}
