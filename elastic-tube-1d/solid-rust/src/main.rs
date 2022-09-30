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

    let mut interface = precice::new("Solid", &config, 0, 1);

    println!("preCICE configured...");

    let dimensions = interface.get_dimensions();
    let mesh_id = interface.get_mesh_id("Solid-Nodes-Mesh");
    let cross_section_length_id = interface.get_data_id("CrossSectionLength", mesh_id);
    let pressure_id = interface.get_data_id("Pressure", mesh_id);

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
        interface
            .pin_mut()
            .set_mesh_vertices(mesh_id, &grid[..], &mut ids[..]);
        ids
    };

    println!("Initializing preCICE...");

    let mut t = 0.0;
    let dt = interface.pin_mut().initialize();

    if interface.is_action_required(&precice::action_write_initial_data()) {
        interface.pin_mut().write_block_scalar_data(
            cross_section_length_id,
            &vertex_ids[..],
            &cross_section_length[..],
        );
        interface
            .pin_mut()
            .mark_action_fulfilled(&precice::action_write_initial_data());
    }

    interface.pin_mut().initialize_data();

    while interface.is_coupling_ongoing() {
        if interface.is_action_required(&precice::action_write_iteration_checkpoint()) {
            interface
                .pin_mut()
                .mark_action_fulfilled(&precice::action_write_iteration_checkpoint());
        }

        interface.read_block_scalar_data(pressure_id, &vertex_ids[..], &mut pressure[..]);

        solid_compute_solution(&pressure, &mut cross_section_length);

        interface.pin_mut().write_block_scalar_data(
            cross_section_length_id,
            &vertex_ids[..],
            &cross_section_length[..],
        );

        interface.pin_mut().advance(dt);

        if interface.is_action_required(&precice::action_read_iteration_checkpoint()) {
            // i.e. fluid not yet converged
            interface
                .pin_mut()
                .mark_action_fulfilled(&precice::action_read_iteration_checkpoint());
        } else {
            t += dt;
        }
    }

    println!("Exiting SolidSolver at t={}", t);
    interface.pin_mut().finalize();

    return ExitCode::SUCCESS;
}
