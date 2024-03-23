use nalgebra::{DMatrix, DVector};

pub fn fluid_compute_solution(
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
    pressure.copy_from_slice(pressure_old);

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
            let mut res = DVector::from_element(equations, 0.0);
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
            let mut lhs = DMatrix::from_element(equations, equations, 0.0);
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
