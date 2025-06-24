pub fn solid_compute_solution(pressure: &[f64], cross_section_length: &mut [f64]) {
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
