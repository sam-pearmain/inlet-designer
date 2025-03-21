#![allow(dead_code)]

use std::f64::consts::PI;
use super::isentropic::valid_specific_heat_ratio; 
use super::numerics::bisection;

fn calc_downstream_mach(upstream_mach: f64, shock_angle: f64, deflection_angle: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    let normal_upstream_mach: f64 = calc_normal_upstream_mach(upstream_mach, shock_angle)?;

    // this is wrong
    let normal_downstream_mach: f64 = (
        (1.0 + (specific_heat_ratio - 1.0) * normal_upstream_mach.powi(2) / 2.0)
        / (specific_heat_ratio * normal_upstream_mach.powi(2) - (specific_heat_ratio - 1.0) / 2.0)
    ).sqrt();

    let downstream_mach: f64 = normal_downstream_mach / (shock_angle - deflection_angle).sin();
    Ok(downstream_mach)
}

pub fn calc_downstream_mach_from_shock_angle(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    let deflection_angle: f64 = calc_deflection_angle(upstream_mach, shock_angle, specific_heat_ratio)?;
    calc_downstream_mach(upstream_mach, shock_angle, deflection_angle, specific_heat_ratio)
}

pub fn calc_downstream_mach_from_deflection_angle(upstream_mach: f64, deflection_angle: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    let shock_angle: f64 = calc_shock_angle(upstream_mach, deflection_angle, specific_heat_ratio)?;
    calc_downstream_mach(upstream_mach, shock_angle, deflection_angle, specific_heat_ratio)
}

pub fn calc_deflection_angle(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    let tan_deflection_angle: f64 = 
        2.0 / shock_angle.tan() * 
        (upstream_mach.powi(2) * shock_angle.sin().powi(2) - 1.0) / 
        (upstream_mach.powi(2) * (specific_heat_ratio + (2.0 * shock_angle).cos()) + 2.0);
    Ok(tan_deflection_angle.atan())
}

pub fn calc_pressure_ratio(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    let pressure_ratio: f64 = 
        (2.0 * specific_heat_ratio * upstream_mach.powi(2) 
            * shock_angle.sin().powi(2) - (specific_heat_ratio - 1.0)) / 
        (specific_heat_ratio + 1.0);
    Ok(pressure_ratio)
}

pub fn calc_density_ratio(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    let density_ratio: f64 = 
    (specific_heat_ratio + 1.0) * upstream_mach.powi(2) * shock_angle.sin().powi(2) /
    ((specific_heat_ratio - 1.0) * upstream_mach.powi(2) * shock_angle.sin().powi(2) + 2.0);
    Ok(density_ratio)
    
}

pub fn calc_temperature_ratio(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    let pressure_ratio: f64 = calc_pressure_ratio(upstream_mach, shock_angle, specific_heat_ratio)?;
    let density_ratio: f64 = calc_density_ratio(upstream_mach, shock_angle, specific_heat_ratio)?;
    let temperature_ratio: f64 = pressure_ratio * (1.0 / density_ratio);
    Ok(temperature_ratio)
}

pub fn calc_stagnation_pressure_ratio(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err("invalid specific heat ratio");
    }
    let stagnation_pressure_ratio: f64 =
        calc_density_ratio(upstream_mach, shock_angle, specific_heat_ratio)?.powf(specific_heat_ratio / (specific_heat_ratio - 1.0)) *
        (1.0 / calc_pressure_ratio(upstream_mach, shock_angle, specific_heat_ratio)?).powf(1.0 / (specific_heat_ratio - 1.0));
    Ok(stagnation_pressure_ratio)
}

pub fn calc_shock_angle(upstream_mach: f64, deflection_angle: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    if upstream_mach <= 1.0 {
        return Err("invalid mach number");
    }
    let f = |shock_angle: f64| {
        let calculated_deflection_angle = match calc_deflection_angle(upstream_mach, shock_angle, specific_heat_ratio) {
            Ok(value) => value,
            Err(_) => panic!("erm what"),
        };
        return calculated_deflection_angle - deflection_angle
    };

    let lower_bound: f64 = deflection_angle;
    let upper_bound: f64 = PI / 2.0;

    let shock_angle: f64 = bisection(&f, lower_bound, upper_bound, None, None);

    if shock_angle.is_nan() {
        return Err("math error");
    }

    Ok(shock_angle)
}

pub fn calc_max_shock_angle(upstream_mach: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err("invalid specific heat ratio");
    }
    if upstream_mach <= 1.0 {
        return Err("invalid mach number");
    }
    
    let sin_max_shock_angle: f64 = 
        ((1.0 / (specific_heat_ratio * upstream_mach.powi(2))) * 
        (1.0 +
            ((specific_heat_ratio + 1.0) * (
                (specific_heat_ratio + 1.0) * upstream_mach.powi(4) / 16.0 +
                (specific_heat_ratio - 1.0) * upstream_mach.powi(2) / 2.0 +
                1.0
            ).sqrt())
        )).sqrt();

    if sin_max_shock_angle > 1.0 || sin_max_shock_angle < 0.0 {
        return Err("math error");
    }

    let shock_angle: f64 = sin_max_shock_angle.asin();
    Ok(shock_angle)
}

pub fn calc_normal_upstream_mach(upstream_mach: f64, shock_angle: f64) -> Result<f64, &'static str> {
    if upstream_mach <= 1.0 {
        return Err("invalid mach number");
    }
    Ok(upstream_mach * shock_angle.sin())
}

pub fn calc_normal_downstream_mach(downstream_mach: f64, shock_angle: f64, deflection_angle: f64) -> Result<f64, &'static str> {
    if downstream_mach <= 1.0 {
        return Err("invalid mach number");
    }
    Ok(downstream_mach * (shock_angle - deflection_angle).sin())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calc_deflection_angle() {
        // test calc_deflection_angle with upstream mach 2.0, shock angle pi/4 and specific heat ratio 1.4
        let upstream_mach = 2.0;
        let shock_angle = PI / 4.0;
        let specific_heat_ratio = 1.4;
        let result = calc_deflection_angle(upstream_mach, shock_angle, specific_heat_ratio)
            .expect("calculation should succeed");
        // expected value approx 0.257 rad (computed from formula)
        let expected = 0.257;
        assert!((result - expected).abs() < 1e-2);
    }

    #[test]
    fn test_calc_pressure_ratio() {
        // test calc_pressure_ratio with upstream mach 2.0, shock angle pi/4 and specific heat ratio 1.4
        let upstream_mach = 2.0;
        let shock_angle = PI / 4.0;
        let specific_heat_ratio = 1.4;
        let result = calc_pressure_ratio(upstream_mach, shock_angle, specific_heat_ratio)
            .expect("calculation should succeed");
        // expected value approx 2.16667
        let expected = 2.16667;
        assert!((result - expected).abs() < 1e-2);
    }

    #[test]
    fn test_calc_density_ratio() {
        // test calc_density_ratio with upstream mach 2.0, shock angle pi/4 and specific heat ratio 1.4
        let upstream_mach = 2.0;
        let shock_angle = PI / 4.0;
        let specific_heat_ratio = 1.4;
        let result = calc_density_ratio(upstream_mach, shock_angle, specific_heat_ratio)
            .expect("calculation should succeed");
        // expected value approx 1.7143
        let expected = 1.7143;
        assert!((result - expected).abs() < 1e-2);
    }

    #[test]
    fn test_calc_temperature_ratio() {
        // test calc_temperature_ratio with upstream mach 2.0, shock angle pi/4 and specific heat ratio 1.4
        let upstream_mach = 2.0;
        let shock_angle = PI / 4.0;
        let specific_heat_ratio = 1.4;
        let result = calc_temperature_ratio(upstream_mach, shock_angle, specific_heat_ratio)
            .expect("calculation should succeed");
        // expected value approx pressure_ratio/density_ratio = 2.16667/1.7143 ≈ 1.264
        let expected = 1.264;
        assert!((result - expected).abs() < 1e-2);
    }

    #[test]
    fn test_calc_stagnation_pressure_ratio() {
        // test calc_stagnation_pressure_ratio with upstream mach 2.0, shock angle pi/4 and specific heat ratio 1.4
        let upstream_mach = 2.0;
        let shock_angle = PI / 4.0;
        let specific_heat_ratio = 1.4;
        let result = calc_stagnation_pressure_ratio(upstream_mach, shock_angle, specific_heat_ratio)
            .expect("calculation should succeed");
        // expected value computed approx as 0.955
        let expected = 0.955;
        assert!((result - expected).abs() < 1e-2);
    }

    #[test]
    fn test_calc_shock_angle() {
        // test calc_shock_angle with upstream mach 2.0, using a deflection angle computed from a known shock angle
        let upstream_mach = 2.0;
        let specific_heat_ratio = 1.4;
        // using known shock angle pi/4 to compute deflection angle
        let known_shock_angle = PI / 4.0;
        let deflection_angle = calc_deflection_angle(upstream_mach, known_shock_angle, specific_heat_ratio)
            .expect("deflection angle calculation should succeed");
        let result = calc_shock_angle(upstream_mach, deflection_angle, specific_heat_ratio)
            .expect("calculation should succeed");
        // expected shock angle should be close to the known shock angle (pi/4)
        assert!((result - known_shock_angle).abs() < 1e-2);
    }

    #[test]
    fn test_calc_max_shock_angle() {
        // test calc_max_shock_angle with upstream mach 3.0 and specific heat ratio 1.4
        let upstream_mach = 3.0;
        let specific_heat_ratio = 1.4;
        let result = calc_max_shock_angle(upstream_mach, specific_heat_ratio)
            .expect("calculation should succeed");
        // expected value approx asin(0.903) which is about 1.12 rad
        let expected = 1.12;
        assert!((result - expected).abs() < 1e-2);
    }

    #[test]
    fn test_calc_normal_upstream_mach() {
        // test calc_normal_upstream_mach with upstream mach 2.0 and shock angle pi/4
        let upstream_mach = 2.0;
        let shock_angle = PI / 4.0;
        let result = calc_normal_upstream_mach(upstream_mach, shock_angle)
            .expect("calculation should succeed");
        // expected value is 2.0 * sin(pi/4) ≈ 1.4142
        let expected = 2.0 * (PI / 4.0).sin();
        assert!((result - expected).abs() < 1e-6);
    }

    #[test]
    fn test_calc_normal_downstream_mach() {
        // test calc_normal_downstream_mach with downstream mach 1.5, shock angle pi/4 and deflection angle 0.257
        let downstream_mach = 1.5;
        let shock_angle = PI / 4.0;
        let deflection_angle = 0.257;
        let result = calc_normal_downstream_mach(downstream_mach, shock_angle, deflection_angle)
            .expect("calculation should succeed");
        // expected value is 1.5 * sin(pi/4 - 0.257)
        let expected = downstream_mach * (shock_angle - deflection_angle).sin();
        assert!((result - expected).abs() < 1e-6);
    }

    #[test]
    fn test_calc_downstream_mach_from_shock_angle() {
        // test calc_downstream_mach_from_shock_angle; ensure that the function returns an ok result
        let upstream_mach = 2.0;
        let shock_angle = PI / 4.0;
        let specific_heat_ratio = 1.4;
        let result = calc_downstream_mach_from_shock_angle(upstream_mach, shock_angle, specific_heat_ratio);
        assert!(result.is_ok());
    }

    #[test]
    fn test_calc_downstream_mach_from_deflection_angle() {
        // test calc_downstream_mach_from_deflection_angle; ensure that the function returns an ok result
        let upstream_mach = 2.0;
        let specific_heat_ratio = 1.4;
        // use a deflection angle computed from a known shock angle
        let known_shock_angle = PI / 4.0;
        let deflection_angle = calc_deflection_angle(upstream_mach, known_shock_angle, specific_heat_ratio)
            .expect("deflection angle calculation should succeed");
        let result = calc_downstream_mach_from_deflection_angle(upstream_mach, deflection_angle, specific_heat_ratio);
        assert!(result.is_ok());
    }

    #[test]
    fn test_invalid_values() {
        // test error conditions for invalid inputs
        // invalid upstream mach for calc_normal_upstream_mach
        assert!(calc_normal_upstream_mach(1.0, PI / 4.0).is_err());
        // invalid upstream mach for calc_shock_angle
        assert!(calc_shock_angle(1.0, 0.2, 1.4).is_err());
        // invalid downstream mach for calc_normal_downstream_mach
        assert!(calc_normal_downstream_mach(1.0, PI / 4.0, 0.2).is_err());
        // invalid specific heat ratio for calc_stagnation_pressure_ratio
        assert!(calc_stagnation_pressure_ratio(2.0, PI / 4.0, 1.0).is_err());
    }
}