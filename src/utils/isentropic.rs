#![allow(dead_code)]

use std::f64::consts::PI;
use super::numerics::*;

pub fn calc_mach_angle_from_mach(mach_number: f64) -> Result<f64, &'static str> {
    if mach_number < 0.0 {
        return Err("invalid mach number");
    }
    let mach_angle = (1.0 / mach_number).asin();
    Ok(mach_angle)
}

pub fn calc_pressure_ratio_from_mach(mach_number: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err("invalid specific heat ratio");
    }
    let pressure_ratio: f64 = (1.0 + (specific_heat_ratio - 1.0) / 2.0 * mach_number.powi(2)).powf(-specific_heat_ratio / (specific_heat_ratio - 1.0));
    Ok(pressure_ratio)
}

pub fn calc_temperature_ratio_from_mach(mach_number: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err("invalid specific heat ratio");
    }
    let temperature_ratio: f64 = (1.0 + (specific_heat_ratio - 1.0) / 2.0 * mach_number.powi(2)).powi(-1);
    Ok(temperature_ratio)
}

pub fn calc_density_ratio_from_mach(mach_number: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err("invalid specific heat ratio");
    }
    let density_ratio: f64 = (1.0 + (specific_heat_ratio - 1.0) / 2.0 * mach_number.powi(2)).powf(-1.0 / (specific_heat_ratio - 1.0));
    Ok(density_ratio)
}

pub fn calc_mach_from_speed_of_sound(velocity: f64, speed_of_sound: f64) -> Result<f64, &'static str> {
    Ok(velocity / speed_of_sound)
}

pub fn prandtl_meyer_function(mach_number: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err("invalid specific heat ratio");
    }
    if mach_number <= 1.0 {
        return Err("invalid mach number");
    }
    let gamma_ratio = (specific_heat_ratio - 1.0) / (specific_heat_ratio + 1.0);
    let sqrt_gamma_ratio = gamma_ratio.sqrt();
    let prandtl_meyer_angle: f64 = (
        (1.0 / sqrt_gamma_ratio) // 1st term
        * (sqrt_gamma_ratio * (mach_number.powi(2) - 1.0).sqrt()).atan()) // 2nd term
        - (mach_number.powi(2) - 1.0).sqrt().atan(); // 3rd term
    Ok(prandtl_meyer_angle)
}

pub fn calc_mach_from_mach_angle(mach_angle: f64) -> Result<f64, &'static str> {
    if mach_angle < 0.0 || mach_angle > PI / 2.0 {
        // check valid mach angle in radians
        return Err("invalid mach angle")
    }
    let mach_number: f64 = 1.0 / mach_angle.sin();
    Ok(mach_number)
}

pub fn calc_mach_from_temperature_ratio(temperature_ratio: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err("invalid specific heat ratio");
    }
    if temperature_ratio <= 0.0 || temperature_ratio > 1.0 {
        // check valid temperature ratio
        return Err("invalid temperature ratio");
    }
    let mach_number: f64 = (2.0 * ((1.0 / temperature_ratio) - 1.0) / (specific_heat_ratio - 1.0)).sqrt();
    Ok(mach_number)
}

pub fn calc_mach_from_pressure_ratio(pressure_ratio: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err("invalid specific heat ratio");
    }
    if pressure_ratio <= 0.0 || pressure_ratio > 1.0 {
        // check valid pressure ratio
        return Err("invalid pressure ratio");
    }
    let mach_number: f64 = (2.0 * ((1.0 / pressure_ratio.powf((specific_heat_ratio - 1.0) / specific_heat_ratio)) - 1.0) / (specific_heat_ratio - 1.0)).sqrt();
    Ok(mach_number)
}

pub fn calc_mach_from_density_ratio(density_ratio: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err("invalid specific heat ratio");
    }
    if density_ratio <= 0.0 || density_ratio > 1.0 {
        // check valid density ratio
        return Err("invalid density ratio");
    }
    let mach_number: f64 = ((2.0 * ((1.0 / density_ratio.powf(specific_heat_ratio - 1.0)) - 1.0)) / (specific_heat_ratio - 1.0)).sqrt();
    Ok(mach_number)
}

pub fn calc_mach_from_prandtl_meyer_angle(prandtl_meyer_angle: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    // eta is the value of (m^2 - 1).sqrt()
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err("invalid specific heat ratio");
    }
    let alpha = ((specific_heat_ratio + 1.0) / (specific_heat_ratio - 1.0)).sqrt(); // just a constant to make things easier
    let f = |eta: f64| {
        alpha * (eta / alpha).atan()
        - eta.atan() - prandtl_meyer_angle    
    };
    let df = |eta: f64| {
        1.0 / ((eta / alpha).powi(2) + 1.0)
        - 1.0 / (eta.powi(2) + 1.0)
    };
    let eta: f64 = newton_raphson(&f, &df, 1.5, None, None);
    let mach_number: f64 = (eta.powi(2) + 1.0).sqrt();
    Ok(mach_number)
}

pub fn valid_specific_heat_ratio(specific_heat_ratio: f64) -> bool {
    // specific heat ratio must be greater than 1
    specific_heat_ratio > 1.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_specific_heat_ratio() {
        // check that valid specific heat ratio returns true for > 1
        assert!(valid_specific_heat_ratio(1.4));
        assert!(!valid_specific_heat_ratio(1.0));
    }

    #[test]
    fn test_calc_mach_angle_from_mach() {
        // test mach angle calculation for a valid mach number
        let mach_number = 2.0;
        let result = calc_mach_angle_from_mach(mach_number).expect("valid mach number");
        let expected = (1.0 / mach_number).asin();
        assert!((result - expected).abs() < 1e-6);
    }

    #[test]
    fn test_calc_pressure_ratio_from_mach() {
        // test pressure ratio calculation
        let mach_number = 2.0;
        let specific_heat_ratio = 1.4;
        let result = calc_pressure_ratio_from_mach(mach_number, specific_heat_ratio)
            .expect("valid pressure ratio");
        // expected = (1 + 0.2 * mach_number^2)^(-specific_heat_ratio/(specific_heat_ratio-1))
        let base = 1.0 + (specific_heat_ratio - 1.0) / 2.0 * mach_number.powi(2);
        let expected = base.powf(-specific_heat_ratio / (specific_heat_ratio - 1.0));
        assert!((result - expected).abs() < 1e-6);
    }

    #[test]
    fn test_calc_temperature_ratio_from_mach() {
        // test temperature ratio calculation
        let mach_number = 2.0;
        let specific_heat_ratio = 1.4;
        let result = calc_temperature_ratio_from_mach(mach_number, specific_heat_ratio)
            .expect("valid temperature ratio");
        let base = 1.0 + (specific_heat_ratio - 1.0) / 2.0 * mach_number.powi(2);
        let expected = base.powi(-1);
        assert!((result - expected).abs() < 1e-6);
    }

    #[test]
    fn test_calc_density_ratio_from_mach() {
        // test density ratio calculation
        let mach_number = 2.0;
        let specific_heat_ratio = 1.4;
        let result = calc_density_ratio_from_mach(mach_number, specific_heat_ratio)
            .expect("valid density ratio");
        let base = 1.0 + (specific_heat_ratio - 1.0) / 2.0 * mach_number.powi(2);
        let expected = base.powf(-1.0 / (specific_heat_ratio - 1.0));
        assert!((result - expected).abs() < 1e-6);
    }

    #[test]
    fn test_calc_mach_from_speed_of_sound() {
        // test mach calculation from speed of sound
        let velocity = 680.0;
        let speed_of_sound = 340.0;
        let result = calc_mach_from_speed_of_sound(velocity, speed_of_sound)
            .expect("calculation should succeed");
        let expected = velocity / speed_of_sound;
        assert!((result - expected).abs() < 1e-6);
    }

    #[test]
    fn test_prandtl_meyer_function() {
        // test prandtl meyer function for valid mach number > 1
        let mach_number = 2.0;
        let specific_heat_ratio = 1.4;
        let result = prandtl_meyer_function(mach_number, specific_heat_ratio)
            .expect("valid prandtl meyer calculation");
        // expected value calculated externally (approximate)
        // note: due to iterative nature, a small tolerance is allowed
        let expected = {
            let gamma_ratio = (specific_heat_ratio - 1.0) / (specific_heat_ratio + 1.0);
            let sqrt_gamma_ratio = gamma_ratio.sqrt();
            let term1 = 1.0 / sqrt_gamma_ratio * ((sqrt_gamma_ratio * (mach_number.powi(2) - 1.0).sqrt()).atan());
            let term2 = (mach_number.powi(2) - 1.0).sqrt().atan();
            term1 - term2
        };
        assert!((result - expected).abs() < 1e-6);
    }

    #[test]
    fn test_calc_mach_from_mach_angle() {
        // test mach number calculation from a valid mach angle
        let mach_angle = (1.0f64 / 2.0f64).asin();
        let result = calc_mach_from_mach_angle(mach_angle).expect("valid mach angle");
        let expected = 2.0;
        assert!((result - expected).abs() < 1e-6);
    }

    #[test]
    fn test_calc_mach_from_temperature_ratio() {
        // test mach calculation from temperature ratio
        let specific_heat_ratio = 1.4;
        // for mach=2.0, temperature_ratio = (1 + (0.4/2 * 4))^-1 = 1/1.8
        let temperature_ratio = 1.0 / 1.8;
        let result = calc_mach_from_temperature_ratio(temperature_ratio, specific_heat_ratio)
            .expect("valid temperature ratio");
        let expected = 2.0;
        assert!((result - expected).abs() < 1e-6);
    }

    #[test]
    fn test_calc_mach_from_pressure_ratio() {
        // test mach calculation from pressure ratio
        let specific_heat_ratio = 1.4;
        // for mach=2.0, pressure_ratio = 1.8^(-3.5)
        let base = 1.0 + (specific_heat_ratio - 1.0) / 2.0 * 2.0f64.powi(2);
        let pressure_ratio = base.powf(-specific_heat_ratio / (specific_heat_ratio - 1.0));
        let result = calc_mach_from_pressure_ratio(pressure_ratio, specific_heat_ratio)
            .expect("valid pressure ratio");
        let expected = 2.0;
        assert!((result - expected).abs() < 1e-5);
    }

    #[test]
    fn test_calc_mach_from_density_ratio() {
        // test mach calculation from density ratio
        let specific_heat_ratio = 1.4;
        // for mach=2.0, density_ratio = 1.8^(-2.5)
        let base = 1.0 + (specific_heat_ratio - 1.0) / 2.0 * 2.0f64.powi(2);
        let density_ratio = base.powf(-1.0 / (specific_heat_ratio - 1.0));
        let result = calc_mach_from_density_ratio(density_ratio, specific_heat_ratio)
            .expect("valid density ratio");
        let expected = 2.0;
        assert!((result - expected).abs() < 1e-5);
    }

    #[test]
    fn test_calc_mach_from_prandtl_meyer_angle() {
        // test mach calculation from prandtl meyer angle
        let specific_heat_ratio = 1.4;
        // calculate prandtl meyer angle for mach=2.0 as input
        let prandtl_angle = prandtl_meyer_function(2.0, specific_heat_ratio)
            .expect("valid prandtl meyer function");
        let result = calc_mach_from_prandtl_meyer_angle(prandtl_angle, specific_heat_ratio)
            .expect("valid mach from prandtl meyer angle calculation");
        let expected = 2.0;
        assert!((result - expected).abs() < 1e-5);
    }

    #[test]
    fn test_invalid_values() {
        // test error conditions for invalid inputs
        // invalid mach number for calc_mach_angle_from_mach (negative)
        assert!(calc_mach_angle_from_mach(-1.0).is_err());
        // invalid specific heat ratio for calc_pressure_ratio_from_mach
        assert!(calc_pressure_ratio_from_mach(2.0, 1.0).is_err());
        // invalid mach number for prandtl_meyer_function (<=1)
        assert!(prandtl_meyer_function(1.0, 1.4).is_err());
        // invalid mach angle for calc_mach_from_mach_angle (out of range)
        assert!(calc_mach_from_mach_angle(PI).is_err());
    }
}