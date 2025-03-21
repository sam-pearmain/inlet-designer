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