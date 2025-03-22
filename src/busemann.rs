#![allow(dead_code)]

use crate::{inlet::Inlet, utils};

pub fn calc_contour_from_machs(freestream_mach: f64, exit_mach: f64) -> Inlet {
    let gamma: f64 = 1.4;


    todo!()
}

pub fn calc_total_pressure_ratio(freestream_mach: f64, exit_mach: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    todo!()
}

pub fn calc_static_temperature_ratio(freestream_mach: f64, exit_mach: f64, specific_heat_ratio: f64) -> Result<f64, &'static str> {
    let freestream_temperature_ratio: f64 = utils::isentropic::calc_temperature_ratio_from_mach(freestream_mach, specific_heat_ratio)?;
    let exit_temperature_ratio: f64 = utils::isentropic::calc_temperature_ratio_from_mach(exit_mach, specific_heat_ratio)?;
    let static_temperature_ratio: f64 = freestream_temperature_ratio * (1.0 / exit_temperature_ratio);
    Ok(static_temperature_ratio)
}