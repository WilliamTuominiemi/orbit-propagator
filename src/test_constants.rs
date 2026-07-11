use crate::constants;
use crate::types;

// Values are NORAD SPACETRACK REPORT NO. 3 SGP4 sample test case input parameters
pub const XNO: f64 = 16.05824518 * (constants::TWOPI / constants::XMNPDA);
pub const XINCL: f64 = 72.8435 * constants::DE2RA;
pub const EO: f64 = 0.0086731;
pub const BSTAR: f64 = 0.000066816;
pub const OMEGAO: f64 = 52.6988 * constants::DE2RA;
pub const XMO: f64 = 110.5714 * constants::DE2RA;
pub const TSINCE: f64 = 0.0;
pub const XNODEO: f64 = 115.9689 * constants::DE2RA;
pub const E6A: f64 = 0.000001;

// Values from SGP4 sample test case output values
pub const POSITION_AND_VELOCITY_0: types::PositionAndVelocity = types::PositionAndVelocity {
    x: 2328.97048951,
    y: -5995.22076416,
    z: 1719.97067261,
    xdot: 2.91207230,
    ydot: -0.98341546,
    zdot: -7.09081703,
};

pub const TOLERANCE: f64 = 1e-3;
