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
pub const MID_TOLERANCE: f64 = 1e-9;
pub const SMALL_TOLERANCE: f64 = 1e-12;

pub const SGP4TLE: &str = "SGP4 (SGP4)
1 88888U 80081S 80275.98708465 .00073094 13844-3 66816-4 0 6
2 88888 72.8435 115.9689 0086731 52.6988 110.5714 16.05824718105";
