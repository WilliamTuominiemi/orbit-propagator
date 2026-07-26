pub struct MeanMotionAndSemimajorAxisOutput {
    pub xnodp: f64,
    pub aodp: f64,
    pub betao2: f64,
    pub betao: f64,
    pub x3thm1: f64,
    pub theta2: f64,
    pub cosio: f64,
}

#[derive(Debug)]
pub struct SecularGravityAndAtmosphericDragUpdateOutput {
    pub e: f64,
    pub a: f64,
    pub xl: f64,
    pub beta: f64,
    pub xn: f64,
    pub xnode: f64,
    pub omega: f64,
}

pub struct KeplersEquationOutput {
    pub sinepw: f64,
    pub cosepw: f64,
    pub ecose: f64,
    pub esine: f64,
}

#[derive(Debug, PartialEq)]
pub struct ShortPeriodicsOutput {
    pub rk: f64,
    pub uk: f64,
    pub xnodek: f64,
    pub xinck: f64,
    pub rdotk: f64,
    pub rfdotk: f64,
}

#[derive(Debug, PartialEq)]
pub struct LongPeriodicsOutput {
    pub xlt: f64,
    pub ayn: f64,
    pub axn: f64,
}

#[derive(Debug, PartialEq)]
pub struct ShortPeriodPrelimenaryQuantities {
    pub r: f64,
    pub rdot: f64,
    pub rfdot: f64,
    pub temp2: f64,
    pub betal: f64,
    pub temp1: f64,
    pub cos2u: f64,
    pub u: f64,
    pub sin2u: f64,
}

#[derive(Debug, PartialEq)]
pub struct CConstants {
    pub c1: f64,
    pub c2: f64,
    pub c3: f64,
    pub c4: f64,
    pub c5: f64,
}
#[derive(Debug, PartialEq)]
pub struct DConstants {
    pub d2: f64,
    pub d3: f64,
    pub d4: f64,
}

#[derive(Debug, PartialEq)]
pub struct OrientationVectors {
    pub ux: f64,
    pub uy: f64,
    pub uz: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
}

#[derive(Debug, PartialEq)]
pub struct PositionAndVelocity {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub xdot: f64,
    pub ydot: f64,
    pub zdot: f64,
}

#[derive(Debug, PartialEq)]
pub struct RotationMatrix {
    pub m0: f64,
    pub m1: f64,
    pub m2: f64,
    pub m3: f64,
    pub m4: f64,
    pub m5: f64,
    pub m6: f64,
    pub m7: f64,
    pub m8: f64,
}

#[derive(Debug)]
pub struct EcefPosition {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub xdot: f64,
    pub ydot: f64,
    pub zdot: f64,
}

#[derive(Debug)]
pub struct GeodeticPosition {
    pub lat: f64,
    pub lon: f64,
    pub alt: f64,
    pub vel_east: f64,
    pub vel_north: f64,
    pub vel_up: f64,
}

pub struct GraphDataPoint {
    pub point: [f64; 2],
    pub altitude: f64,
    pub velocity: f64,
}

pub struct TLE {
    pub name: String,
    pub number: u32,
    pub international_designator: String,
    pub epoch_year_julian_fraction: f64,
    pub first_derivative_of_mean_motion: f64,
    pub second_derivative_of_mean_motion: f64,
    pub drag_term: f64,
    pub ephemeris_type: u32,
    pub element_number_check_sum: u32,
    pub inclination: f64,
    pub right_ascension_of_ascending_node: f64,
    pub eccentricity: f64,
    pub argument_of_perigee: f64,
    pub mean_anomaly: f64,
    pub mean_motion: f64,
    pub revolution_number_check_sum: u32,
}
