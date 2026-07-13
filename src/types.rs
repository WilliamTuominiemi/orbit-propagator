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
    pub temp2: f64,
    pub temp3: f64,
    pub temp4: f64,
    pub temp5: f64,
    pub temp6: f64,
    pub sinepw: f64,
    pub cosepw: f64,
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
}

pub struct GeodeticPosition {
    pub lat: f64,
    pub lon: f64,
    pub alt: f64,
}
