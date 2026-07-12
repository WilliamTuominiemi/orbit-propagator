use crate::{
    constants::{self, DE2RA},
    types::{self, PositionAndVelocity, RotationMatrix},
};

pub struct GroundTrack {}

impl GroundTrack {
    pub fn new() -> Self {
        GroundTrack {}
    }

    pub fn eci_to_ecef(&self, ut1: f64, pav: PositionAndVelocity) -> types::EcefPosition {
        let rotation_matrix = self.calculate_rotation_matrix(ut1);

        let x =
            rotation_matrix.m0 * pav.x + rotation_matrix.m1 * pav.y + rotation_matrix.m2 * pav.z;
        let y =
            rotation_matrix.m3 * pav.x + rotation_matrix.m4 * pav.y + rotation_matrix.m5 * pav.z;
        let z =
            rotation_matrix.m6 * pav.x + rotation_matrix.m7 * pav.y + rotation_matrix.m8 * pav.z;

        types::EcefPosition { x, y, z }
    }

    fn calculate_rotation_matrix(&self, ut1: f64) -> types::RotationMatrix {
        // precession (gm2000 to mod)
        let mut t = (ut1 - 0.5) / 36525.0;
        let zeta = t * (0.6406161 + t * (0.0000839 + 0.0000050 * t));
        let z = t * (0.6406161 + t * (0.0003041 + 0.0000051 * t));
        let theta = t * (0.5567530 - t * (0.0001185 + 0.0000116 * t));
        let mut a = self.calc_rz(-constants::PI * 0.5 - z * constants::DE2RA);
        let mut b = self.calc_rx(theta * constants::DE2RA);
        let mut c = self.calc_rz(constants::PI * 0.5 - zeta * constants::DE2RA);
        let m1 = self.multiply_matrix(a, self.multiply_matrix(b, c));

        // nutation (mod to tod)
        t = ut1 - 0.5;
        let a1 = (125.0 - 0.05295 * t) * constants::DE2RA;
        let b1 = (200.9 + 1.97129 * t) * constants::DE2RA;
        let dpsi = (-0.0048 * a1.sin() - 0.0004 * b1.sin()) * constants::DE2RA;
        let deps = (0.0026 * a1.cos() + 0.0002 * b1.cos()) * constants::DE2RA;
        let eps = 23.439291 * constants::DE2RA;
        let dmu = dpsi * eps.cos();
        let dnu = dpsi * eps.sin();
        a = self.calc_rz(-dmu);
        b = self.calc_rx(-deps);
        c = self.calc_ry(dnu);
        let m2 = self.multiply_matrix(a, self.multiply_matrix(b, c));

        // earth rotation (tod to pef)
        t = ut1;
        let g = 99.96779469 + t * (360.9856473662860 + 0.29079e-12 * t);
        let h = g * constants::DE2RA + dmu;
        let m3 = self.calc_rz(h);

        self.multiply_matrix(self.multiply_matrix(m3, m2), m1)
    }

    fn calc_rx(&self, theta: f64) -> types::RotationMatrix {
        let cs = theta.cos();
        let sn = theta.sin();

        types::RotationMatrix {
            m0: 1.0,
            m1: 0.0,
            m2: 0.0,
            m3: 0.0,
            m4: cs,
            m5: sn,
            m6: 0.0,
            m7: -sn,
            m8: cs,
        }
    }

    fn calc_ry(&self, theta: f64) -> types::RotationMatrix {
        let cs = theta.cos();
        let sn = theta.sin();

        types::RotationMatrix {
            m0: cs,
            m1: 0.0,
            m2: -sn,
            m3: 0.0,
            m4: 1.0,
            m5: 0.0,
            m6: sn,
            m7: 0.0,
            m8: cs,
        }
    }
    fn calc_rz(&self, theta: f64) -> types::RotationMatrix {
        let cs = theta.cos();
        let sn = theta.sin();

        types::RotationMatrix {
            m0: cs,
            m1: sn,
            m2: 0.0,
            m3: -sn,
            m4: cs,
            m5: 0.0,
            m6: 0.0,
            m7: 0.0,
            m8: 1.0,
        }
    }

    fn multiply_matrix(
        &self,
        first: types::RotationMatrix,
        second: types::RotationMatrix,
    ) -> types::RotationMatrix {
        let mut x = first.m0;
        let mut y = first.m1;
        let mut z = first.m2;
        let m0 = x * second.m0 + y * second.m3 + z * second.m6;
        let m1 = x * second.m1 + y * second.m4 + z * second.m7;
        let m2 = x * second.m2 + y * second.m5 + z * second.m8;

        x = first.m3;
        y = first.m4;
        z = first.m5;
        let m3 = x * second.m0 + y * second.m3 + z * second.m6;
        let m4 = x * second.m1 + y * second.m4 + z * second.m7;
        let m5 = x * second.m2 + y * second.m5 + z * second.m8;

        x = first.m6;
        y = first.m7;
        z = first.m8;
        let m6 = x * second.m0 + y * second.m3 + z * second.m6;
        let m7 = x * second.m1 + y * second.m4 + z * second.m7;
        let m8 = x * second.m2 + y * second.m5 + z * second.m8;

        types::RotationMatrix {
            m0,
            m1,
            m2,
            m3,
            m4,
            m5,
            m6,
            m7,
            m8,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_rotation_matrix() {
        let ground_track = GroundTrack::new();

        let ut1 = 0.5;

        let eterm = ground_track.calculate_rotation_matrix(ut1);

        let expected_output = types::RotationMatrix {
            m0: 0.18155964819521067,
            m1: -0.9833799334794167,
            m2: -2.40260147885129e-5,
            m3: 0.9833799334210463,
            m4: 0.18155964742246256,
            m5: 3.118733126373619e-5,
            m6: -2.630684096956415e-5,
            m7: -2.9289061715492326e-5,
            m8: 0.9999999992250505,
        };

        assert_eq!(eterm, expected_output);
    }

    #[test]
    fn test_multiply_matrix() {
        let first = types::RotationMatrix {
            m0: 1.0,
            m1: 0.0,
            m2: 1.0,
            m3: 0.0,
            m4: 1.0,
            m5: 0.0,
            m6: 1.0,
            m7: 0.0,
            m8: 1.0,
        };

        let second = types::RotationMatrix {
            m0: 1.0,
            m1: 1.0,
            m2: 1.0,
            m3: 2.0,
            m4: 2.0,
            m5: 2.0,
            m6: 3.0,
            m7: 3.0,
            m8: 3.0,
        };

        let ground_track = GroundTrack::new();

        let result = ground_track.multiply_matrix(first, second);

        assert_eq!(
            result,
            types::RotationMatrix {
                m0: 4.0,
                m1: 4.0,
                m2: 4.0,
                m3: 2.0,
                m4: 2.0,
                m5: 2.0,
                m6: 4.0,
                m7: 4.0,
                m8: 4.0,
            }
        )
    }
}
