use crate::constants;
use crate::types;

pub struct GroundTrack {
    pub epoch_year_julian_fraction: String,
}

impl GroundTrack {
    pub fn new(epoch_year_julian_fraction: String) -> Self {
        GroundTrack {
            epoch_year_julian_fraction,
        }
    }

    pub fn eci_to_geodetic(
        &self,
        tsince: f64,
        pav: types::PositionAndVelocity,
    ) -> types::GeodeticPosition {
        let ecef = self.eci_to_ecef(self.tsince_to_ut1(tsince), pav);

        let a = constants::XKMPER * 1000.0;
        let f: f64 = 1.0 / 298.257223563;
        let b = a - f * a;

        let e2 = f * (2.0 - f);

        let clambda = ecef.y.atan2(ecef.x);
        let p = (ecef.x * ecef.x + ecef.y * ecef.y).sqrt();

        let mut h_old = 0.0;
        let mut theta = ecef.z.atan2(p * (1.0 - e2));
        let mut cs = theta.cos();
        let mut sn = theta.sin();

        let mut n = (a * a) / ((a * cs) * (a * cs) + (b * sn) * (b * sn)).sqrt();
        let mut h = p / cs - n;
        let mut iterations = 0;

        while (h - h_old).abs() > 1.0e-6 && iterations < 100 {
            h_old = h;
            theta = ecef.z.atan2(p * (1.0 - e2 * n / (n + h)));
            cs = theta.cos();
            sn = theta.sin();
            n = (a * a) / ((a * cs) * (a * cs) + (b * sn) * (b * sn)).sqrt();
            h = p / cs - n;
            iterations += 1;
        }

        let (vel_east, vel_north, vel_up) = self.ecef_velocity_to_enu(&ecef, theta, clambda);

        types::GeodeticPosition {
            lat: theta,
            lon: clambda,
            alt: h,
            vel_east,
            vel_north,
            vel_up,
        }
    }

    fn eci_to_ecef(&self, ut1: f64, pav: types::PositionAndVelocity) -> types::EcefPosition {
        let rotation_matrix = self.calculate_rotation_matrix(ut1);

        let x =
            rotation_matrix.m0 * pav.x + rotation_matrix.m1 * pav.y + rotation_matrix.m2 * pav.z;
        let y =
            rotation_matrix.m3 * pav.x + rotation_matrix.m4 * pav.y + rotation_matrix.m5 * pav.z;
        let z =
            rotation_matrix.m6 * pav.x + rotation_matrix.m7 * pav.y + rotation_matrix.m8 * pav.z;

        let vx_rot = rotation_matrix.m0 * pav.xdot
            + rotation_matrix.m1 * pav.ydot
            + rotation_matrix.m2 * pav.zdot;
        let vy_rot = rotation_matrix.m3 * pav.xdot
            + rotation_matrix.m4 * pav.ydot
            + rotation_matrix.m5 * pav.zdot;
        let vz_rot = rotation_matrix.m6 * pav.xdot
            + rotation_matrix.m7 * pav.ydot
            + rotation_matrix.m8 * pav.zdot;

        let xdot = vx_rot + constants::EARTH_ROTATION_RATE * y;
        let ydot = vy_rot - constants::EARTH_ROTATION_RATE * x;
        let zdot = vz_rot;

        types::EcefPosition {
            x,
            y,
            z,
            xdot,
            ydot,
            zdot,
        }
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
        let g = 99.96779469 + t * (360.985_647_366_286 + 0.29079e-12 * t);
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

    fn tsince_to_ut1(&self, tsince: f64) -> f64 {
        let ut1 = self.epoch_into_days_since_J2000();
        ut1 + (tsince / 1440.0)
    }

    fn ecef_velocity_to_enu(
        &self,
        ecef: &types::EcefPosition,
        lat: f64,
        lon: f64,
    ) -> (f64, f64, f64) {
        let sin_lat = lat.sin();
        let cos_lat = lat.cos();
        let sin_lon = lon.sin();
        let cos_lon = lon.cos();

        let east = -sin_lon * ecef.xdot + cos_lon * ecef.ydot;
        let north =
            -sin_lat * cos_lon * ecef.xdot - sin_lat * sin_lon * ecef.ydot + cos_lat * ecef.zdot;
        let up =
            cos_lat * cos_lon * ecef.xdot + cos_lat * sin_lon * ecef.ydot + sin_lat * ecef.zdot;

        (east, north, up)
    }

    fn epoch_into_days_since_J2000(&self) -> f64 {
        let year_str = &self.epoch_year_julian_fraction[0..2];
        let day_of_year_str = &self.epoch_year_julian_fraction[2..];

        let mut year: i32 = year_str.parse().unwrap();
        let day_of_year: f64 = day_of_year_str.parse().unwrap();

        let j2000 = 2451545.0;

        if year < 57 {
            year += 2000;
        } else {
            year += 1900;
        }

        let y = (year - 1) as f64;
        let A = (y / 100.0_f64).floor();
        let B = 2.0_f64 - A + (A / 4.0_f64).floor();

        let jd_jan_0 =
            (365.25_f64 * (y + 4716.0_f64)).floor() + (30.6001_f64 * 14.0_f64).floor() + B
                - 1524.5_f64;

        let jd = jd_jan_0 + day_of_year;

        jd - j2000
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::helpers;

    #[test]
    fn test_eci_to_geodetic() {
        let ground_track = GroundTrack::new("80275.98708465".to_string());

        let pav = types::PositionAndVelocity {
            x: 4263871.9243,
            y: 722591.1075,
            z: 4672986.8878,
            xdot: 0.0,
            ydot: 0.0,
            zdot: 0.0,
        };

        let tsince = 0.0;

        let geodetic = ground_track.eci_to_geodetic(tsince, pav);

        helpers::assert_approx(geodetic.lat.to_degrees(), 47.301461267437865);
        helpers::assert_approx(geodetic.lon.to_degrees(), -176.36028751239033);
        helpers::assert_approx(geodetic.alt, 438.2550053251907);
    }

    #[test]
    fn test_calculate_rotation_matrix() {
        let ut1 = "80275.98708465".to_string();

        let ground_track = GroundTrack::new(ut1);

        let eterm = ground_track.calculate_rotation_matrix(-7030.512915350031);

        let expected_output = types::RotationMatrix {
            m0: -0.9945905437510192,
            m1: -0.10385609314679319,
            m2: -0.0018873784295872415,
            m3: 0.10385597135741614,
            m4: -0.9945923339469838,
            m5: 0.00016268773489335004,
            m6: -0.001894068229873292,
            m7: -3.4207837414618897e-5,
            m8: 0.9999982056660726,
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

        let ground_track = GroundTrack::new("80275.98708465".to_string());

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

    #[test]
    fn test_epoch_into_days_since_J2000() {
        let epoch_year_julian_fraction = "80275.98708465".to_string();
        let ground_track = GroundTrack::new(epoch_year_julian_fraction);
        let days_since_j2000 = ground_track.epoch_into_days_since_J2000();

        assert_eq!(days_since_j2000, -7030.512915350031);
    }
}
