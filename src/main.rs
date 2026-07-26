mod constants;
mod ground_track;
mod helpers;
mod renderer;
mod sgp4;
mod test_constants;
mod types;

fn calculate_points(
    sgp4: &sgp4::Sgp4,
    gt: &ground_track::GroundTrack,
    t_until: i32,
) -> Vec<types::GraphDataPoint> {
    let mut last_lon = 0.0;

    let mut data_points = Vec::new();

    for i in 0..t_until {
        let tsince = i as f64;
        let pav = sgp4.propagate(tsince);
        let pav_m = types::PositionAndVelocity {
            x: pav.x * 1000.0,
            y: pav.y * 1000.0,
            z: pav.z * 1000.0,
            xdot: pav.xdot * 1000.0,
            ydot: pav.ydot * 1000.0,
            zdot: pav.zdot * 1000.0,
        };
        let geodetic = gt.eci_to_geodetic(tsince, pav_m);

        let lon = geodetic.lon.to_degrees();
        let lat = geodetic.lat.to_degrees();

        if i > 0 && (lon - last_lon).abs() > 180.0 {
            data_points.push(types::GraphDataPoint {
                point: [f64::NAN, f64::NAN],
                altitude: 0.0,
                velocity: 0.0,
            });
        }

        data_points.push(types::GraphDataPoint {
            point: [lon, lat],
            altitude: geodetic.alt,
            velocity: (geodetic.vel_east.powi(2)
                + geodetic.vel_north.powi(2)
                + geodetic.vel_up.powi(2))
            .sqrt(),
        });
        last_lon = lon;
    }

    data_points
}

fn compute_points(
    eo: f64,
    bstar: f64,
    xincl: f64,
    omegao: f64,
    xmo: f64,
    xno: f64,
    xnodeo: f64,
    t_until: i32,
) -> Vec<types::GraphDataPoint> {
    let norad_tle = types::TLE {
        name: "SGP4".to_string(),
        number: 88888,
        international_designator: "SGP4".to_string(),
        epoch_year_julian_fraction: 80275.98708465,
        first_derivative_of_mean_motion: 0.00073094,
        second_derivative_of_mean_motion: (13844_f64).powi(-3),
        drag_term: (66816_f64).powi(-4),
        ephemeris_type: 0,
        element_number_check_sum: 8,
        inclination: 72.8435,
        right_ascension_of_ascending_node: 115.9689,
        eccentricity: 0.0086731,
        argument_of_perigee: 52.6988,
        mean_anomaly: 110.5714,
        mean_motion: 16.05824518,
        revolution_number_check_sum: 105,
    };

    let sgp4 = sgp4::Sgp4::new(
        eo,
        bstar,
        xincl,
        omegao,
        xmo,
        xno,
        xnodeo,
        test_constants::E6A,
    );

    let test_epoch = -7030.01291535; // Spacetrack Report No. 3 base epoch
    let gt = ground_track::GroundTrack::new(test_epoch);
    calculate_points(&sgp4, &gt, t_until)
}

fn main() -> eframe::Result {
    let renderer = renderer::Renderer::new(compute_points);

    renderer.run()
}
