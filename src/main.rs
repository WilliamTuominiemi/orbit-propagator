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

fn compute_points(tle: &types::TLE, t_until: i32) -> Vec<types::GraphDataPoint> {
    let sgp4 = sgp4::Sgp4::new(tle);

    let test_epoch = -7030.01291535; // Spacetrack Report No. 3 base epoch
    let gt = ground_track::GroundTrack::new(test_epoch);
    calculate_points(&sgp4, &gt, t_until)
}

fn main() -> eframe::Result {
    let renderer = renderer::Renderer::new(compute_points);

    renderer.run()
}
