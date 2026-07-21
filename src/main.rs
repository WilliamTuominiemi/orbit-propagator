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
    t_since: i32,
) -> Vec<[f64; 2]> {
    let mut points = Vec::new();
    let mut last_lon = 0.0;

    for i in 0..t_since {
        let tsince = i as f64 * 0.01;
        let pav = sgp4.propagate(tsince);
        let geodetic = gt.eci_to_geodetic(tsince, pav);

        let lon = geodetic.lon.to_degrees();
        let lat = geodetic.lat.to_degrees();

        if i > 0 && (lon - last_lon).abs() > 180.0 {
            points.push([f64::NAN, f64::NAN]);
        }

        points.push([lon, lat]);
        last_lon = lon;
    }

    points
}

fn compute_points(
    eo: f64,
    bstar: f64,
    xincl: f64,
    omegao: f64,
    xmo: f64,
    xno: f64,
    xnodeo: f64,
    t_since: i32,
) -> Vec<[f64; 2]> {
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
    calculate_points(&sgp4, &gt, t_since)
}

fn main() -> eframe::Result {
    let renderer = renderer::Renderer::new(compute_points);

    renderer.run()
}
