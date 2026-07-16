mod constants;
mod ground_track;
mod helpers;
mod renderer;
mod sgp4;
mod test_constants;
mod types;

fn calculate_points(sgp4: &sgp4::Sgp4, gt: &ground_track::GroundTrack) -> Vec<[f64; 2]> {
    let mut points = Vec::new();
    let mut last_lon = 0.0;

    for i in 0..9000 {
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

fn main() -> eframe::Result {
    let eo = test_constants::EO;
    let xno = test_constants::XNO;
    let xmo = test_constants::XMO;
    let xincl = test_constants::XINCL;
    let xnodeo = test_constants::XNODEO;
    let omegao = test_constants::OMEGAO;
    let bstar = test_constants::BSTAR;

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
    let points = calculate_points(&sgp4, &gt);

    let mut renderer = renderer::Renderer {
        eo,
        bstar,
        xincl,
        omegao,
        xmo,
        xno,
        xnodeo,
    };

    renderer.render(points)
}
