use eframe::egui;
use egui_plot::Line;
use egui_plot::Plot;
use egui_plot::PlotPoints;

mod constants;
mod ground_track;
mod helpers;
mod sgp4;
mod test_constants;
mod types;

fn init_points(sgp4: &sgp4::Sgp4, gt: &ground_track::GroundTrack) -> Vec<[f64; 2]> {
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
    let sgp4 = sgp4::Sgp4::new(
        test_constants::EO,
        test_constants::BSTAR,
        test_constants::XINCL,
        test_constants::OMEGAO,
        test_constants::XMO,
        test_constants::XNO,
        test_constants::XNODEO,
        test_constants::E6A,
    );

    // Initialize with the Spacetrack Report No. 3 base epoch
    let test_epoch = -7030.01291535;
    let gt = ground_track::GroundTrack::new(test_epoch);
    let points = init_points(&sgp4, &gt);

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([500.0, 240.0]),
        ..Default::default()
    };
    eframe::run_ui_native("Janus", options, move |ctx, _frame| {
        egui::CentralPanel::default().show_inside(ctx, |ui| {
            let orbit = PlotPoints::new(points.clone());
            let line = Line::new("orbit", orbit);

            Plot::new("my_plot")
                .view_aspect(2.0)
                .show(ui, |plot_ui| plot_ui.line(line));
        });
    })
}
