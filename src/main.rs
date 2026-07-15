use eframe::egui;
use eframe::egui::Color32;
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

    render(points)
}

fn render(points: Vec<[f64; 2]>) -> eframe::Result {
    let window_size = egui::vec2(1400.0, 600.0);
    let control_pane_width = 200.0;

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size(window_size),
        ..Default::default()
    };
    eframe::run_ui_native("Janus", options, move |ui, _frame| {
        egui_extras::install_image_loaders(ui.ctx());

        egui::Panel::left("my_left_panel")
            .exact_size(control_pane_width)
            .show_inside(ui, |ui| {
                ui.label(
                    egui::RichText::new("NORAD SPACETRACK REPORT No.3 SGP4 Sample test case")
                        .underline(),
                );
                ui.label(format!("Eccentricity (EO): {}", test_constants::EO));
                ui.label(format!("Mean Motion (XNO): {}", test_constants::XNO));
                ui.label(format!("Mean Anomaly (XMO): {}", test_constants::XMO));
                ui.label(format!("Inclination (XINCL): {}", test_constants::XINCL));
                ui.label(format!(
                    "Right Ascension of the Ascending Node (XNODEO): {}",
                    test_constants::XNODEO
                ));
                ui.label(format!(
                    "Argument of Perigee (OMEGAO): {}",
                    test_constants::OMEGAO
                ));
                ui.label(format!(
                    "B-Star Drag Term (BSTAR): {}",
                    test_constants::BSTAR
                ));
            });

        egui::CentralPanel::default().show_inside(ui, |ui| {
            let panel_rect = ui.available_rect_before_wrap();

            egui::Image::new(egui::include_image!(".././images/map.png")).paint_at(ui, panel_rect);

            let orbit = PlotPoints::new(points.clone());
            let line = Line::new("orbit", orbit).width(3.0).color(Color32::ORANGE);

            Plot::new("my_plot")
                .show_background(false)
                .allow_drag(false)
                .allow_zoom(false)
                .allow_scroll(false)
                .grid_color(Color32::WHITE)
                .show_axes(false)
                .show(ui, |plot_ui| plot_ui.line(line));
        });
    })
}
