use eframe::egui;
use egui_plot::Line;
use egui_plot::Plot;
use egui_plot::PlotPoints;

mod constants;
mod helpers;
mod sgp4;
mod test_constants;
mod types;

fn main() -> eframe::Result {
    let sgp4 = sgp4::Sgp4::new(
        test_constants::EO,
        test_constants::BSTAR,
        test_constants::XINCL,
        test_constants::OMEGAO,
        constants::CK4,
        test_constants::XMO,
        test_constants::XNO,
        test_constants::XNODEO,
        test_constants::E6A,
    );

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([500.0, 240.0]),
        ..Default::default()
    };
    eframe::run_ui_native("Janus", options, move |ctx, _frame| {
        egui::CentralPanel::default().show_inside(ctx, |ui| {
            let orbit: PlotPoints = (0..9000)
                .map(|i| {
                    let tsince = i as f64 * 0.01;
                    [tsince, sgp4.propagate(tsince).x]
                })
                .collect();
            let line = Line::new("orbit", orbit);
            Plot::new("my_plot")
                .view_aspect(2.0)
                .show(ui, |plot_ui| plot_ui.line(line));
        });
    })
}
