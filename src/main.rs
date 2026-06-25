use eframe::egui;
use egui_plot::Line;
use egui_plot::Plot;
use egui_plot::PlotPoints;

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([320.0, 240.0]),
        ..Default::default()
    };
    eframe::run_ui_native("Orbit propagator", options, move |ctx, _frame| {
        egui::CentralPanel::default().show_inside(ctx, |ui| {
            ui.heading("Orbit propagator");

            let sin: PlotPoints = (0..1000)
                .map(|i| {
                    let x = i as f64 * 0.01;
                    [x, x.sin()]
                })
                .collect();
            let line = Line::new("sin", sin);
            Plot::new("my_plot")
                .view_aspect(2.0)
                .show(ui, |plot_ui| plot_ui.line(line));
        });
    })
}
