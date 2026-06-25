use eframe::egui;
fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([320.0, 240.0]),
        ..Default::default()
    };
    eframe::run_ui_native("Orbit propagator", options, move |ui, _frame| {
        egui::CentralPanel::default().show_inside(ui, |ui| {
            ui.heading("Orbit propagator");
        });
    })
}
