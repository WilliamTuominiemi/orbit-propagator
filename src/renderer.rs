use eframe::egui;
use eframe::egui::Color32;
use egui_plot::Line;
use egui_plot::Plot;
use egui_plot::PlotPoints;

pub struct Renderer {
    pub eo: f64,
    pub xno: f64,
    pub xmo: f64,
    pub xincl: f64,
    pub xnodeo: f64,
    pub omegao: f64,
    pub bstar: f64,
}

impl Renderer {
    pub fn render(&mut self, points: Vec<[f64; 2]>) -> eframe::Result {
        let window_size = egui::vec2(1400.0, 600.0);
        let control_pane_width = 200.0;

        let options = eframe::NativeOptions {
            viewport: egui::ViewportBuilder::default().with_inner_size(window_size),
            ..Default::default()
        };

        let mut eo = self.eo.to_string();
        let mut xno = self.xno.to_string();
        let mut xmo = self.xmo.to_string();
        let mut xincl = self.xincl.to_string();
        let mut xnodeo = self.xnodeo.to_string();
        let mut omegao = self.omegao.to_string();
        let mut bstar = self.bstar.to_string();

        eframe::run_ui_native("Janus", options, move |ui, _frame| {
            egui_extras::install_image_loaders(ui.ctx());

            egui::Panel::left("my_left_panel")
                .exact_size(control_pane_width)
                .show_inside(ui, |ui| {
                    ui.label(
                        egui::RichText::new("NORAD SPACETRACK REPORT No.3 SGP4 Sample test case")
                            .underline(),
                    );

                    ui.label("Eccentricity (EO)");
                    let eo_input = ui.add(egui::TextEdit::singleline(&mut eo));

                    ui.label("Mean Motion (XNO)");
                    let xno_input = ui.add(egui::TextEdit::singleline(&mut xno));

                    ui.label("Mean Anomaly (XMO)");
                    let xno_input = ui.add(egui::TextEdit::singleline(&mut xmo));

                    ui.label("Inclination (XINCL)");
                    let xno_input = ui.add(egui::TextEdit::singleline(&mut xincl));

                    ui.label("Right Ascension of the Ascending Node (XNODEO)");
                    let xno_input = ui.add(egui::TextEdit::singleline(&mut xnodeo));

                    ui.label("Argument of Perigee (OMEGAO)");
                    let xno_input = ui.add(egui::TextEdit::singleline(&mut omegao));

                    ui.label("B-Star Drag Term (BSTAR)");
                    let xno_input = ui.add(egui::TextEdit::singleline(&mut bstar));
                });

            egui::CentralPanel::default().show_inside(ui, |ui| {
                let panel_rect = ui.available_rect_before_wrap();

                egui::Image::new(egui::include_image!(".././images/map.png"))
                    .paint_at(ui, panel_rect);

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
}
