use eframe::egui;
use eframe::egui::Color32;
use egui::Vec2;
use egui_plot::Line;
use egui_plot::{Plot, PlotPoints, Points};

use crate::test_constants;
use crate::types;

pub struct Renderer {
    pub eo: f64,
    pub xno: f64,
    pub xmo: f64,
    pub xincl: f64,
    pub xnodeo: f64,
    pub omegao: f64,
    pub bstar: f64,
    pub t_until: i32,
    eo_str: String,
    xno_str: String,
    xmo_str: String,
    xincl_str: String,
    xnodeo_str: String,
    omegao_str: String,
    bstar_str: String,
    t_until_str: String,
    points: Vec<[f64; 2]>,
    pub compute_points: fn(f64, f64, f64, f64, f64, f64, f64, i32) -> types::GraphData,
    t_since: i32,
    altitude: f64,
    velocity: f64,
}

impl Renderer {
    pub fn new(
        compute_points: fn(f64, f64, f64, f64, f64, f64, f64, i32) -> types::GraphData,
    ) -> Self {
        let eo = test_constants::EO;
        let xno = test_constants::XNO;
        let xmo = test_constants::XMO;
        let xincl = test_constants::XINCL;
        let xnodeo = test_constants::XNODEO;
        let omegao = test_constants::OMEGAO;
        let bstar = test_constants::BSTAR;
        let t_until = 27000;
        let t_since = t_until / 3;

        let gd = compute_points(eo, bstar, xincl, omegao, xmo, xno, xnodeo, t_until);

        Self {
            eo_str: eo.to_string(),
            xno_str: xno.to_string(),
            xmo_str: xmo.to_string(),
            xincl_str: xincl.to_string(),
            xnodeo_str: xnodeo.to_string(),
            omegao_str: omegao.to_string(),
            bstar_str: bstar.to_string(),
            t_until_str: t_until.to_string(),
            eo,
            xno,
            xmo,
            xincl,
            xnodeo,
            omegao,
            bstar,
            t_until,
            points: gd.points,
            compute_points,
            t_since,
            altitude: gd.altitude,
            velocity: gd.velocity,
        }
    }

    pub fn run(self) -> eframe::Result {
        let window_size = egui::vec2(1400.0, 650.0);

        let options = eframe::NativeOptions {
            viewport: egui::ViewportBuilder::default().with_inner_size(window_size),
            ..Default::default()
        };

        eframe::run_native(
            "Janus",
            options,
            Box::new(|cc| {
                egui_extras::install_image_loaders(&cc.egui_ctx);
                Ok(Box::new(self))
            }),
        )
    }
}

impl eframe::App for Renderer {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let control_pane_width = 200.0;

        egui::Panel::left("control_pane")
            .exact_size(control_pane_width)
            .show_inside(ui, |ui| {
                egui::Frame::new()
                    .inner_margin(egui::Margin::same(5))
                    .show(ui, |ui| {
                        ui.label("Eccentricity (EO)");
                        ui.add(egui::TextEdit::singleline(&mut self.eo_str));

                        ui.label("Mean Motion (XNO)");
                        ui.add(egui::TextEdit::singleline(&mut self.xno_str));

                        ui.label("Mean Anomaly (XMO)");
                        ui.add(egui::TextEdit::singleline(&mut self.xmo_str));

                        ui.label("Inclination (XINCL)");
                        ui.add(egui::TextEdit::singleline(&mut self.xincl_str));

                        ui.label("Right Ascension of the Ascending Node (XNODEO)");
                        ui.add(egui::TextEdit::singleline(&mut self.xnodeo_str));

                        ui.label("Argument of Perigee (OMEGAO)");
                        ui.add(egui::TextEdit::singleline(&mut self.omegao_str));

                        ui.label("B-Star Drag Term (BSTAR)");
                        ui.add(egui::TextEdit::singleline(&mut self.bstar_str));

                        ui.label("Tracking time");
                        ui.add(egui::TextEdit::singleline(&mut self.t_until_str));
                    });

                egui::Frame::new()
                    .inner_margin(egui::Margin::same(5)) // adjust margin amount as needed
                    .show(ui, |ui| {
                        ui.with_layout(egui::Layout::top_down(egui::Align::Center), |ui| {
                            if ui
                                .add(
                                    egui::Button::new("Update graph")
                                        .min_size(Vec2::new(
                                            control_pane_width * 0.9,
                                            control_pane_width * 0.9 / 5.0,
                                        ))
                                        .stroke(egui::Stroke::new(2.0, egui::Color32::ORANGE)),
                                )
                                .clicked()
                            {
                                let parsed = (
                                    self.eo_str.parse::<f64>(),
                                    self.xno_str.parse::<f64>(),
                                    self.xmo_str.parse::<f64>(),
                                    self.xincl_str.parse::<f64>(),
                                    self.xnodeo_str.parse::<f64>(),
                                    self.omegao_str.parse::<f64>(),
                                    self.bstar_str.parse::<f64>(),
                                    self.t_until_str.parse::<i32>(),
                                );

                                if let (
                                    Ok(eo),
                                    Ok(xno),
                                    Ok(xmo),
                                    Ok(xincl),
                                    Ok(xnodeo),
                                    Ok(omegao),
                                    Ok(bstar),
                                    Ok(t_until),
                                ) = parsed
                                {
                                    self.eo = eo;
                                    self.xno = xno;
                                    self.xmo = xmo;
                                    self.xincl = xincl;
                                    self.xnodeo = xnodeo;
                                    self.omegao = omegao;
                                    self.bstar = bstar;
                                    self.t_until = t_until;

                                    let gd = (self.compute_points)(
                                        eo,
                                        bstar,
                                        xincl,
                                        omegao,
                                        xmo,
                                        xno,
                                        xnodeo,
                                        self.t_until,
                                    );
                                    self.points = gd.points;
                                }
                            }
                        });
                    });

                egui::Panel::bottom("metrics")
                    .frame(egui::Frame::default().outer_margin(12.6))
                    .show_inside(ui, |ui| {
                        ui.label(format!("Altitude: {} km", self.altitude / 1000.0));
                        ui.label(format!("Velocity: {} km/s", self.velocity / 1000.0));
                    });
            });

        egui::CentralPanel::default().show_inside(ui, |ui| {
            let mut panel_rect = ui.available_rect_before_wrap();
            panel_rect.set_height(panel_rect.width() / 2.0);

            egui::Image::new(egui::include_image!(".././images/map.png")).paint_at(ui, panel_rect);

            let orbit = PlotPoints::new(
                self.points
                    .iter()
                    .take(self.t_since as usize)
                    .copied()
                    .collect::<Vec<_>>(),
            );
            let line = Line::new("orbit", orbit).width(4.0).color(Color32::ORANGE);

            let current_position = self
                .points
                .iter()
                .take(self.t_since as usize)
                .last()
                .copied()
                .unwrap_or([0.0, 0.0]);
            let points = Points::new("current_position", PlotPoints::new(vec![current_position]))
                .radius(12.0)
                .shape(egui_plot::MarkerShape::Diamond)
                .color(egui::Color32::LIGHT_RED);

            let max_x = 180.0;
            let max_y = 80.0;

            Plot::new("orbit_plot")
                .show_background(false)
                .allow_drag(false)
                .allow_zoom(false)
                .allow_scroll(false)
                .grid_color(Color32::WHITE)
                .show_axes(false)
                .width(panel_rect.width())
                .height(panel_rect.height())
                .include_x(-max_x)
                .include_x(max_x)
                .include_y(-max_y)
                .include_y(max_y)
                .show(ui, |plot_ui| {
                    plot_ui.line(line);

                    plot_ui.points(points);
                });

            egui::Panel::bottom("time_slider_pane").show_inside(ui, |ui| {
                ui.style_mut().spacing.slider_width = panel_rect.width() - 100.0;
                ui.label("Minutes since start:");
                ui.add(egui::Slider::new(&mut self.t_since, 0..=self.t_until));
            });
        });
    }
}
