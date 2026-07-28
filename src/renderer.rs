use eframe::egui;
use eframe::egui::Color32;
use egui::Vec2;
use egui_plot::Line;
use egui_plot::{Plot, PlotBounds, PlotPoints, Points};

use crate::{helpers, types};

pub struct Renderer {
    t_until: i32,
    tle: types::TLE,
    tle_input: String,
    t_until_str: String,
    data_points: Vec<types::GraphDataPoint>,
    pub compute_points: fn(&types::TLE, i32) -> Vec<types::GraphDataPoint>,
    t_since: i32,
}

impl Renderer {
    pub fn new(compute_points: fn(&types::TLE, i32) -> Vec<types::GraphDataPoint>) -> Self {
        let t_until = 900;
        let t_since = t_until / 10;

        let tle_input = "SGP4 (SGP4)
1 88888U 98067A 80275.98708465 .00073094 13844-3 66816-4 0 8
2 88888 72.8435 115.9689 0086731 52.6988 110.5714 16.05824518 105"
            .to_string();
        let tle = helpers::text_to_tle(&tle_input);

        let data_points = compute_points(&tle, t_until);

        Self {
            tle,
            tle_input,
            t_until_str: t_until.to_string(),
            t_until,
            data_points,
            compute_points,
            t_since,
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
        egui::Panel::left("control_pane").show_inside(ui, |ui| {
            egui::Frame::new()
                .inner_margin(egui::Margin::same(5))
                .show(ui, |ui| {
                    let current_data_point = self
                        .data_points
                        .iter()
                        .take(self.t_since as usize)
                        .last()
                        .unwrap_or(&types::GraphDataPoint {
                            point: [f64::NAN, f64::NAN],
                            altitude: 0.0,
                            velocity: 0.0,
                        });

                    ui.label(format!(
                        "Altitude: {:.3} km",
                        current_data_point.altitude / 1000.0
                    ));
                    ui.label(format!(
                        "Velocity: {:.3} km/s",
                        current_data_point.velocity / 1000.0
                    ));
                    ui.label(format!("Latitude: {:.6}°", current_data_point.point[1]));
                    ui.label(format!("Longitude: {:.6}°", current_data_point.point[0]));
                });

            ui.separator();

            egui::Frame::new()
                .inner_margin(egui::Margin::same(5))
                .show(ui, |ui| {
                    ui.label("Two-line element set");
                    ui.add(egui::TextEdit::multiline(&mut self.tle_input));
                });

            egui::Panel::bottom("update_button")
                .frame(egui::Frame::default().outer_margin(12.6))
                .show_inside(ui, |ui| {
                    ui.with_layout(egui::Layout::top_down(egui::Align::Center), |ui| {
                        if ui
                            .add(
                                egui::Button::new("Update graph")
                                    .min_size(Vec2::new(200.0, 40.0))
                                    .stroke(egui::Stroke::new(2.0, egui::Color32::ORANGE)),
                            )
                            .clicked()
                        {
                            let parsed = (self.t_until_str.parse::<i32>(),);

                            if let (Ok(t_until),) = parsed {
                                self.t_until = t_until;
                                self.tle = helpers::text_to_tle(&self.tle_input);

                                self.data_points = (self.compute_points)(&self.tle, self.t_until);
                            }
                        }
                    });
                });
        });

        egui::CentralPanel::default().show_inside(ui, |ui| {
            let current_data_point = self
                .data_points
                .iter()
                .take(self.t_since as usize)
                .last()
                .unwrap_or(&types::GraphDataPoint {
                    point: [f64::NAN, f64::NAN],
                    altitude: 0.0,
                    velocity: 0.0,
                });

            let mut panel_rect = ui.available_rect_before_wrap();
            panel_rect.set_height(panel_rect.width() / 2.0);

            let offset_x = 30.0;
            let image_rect = panel_rect.translate(egui::vec2(offset_x, 0.0));

            egui::Image::new(egui::include_image!(".././images/map.png")).paint_at(ui, image_rect);

            let orbit = PlotPoints::new(
                self.data_points
                    .iter()
                    .take(self.t_since as usize)
                    .map(|dp| dp.point)
                    .collect::<Vec<_>>(),
            );
            let line = Line::new("orbit", orbit).width(4.0).color(Color32::ORANGE);

            let points = Points::new(
                "current_position",
                PlotPoints::new(vec![current_data_point.point]),
            )
            .radius(12.0)
            .shape(egui_plot::MarkerShape::Asterisk)
            .color(egui::Color32::ORANGE);

            let max_x = 170.0;
            let max_y = 80.0;

            Plot::new("orbit_plot")
                .show_background(false)
                .allow_drag(false)
                .allow_zoom(false)
                .allow_scroll(false)
                .grid_color(Color32::WHITE)
                // .show_axes(false)
                .width(panel_rect.width())
                .height(panel_rect.height())
                .show(ui, |plot_ui| {
                    plot_ui.set_plot_bounds(PlotBounds::from_min_max(
                        [-max_x, -max_y],
                        [max_x, max_y],
                    ));

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
