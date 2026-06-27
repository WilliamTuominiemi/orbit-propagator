use eframe::egui;
use egui_plot::Line;
use egui_plot::Plot;
use egui_plot::PlotPoints;

fn recover_original_mean_motion_and_semimajor_axis(
    xke: f64,
    xno: f64,
    tothrd: f64,
    xincl: f64,
    eo: f64,
    ck2: f64,
) -> (f64, f64) {
    let a1 = (xke / xno).powf(tothrd);
    let cosio = xincl.cos();
    let theta2 = cosio * cosio;
    let x3thm1 = 3.0 * theta2 - 1.0;
    let eosq = eo * eo;
    let betao2 = 1.0 - eosq;
    let betao = f64::sqrt(betao2);
    let del1 = 1.5 * ck2 * x3thm1 / (a1 * a1 * betao * betao2);
    let ao = a1 * (1.0 - del1 * (0.5 * tothrd + del1 * (1.0 + 134.0 / 81.0 * del1)));
    let delo = 1.5 * ck2 * x3thm1 / (ao * ao * betao * betao2);

    println!("xno {:?} | ao {:?} | delo {:?}", xno, ao, delo);

    let xnodp = xno / (1.0 + delo);
    let aodp = ao / (1.0 - delo);

    (xnodp, aodp)
}

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([500.0, 240.0]),
        ..Default::default()
    };
    eframe::run_ui_native("Latveria", options, move |ctx, _frame| {
        egui::CentralPanel::default().show_inside(ctx, |ui| {
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

#[cfg(test)]
mod tests {
    use super::*;

    // Values from NORAD SPACETRACK REPORT NO. 3
    // From constants and sample test case parameters
    const DE2RA: f64 = 0.0174532925;
    const XKE: f64 = 0.074366916;
    const TWOPI: f64 = 6.2831853;
    const XMNPDA: f64 = 1440.0;
    const XNO: f64 = 16.05824518 * (TWOPI / XMNPDA);
    const TOTHRD: f64 = 0.66666667;
    const XINCL: f64 = 72.8435 * DE2RA;
    const EO: f64 = 0.0086731;
    const CK2: f64 = 0.0005413080;

    #[test]
    fn test_recover_original_mean_motion_and_semimajor_axis() {
        let (xnodp, aodp) =
            recover_original_mean_motion_and_semimajor_axis(XKE, XNO, TOTHRD, XINCL, EO, CK2);

        assert_eq!(xnodp, 0.07010615558630984);
        assert_eq!(aodp, 1.040117522759639);
    }
}
