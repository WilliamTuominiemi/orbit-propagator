use eframe::egui;
use egui_plot::Line;
use egui_plot::Plot;
use egui_plot::PlotPoints;

struct mean_motion_and_semimajor_axis_output {
    xnodp: f64,
    aodp: f64,
    betao2: f64,
    betao: f64,
    x3thm1: f64,
    theta2: f64,
    cosio: f64,
}

fn recover_original_mean_motion_and_semimajor_axis(
    xke: f64,
    xno: f64,
    tothrd: f64,
    xincl: f64,
    eo: f64,
    ck2: f64,
) -> mean_motion_and_semimajor_axis_output {
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
    let xnodp = xno / (1.0 + delo);
    let aodp = ao / (1.0 - delo);

    mean_motion_and_semimajor_axis_output {
        xnodp,
        aodp,
        betao2,
        betao,
        x3thm1,
        theta2,
        cosio,
    }
}

fn adjust_atmospheric_drag_for_low_orbit(
    s: f64,
    qoms2t: f64,
    aodp: f64,
    eo: f64,
    ae: f64,
    xkmper: f64,
) -> (f64, f64) {
    let mut s4 = s;
    let mut qoms24 = qoms2t;
    let perigee = (aodp * (1.0 - eo) - ae) * xkmper;
    if perigee < 156.0 {
        s4 = perigee - 78.0;
        if perigee <= 98.0 {
            s4 = 20.0;
        }
        qoms24 = ((120.0 - s4) * ae / xkmper).powf(4.0);
        s4 = s4 / xkmper + ae;
    }

    (s4, qoms24)
}

fn sgp4(
    mmasmao: mean_motion_and_semimajor_axis_output,
    eo: f64,
    ae: f64,
    xkmper: f64,
    s: f64,
    qoms2t: f64,
    ck2: f64,
    bstar: f64,
    xincl: f64,
    xj3: f64,
    omegao: f64,
    ck4: f64,
    tothrd: f64,
    xmo: f64,
) {
    let (xnodp, aodp, betao2, betao, x3thm1, theta2, cosio) = (
        mmasmao.xnodp,
        mmasmao.aodp,
        mmasmao.betao2,
        mmasmao.betao,
        mmasmao.x3thm1,
        mmasmao.theta2,
        mmasmao.cosio,
    );

    let (s4, qoms24) = adjust_atmospheric_drag_for_low_orbit(s, qoms2t, aodp, eo, ae, xkmper);

    let pinvsq = 1.0 / (aodp * aodp * betao2 * betao2);
    let tsi = 1.0 / (aodp - s4);
    let eta = aodp * eo * tsi;
    let etasq = eta * eta;
    let eeta = eo * eta;

    let psisq = (1.0 - etasq).abs();
    let coef = qoms24 * tsi.powf(4.0);
    let coef1 = coef / psisq.powf(3.5);

    let c2 = coef1
        * xnodp
        * (aodp * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
            + 0.75 * ck2 * tsi / psisq * x3thm1 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    let c1 = bstar * c2;
    let sinio = xincl.sin();
    let a3ovk2 = -xj3 / ck2 * ae.powf(3.0);
    let c3 = coef * tsi * a3ovk2 * xnodp * ae * sinio / eo;
    let x1mth2 = 1.0 - theta2;
    let c4 = 2.0
        * xnodp
        * coef1
        * aodp
        * betao2
        * (eta * (2.0 + 0.5 * etasq) + eo * (0.5 + 2.0 * etasq)
            - 2.0 * ck2 * tsi / (aodp * psisq)
                * (-3.0 * x3thm1 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta))
                    + 0.75 * x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * (2.0 * omegao).cos()));
    let c5 = 2.0 * coef1 * aodp * betao2 * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
    let theta4 = theta2 * theta2;
    let temp1 = 3.0 * ck2 * pinvsq * xnodp;
    let temp2 = temp1 * ck2 * pinvsq * xnodp;
    let temp3 = 1.25 * ck4 * pinvsq * pinvsq * xnodp;

    let xmdot = xnodp
        + 0.5 * temp1 * betao * x3thm1
        + 0.625 * temp2 * betao * (13.0 * 78.0 * theta2 + 137.0 * theta4);
    let x1m5th = 1.0 - 5.0 * theta2;
    let omgdot = -0.5 * temp1 * x1m5th
        + 0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4)
        + temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);
    let xhdot1 = -temp1 * cosio;
    let xnodot =
        xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 * (3.0 - 7.0 * theta2)) * cosio;
    let omgcof = bstar * c3 * omegao.cos();
    let xmcof = tothrd * coef * bstar * ae / eeta;
    let xnodcf = 3.5 * betao2 * xhdot1 * c1;
    let t2cof = 1.5 * c1;
    let xlcof = 0.125 * a3ovk2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
    let aycof = 0.25 * a3ovk2 * sinio;
    let delmo = (1.0 + eta * xmo.cos()).powf(3.0);
    let sinmo = xmo.sin();
    let x7thm1 = 7.0 * theta2 - 1.0;
    let c1sq = c1 * c1;
    let d2 = 4.0 * aodp * tsi * c1sq;
    let temp = d2 * tsi * c1 / 3.0;
    let d3 = (17.0 * aodp + s4) * temp;
    let d4 = 0.5 * temp * aodp * tsi * (221.0 * aodp + 31.0 * s4) * c1;
    let t3cof = d2 + 2.0 * c1sq;
    let t4cof = 0.25 * (3.0 * d3 + c1 * (12.0 * d2 + 10.0 * c1sq));
    let t5cof = 0.2 * (3.0 * d4 + 12.0 * c1 * d3 + 6.0 * d2 * d2 + 15.0 * c1sq * (2.0 * d2 + c1sq));
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
        let mmasmao =
            recover_original_mean_motion_and_semimajor_axis(XKE, XNO, TOTHRD, XINCL, EO, CK2);

        assert_eq!(mmasmao.xnodp, 0.07010615558630984);
        assert_eq!(mmasmao.aodp, 1.040117522759639);
        assert_eq!(mmasmao.betao2, 0.99992477733639);
        assert_eq!(mmasmao.betao, 0.9999623879608622);
        assert_eq!(mmasmao.x3thm1, -0.73895561738563);
        assert_eq!(mmasmao.theta2, 0.08701479420478998);
        assert_eq!(mmasmao.cosio, 0.29498270153483575);
    }
}
