use eframe::egui;
use egui_plot::Line;
use egui_plot::Plot;
use egui_plot::PlotPoints;

mod constants;
mod helpers;
mod sgp4;
mod test_constants;
mod types;

fn sgp4(
    tsince: f64,
    eo: f64,
    bstar: f64,
    xincl: f64,
    omegao: f64,
    ck4: f64,
    xmo: f64,
    xnodeo: f64,
    e6a: f64,
    xno: f64,
) -> types::PositionAndVelocity {
    let mmasmao = helpers::recover_original_mean_motion_and_semimajor_axis(xno, xincl, eo);

    let (xnodp, aodp, betao2, betao, x3thm1, theta2, cosio) = (
        mmasmao.xnodp,
        mmasmao.aodp,
        mmasmao.betao2,
        mmasmao.betao,
        mmasmao.x3thm1,
        mmasmao.theta2,
        mmasmao.cosio,
    );

    let (s4, qoms24) = helpers::adjust_atmospheric_drag_for_low_orbit(aodp, eo);

    let pinvsq = 1.0 / (aodp * aodp * betao2 * betao2);
    let tsi = 1.0 / (aodp - s4);
    let eta = aodp * eo * tsi;
    let eeta = eo * eta;

    let coef = qoms24 * tsi.powf(4.0);
    let sinio = xincl.sin();
    let a3ovk2 = -constants::XJ3 / constants::CK2 * constants::AE.powf(3.0);
    let x1mth2 = 1.0 - theta2;

    let c_constants = helpers::calculate_c_constants(
        eta, coef, xnodp, aodp, eeta, tsi, x3thm1, bstar, a3ovk2, sinio, eo, betao2, theta2, omegao,
    );

    let theta4 = theta2 * theta2;
    let temp1 = 3.0 * constants::CK2 * pinvsq * xnodp;
    let temp2 = temp1 * constants::CK2 * pinvsq;
    let temp3 = 1.25 * ck4 * pinvsq * pinvsq * xnodp;

    let xmdot = xnodp
        + 0.5 * temp1 * betao * x3thm1
        + 0.0625 * temp2 * betao * (13.0 - 78.0 * theta2 + 137.0 * theta4);
    let x1m5th = 1.0 - 5.0 * theta2;
    let omgdot = -0.5 * temp1 * x1m5th
        + 0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4)
        + temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);
    let xhdot1 = -temp1 * cosio;
    let xnodot =
        xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 * (3.0 - 7.0 * theta2)) * cosio;
    let omgcof = bstar * c_constants.c3 * omegao.cos();
    let xmcof = -constants::TOTHRD * coef * bstar * constants::AE / eeta;
    let xnodcf = 3.5 * betao2 * xhdot1 * c_constants.c1;
    let t2cof = 1.5 * c_constants.c1;
    let xlcof = 0.125 * a3ovk2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
    let aycof = 0.25 * a3ovk2 * sinio;
    let delmo = (1.0 + eta * xmo.cos()).powf(3.0);
    let sinmo = xmo.sin();
    let x7thm1 = 7.0 * theta2 - 1.0;
    let c1sq = c_constants.c1 * c_constants.c1;

    let d_constants = helpers::calculate_d_constants(aodp, tsi, c1sq, s4, c_constants.c1);

    let t3cof = d_constants.d2 + 2.0 * c1sq;
    let t4cof =
        0.25 * (3.0 * d_constants.d3 + c_constants.c1 * (12.0 * d_constants.d2 + 10.0 * c1sq));
    let t5cof = 0.2
        * (3.0 * d_constants.d4
            + 12.0 * c_constants.c1 * d_constants.d3
            + 6.0 * d_constants.d2 * d_constants.d2
            + 15.0 * c1sq * (2.0 * d_constants.d2 + c1sq));

    let sgaaduo = helpers::update_for_secular_gravity_and_atmospheric_drag(
        tsince,
        xmo,
        xmdot,
        omegao,
        omgdot,
        xnodeo,
        xnodot,
        xnodcf,
        bstar,
        omgcof,
        xmcof,
        eta,
        aodp,
        xnodp,
        eo,
        delmo,
        sinmo,
        t2cof,
        t3cof,
        t4cof,
        t5cof,
        c_constants,
        d_constants,
    );

    let (xlt, ayn, axn) = helpers::long_period_periodics(&sgaaduo, xlcof, aycof);

    let keo = helpers::keplers_equation(xlt, sgaaduo.xnode, axn, ayn, e6a);

    let (r, rdot, rfdot, temp2, betal, temp1, cos2u, u, sin2u) =
        helpers::short_period_prelimenary_quantities(keo, axn, ayn, sgaaduo.a, temp1);

    let spo = helpers::short_periodics(
        r,
        temp2,
        betal,
        x3thm1,
        temp1,
        x1mth2,
        cos2u,
        u,
        x7thm1,
        sin2u,
        sgaaduo.xnode,
        cosio,
        sinio,
        xincl,
        rdot,
        sgaaduo.xn,
        rfdot,
    );

    let ov = helpers::calculate_orientation_vectors(spo.uk, spo.xinck, spo.xnodek);

    helpers::calculate_position_and_velocity(ov, spo)
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
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_sgp4() {
        let ouput = sgp4(
            test_constants::TSINCE,
            test_constants::EO,
            test_constants::BSTAR,
            test_constants::XINCL,
            test_constants::OMEGAO,
            constants::CK4,
            test_constants::XMO,
            test_constants::XNODEO,
            test_constants::E6A,
            test_constants::XNO,
        );

        assert_abs_diff_eq!(
            test_constants::POSITION_AND_VELOCITY_0.x,
            ouput.x,
            epsilon = test_constants::TOLERANCE
        );
        assert_abs_diff_eq!(
            test_constants::POSITION_AND_VELOCITY_0.y,
            ouput.y,
            epsilon = test_constants::TOLERANCE
        );
        assert_abs_diff_eq!(
            test_constants::POSITION_AND_VELOCITY_0.z,
            ouput.z,
            epsilon = test_constants::TOLERANCE
        );
        assert_abs_diff_eq!(
            test_constants::POSITION_AND_VELOCITY_0.xdot,
            ouput.xdot,
            epsilon = test_constants::TOLERANCE
        );
        assert_abs_diff_eq!(
            test_constants::POSITION_AND_VELOCITY_0.ydot,
            ouput.ydot,
            epsilon = test_constants::TOLERANCE
        );
        assert_abs_diff_eq!(
            test_constants::POSITION_AND_VELOCITY_0.zdot,
            ouput.zdot,
            epsilon = test_constants::TOLERANCE
        );
    }
}
