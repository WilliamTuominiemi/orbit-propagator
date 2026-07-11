use eframe::egui;
use egui_plot::Line;
use egui_plot::Plot;
use egui_plot::PlotPoints;

mod constants;
mod types;

fn recover_original_mean_motion_and_semimajor_axis(
    xno: f64,
    xincl: f64,
    eo: f64,
) -> types::MeanMotionAndSemimajorAxisOutput {
    let a1 = (constants::XKE / xno).powf(constants::TOTHRD);
    let cosio = xincl.cos();
    let theta2 = cosio * cosio;
    let x3thm1 = 3.0 * theta2 - 1.0;
    let eosq = eo * eo;
    let betao2 = 1.0 - eosq;
    let betao = f64::sqrt(betao2);
    let del1 = 1.5 * constants::CK2 * x3thm1 / (a1 * a1 * betao * betao2);
    let ao = a1 * (1.0 - del1 * (0.5 * constants::TOTHRD + del1 * (1.0 + 134.0 / 81.0 * del1)));
    let delo = 1.5 * constants::CK2 * x3thm1 / (ao * ao * betao * betao2);
    let xnodp = xno / (1.0 + delo);
    let aodp = ao / (1.0 - delo);

    types::MeanMotionAndSemimajorAxisOutput {
        xnodp,
        aodp,
        betao2,
        betao,
        x3thm1,
        theta2,
        cosio,
    }
}

fn adjust_atmospheric_drag_for_low_orbit(aodp: f64, eo: f64) -> (f64, f64) {
    let mut s4 = constants::S;
    let mut qoms24 = constants::QOMS2T;
    let perigee = (aodp * (1.0 - eo) - constants::AE) * constants::XKMPER;
    if perigee < 156.0 {
        s4 = perigee - 78.0;
        if perigee <= 98.0 {
            s4 = 20.0;
        }
        qoms24 = ((120.0 - s4) * constants::AE / constants::XKMPER).powf(4.0);
        s4 = s4 / constants::XKMPER + constants::AE;
    }

    (s4, qoms24)
}

fn calculate_c_constants(
    eta: f64,
    coef: f64,
    xnodp: f64,
    aodp: f64,
    eeta: f64,
    tsi: f64,
    x3thm1: f64,
    bstar: f64,
    a3ovk2: f64,
    sinio: f64,
    eo: f64,
    betao2: f64,
    theta2: f64,
    omegao: f64,
) -> types::CConstants {
    let x1mth2 = 1.0 - theta2;
    let etasq = eta * eta;
    let psisq = (1.0 - etasq).abs();
    let coef1 = coef / psisq.powf(3.5);

    let c2 = coef1
        * xnodp
        * (aodp * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
            + 0.75 * constants::CK2 * tsi / psisq * x3thm1 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    let c1 = bstar * c2;
    let c3 = coef * tsi * a3ovk2 * xnodp * constants::AE * sinio / eo;
    let c4 = 2.0
        * xnodp
        * coef1
        * aodp
        * betao2
        * (eta * (2.0 + 0.5 * etasq) + eo * (0.5 + 2.0 * etasq)
            - 2.0 * constants::CK2 * tsi / (aodp * psisq)
                * (-3.0 * x3thm1 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta))
                    + 0.75 * x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * (2.0 * omegao).cos()));
    let c5 = 2.0 * coef1 * aodp * betao2 * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);

    types::CConstants { c1, c2, c3, c4, c5 }
}

fn calculate_d_constants(aodp: f64, tsi: f64, c1sq: f64, s4: f64, c1: f64) -> types::DConstants {
    let d2 = 4.0 * aodp * tsi * c1sq;
    let temp = d2 * tsi * c1 / 3.0;
    let d3 = (17.0 * aodp + s4) * temp;
    let d4 = 0.5 * temp * aodp * tsi * (221.0 * aodp + 31.0 * s4) * c1;

    types::DConstants { d2, d3, d4 }
}

fn update_for_secular_gravity_and_atmospheric_drag(
    tsince: f64,
    xmo: f64,
    xmdot: f64,
    omegao: f64,
    omgdot: f64,
    xnodeo: f64,
    xnodot: f64,
    xnodcf: f64,
    bstar: f64,
    omgcof: f64,
    xmcof: f64,
    eta: f64,
    aodp: f64,
    xnodp: f64,
    eo: f64,
    delmo: f64,
    sinmo: f64,
    t2cof: f64,
    t3cof: f64,
    t4cof: f64,
    t5cof: f64,
    c_constants: types::CConstants,
    d_constants: types::DConstants,
) -> types::SecularGravityAndAtmosphericDragUpdateOutput {
    let xmdf = xmo + xmdot * tsince;
    let omgadf = omegao + omgdot * tsince;
    let xnoddf = xnodeo + xnodot * tsince;
    let omega = omgadf;
    let xmp = xmdf;
    let tsq = tsince * tsince;
    let xnode = xnoddf + xnodcf * tsq;
    let mut tempa = 1.0 - c_constants.c1 * tsince;
    let mut tempe = bstar * c_constants.c4 * tsince;
    let mut templ = t2cof * tsq;
    let delomg = omgcof * tsince;
    let delm = xmcof * ((1.0 + eta * xmdf.cos()).powf(3.0) - delmo);
    let temp = delomg + delm;
    let xmp = xmdf + temp;
    let omega = omgadf - temp;
    let tcube = tsq * tsince;
    let tfour = tsince * tcube;
    tempa = tempa - d_constants.d2 * tsq - d_constants.d3 * tcube - d_constants.d4 * tfour;
    tempe = tempe + bstar * c_constants.c5 * (xmp.sin() - sinmo);
    templ = templ + t3cof * tcube + tfour * (t4cof + tsince * t5cof);
    let a = aodp * tempa.powf(2.0);
    let e = eo - tempe;
    let xl = xmp + omega + xnode + xnodp * templ;
    let beta = (1.0 - e * e).sqrt();
    let xn = constants::XKE / a.powf(1.5);

    types::SecularGravityAndAtmosphericDragUpdateOutput {
        e,
        a,
        xl,
        beta,
        xn,
        xnode,
        omega,
    }
}

fn long_period_periodics(
    sgaaduo: &types::SecularGravityAndAtmosphericDragUpdateOutput,
    xlcof: f64,
    aycof: f64,
) -> (f64, f64, f64) {
    let axn = sgaaduo.e * sgaaduo.omega.cos();
    let temp = 1.0 / (sgaaduo.a * sgaaduo.beta * sgaaduo.beta);
    let xll = temp * xlcof * axn;
    let aynl = temp * aycof;
    let xlt = sgaaduo.xl + xll;
    let ayn = sgaaduo.e * sgaaduo.omega.sin() + aynl;

    (xlt, ayn, axn)
}

fn short_periodics(
    r: f64,
    temp2: f64,
    betal: f64,
    x3thm1: f64,
    temp1: f64,
    x1mth2: f64,
    cos2u: f64,
    u: f64,
    x7thm1: f64,
    sin2u: f64,
    xnode: f64,
    cosio: f64,
    sinio: f64,
    xincl: f64,
    rdot: f64,
    xn: f64,
    rfdot: f64,
) -> types::ShortPeriodicsOutput {
    let rk = (r * (1.0 - 1.5 * temp2 * betal * x3thm1) + 0.5 * temp1 * x1mth2 * cos2u)
        * constants::XKMPER;
    let uk = u - 0.25 * temp2 * x7thm1 * sin2u;
    let xnodek = xnode + 1.5 * temp2 * cosio * sin2u;
    let xinck = xincl + 1.5 * temp2 * cosio * sinio * cos2u;
    let rdotk = (rdot - xn * temp1 * x1mth2 * sin2u) * constants::XKMPER / 60.0;
    let rfdotk = (rfdot + xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1)) * constants::XKMPER / 60.0;

    types::ShortPeriodicsOutput {
        rk,
        uk,
        xnodek,
        xinck,
        rdotk,
        rfdotk,
    }
}

fn fmod2p(x: f64) -> f64 {
    let rev = x / constants::TWOPI;
    let mut temp = x - (rev.trunc()) * constants::TWOPI;
    if temp < 0.0 {
        temp += constants::TWOPI;
    }
    temp
}

fn actan(sinx: f64, cosx: f64) -> f64 {
    let angle = sinx.atan2(cosx);
    if angle < 0.0 {
        angle + constants::TWOPI
    } else {
        angle
    }
}

fn keplers_equation(
    xlt: f64,
    xnode: f64,
    axn: f64,
    ayn: f64,
    e6a: f64,
) -> types::KeplersEquationOutput {
    let capu = fmod2p(xlt - xnode);
    let mut temp2 = capu;
    let mut temp3 = 0.0;
    let mut temp4 = 0.0;
    let mut temp5 = 0.0;
    let mut temp6 = 0.0;
    let mut epw = temp2;

    let mut sinepw = 0.0;
    let mut cosepw = 0.0;

    for _ in 0..10 {
        sinepw = temp2.sin();
        cosepw = temp2.cos();
        temp3 = axn * sinepw;
        temp4 = ayn * cosepw;
        temp5 = axn * cosepw;
        temp6 = ayn * sinepw;
        epw = (capu - temp4 + temp3 - temp2) / (1.0 - temp5 - temp6) + temp2;

        if (epw - temp2).abs() <= e6a {
            break;
        }

        temp2 = epw;
    }

    types::KeplersEquationOutput {
        temp2,
        temp3,
        temp4,
        temp5,
        temp6,
        sinepw,
        cosepw,
    }
}

fn short_period_prelimenary_quantities(
    keo: types::KeplersEquationOutput,
    axn: f64,
    ayn: f64,
    a: f64,
    temp1: f64,
) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    let ecose = keo.temp5 + keo.temp6;
    let esine = keo.temp3 - keo.temp4;
    let elsq = axn * axn + ayn * ayn;
    let temp = 1.0 - elsq;
    let pl = a * temp;
    let r = a * (1.0 - ecose);
    let temp1 = 1.0 / r;
    let rdot = constants::XKE * a.sqrt() * esine * temp1;
    let rfdot = constants::XKE * pl.sqrt() * temp1;
    let mut temp2 = a * temp1;
    let betal = temp.sqrt();
    let temp3 = 1.0 / (1.0 + betal);
    let cosu = temp2 * (keo.cosepw - axn + ayn * esine * temp3);
    let sinu = temp2 * (keo.sinepw - ayn - axn * esine * temp3);
    let u = actan(sinu, cosu);
    let sin2u = 2.0 * sinu * cosu;
    let cos2u = 2.0 * cosu * cosu - 1.0;
    let temp = 1.0 / pl;
    let temp1 = constants::CK2 * temp;
    temp2 = temp1 * temp;

    (r, rdot, rfdot, temp2, betal, temp1, cos2u, u, sin2u)
}

fn calculate_orientation_vectors(uk: f64, xinck: f64, xnodek: f64) -> types::OrientationVectors {
    let sinuk = uk.sin();
    let cosuk = uk.cos();
    let sinik = xinck.sin();
    let cosik = xinck.cos();
    let sinnok = xnodek.sin();
    let cosnok = xnodek.cos();
    let xmx = -sinnok * cosik;
    let xmy = cosnok * cosik;
    let ux = xmx * sinuk + cosnok * cosuk;
    let uy = xmy * sinuk + sinnok * cosuk;
    let uz = sinik * sinuk;
    let vx = xmx * cosuk - cosnok * sinuk;
    let vy = xmy * cosuk - sinnok * sinuk;
    let vz = sinik * cosuk;

    types::OrientationVectors {
        ux,
        uy,
        uz,
        vx,
        vy,
        vz,
    }
}

fn calculate_position_and_velocity(
    ov: types::OrientationVectors,
    spo: types::ShortPeriodicsOutput,
) -> types::PositionAndVelocity {
    let x = spo.rk * ov.ux;
    let y = spo.rk * ov.uy;
    let z = spo.rk * ov.uz;
    let xdot = spo.rdotk * ov.ux + spo.rfdotk * ov.vx;
    let ydot = spo.rdotk * ov.uy + spo.rfdotk * ov.vy;
    let zdot = spo.rdotk * ov.uz + spo.rfdotk * ov.vz;

    types::PositionAndVelocity {
        x,
        y,
        z,
        xdot,
        ydot,
        zdot,
    }
}

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
    let mmasmao = recover_original_mean_motion_and_semimajor_axis(xno, xincl, eo);

    let (xnodp, aodp, betao2, betao, x3thm1, theta2, cosio) = (
        mmasmao.xnodp,
        mmasmao.aodp,
        mmasmao.betao2,
        mmasmao.betao,
        mmasmao.x3thm1,
        mmasmao.theta2,
        mmasmao.cosio,
    );

    let (s4, qoms24) = adjust_atmospheric_drag_for_low_orbit(aodp, eo);

    let pinvsq = 1.0 / (aodp * aodp * betao2 * betao2);
    let tsi = 1.0 / (aodp - s4);
    let eta = aodp * eo * tsi;
    let eeta = eo * eta;

    let coef = qoms24 * tsi.powf(4.0);
    let sinio = xincl.sin();
    let a3ovk2 = -constants::XJ3 / constants::CK2 * constants::AE.powf(3.0);
    let x1mth2 = 1.0 - theta2;

    let c_constants = calculate_c_constants(
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

    let d_constants = calculate_d_constants(aodp, tsi, c1sq, s4, c_constants.c1);

    let t3cof = d_constants.d2 + 2.0 * c1sq;
    let t4cof =
        0.25 * (3.0 * d_constants.d3 + c_constants.c1 * (12.0 * d_constants.d2 + 10.0 * c1sq));
    let t5cof = 0.2
        * (3.0 * d_constants.d4
            + 12.0 * c_constants.c1 * d_constants.d3
            + 6.0 * d_constants.d2 * d_constants.d2
            + 15.0 * c1sq * (2.0 * d_constants.d2 + c1sq));

    let sgaaduo = update_for_secular_gravity_and_atmospheric_drag(
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

    let (xlt, ayn, axn) = long_period_periodics(&sgaaduo, xlcof, aycof);

    let keo = keplers_equation(xlt, sgaaduo.xnode, axn, ayn, e6a);

    let (r, rdot, rfdot, temp2, betal, temp1, cos2u, u, sin2u) =
        short_period_prelimenary_quantities(keo, axn, ayn, sgaaduo.a, temp1);

    let spo = short_periodics(
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

    let ov = calculate_orientation_vectors(spo.uk, spo.xinck, spo.xnodek);

    calculate_position_and_velocity(ov, spo)
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

    // Values are NORAD SPACETRACK REPORT NO. 3 SGP4 sample test case input parameters
    const XNO: f64 = 16.05824518 * (constants::TWOPI / constants::XMNPDA);
    const XINCL: f64 = 72.8435 * constants::DE2RA;
    const EO: f64 = 0.0086731;
    const BSTAR: f64 = 0.000066816;
    const OMEGAO: f64 = 52.6988 * constants::DE2RA;
    const XMO: f64 = 110.5714 * constants::DE2RA;
    const TSINCE: f64 = 0.0;
    const XNODEO: f64 = 115.9689 * constants::DE2RA;
    const E6A: f64 = 0.000001;

    // Values from SGP4 sample test case output values
    const POSITION_AND_VELOCITY_0: types::PositionAndVelocity = types::PositionAndVelocity {
        x: 2328.97048951,
        y: -5995.22076416,
        z: 1719.97067261,
        xdot: 2.91207230,
        ydot: -0.98341546,
        zdot: -7.09081703,
    };

    const TOLERANCE: f64 = 1e-3;

    #[test]
    fn test_recover_original_mean_motion_and_semimajor_axis() {
        let mmasmao = recover_original_mean_motion_and_semimajor_axis(XNO, XINCL, EO);

        assert_eq!(mmasmao.xnodp, 0.07010615558630984);
        assert_eq!(mmasmao.aodp, 1.040117522759639);
        assert_eq!(mmasmao.betao2, 0.99992477733639);
        assert_eq!(mmasmao.betao, 0.9999623879608622);
        assert_eq!(mmasmao.x3thm1, -0.73895561738563);
        assert_eq!(mmasmao.theta2, 0.08701479420478998);
        assert_eq!(mmasmao.cosio, 0.29498270153483575);
    }

    #[test]
    fn test_adjust_atmospheric_drag_for_low_orbit() {
        let aodp = 1.040117522759639; // From previous test, same sample test case
        let (s4, qoms24) = adjust_atmospheric_drag_for_low_orbit(aodp, EO);

        assert_eq!(s4, constants::S);
        assert_eq!(qoms24, constants::QOMS2T);
    }

    #[test]
    fn test_calculate_c_constants() {
        // From recover_original_mean_motion_and_semimajor_axis function
        let xnodp = 0.07010615558630984;
        let aodp = 1.040117522759639;
        let x3thm1 = -0.73895561738563;

        let eta = 0.3234711976798404;
        let coef = 0.003108405951369967;
        let eeta = 0.0028054980445970236;
        let tsi = 35.85740444884659;
        let a3ovk2 = 0.004690139440023056;
        let sinio = 0.9555025932959105;
        let betao2 = 0.99992477733639;
        let theta2 = 0.08701479420479;

        let c_constants = calculate_c_constants(
            eta, coef, xnodp, aodp, eeta, tsi, x3thm1, BSTAR, a3ovk2, sinio, EO, betao2, theta2,
            OMEGAO,
        );

        assert_eq!(c_constants.c1, 2.3338044215116538e-8);
        assert_eq!(c_constants.c2, 0.0003492882575298811);
        assert_eq!(c_constants.c3, 0.004037532255765166);
        assert_eq!(c_constants.c4, 0.000377201121554739);
        assert_eq!(c_constants.c5, 0.012334919304344908);
    }

    #[test]
    fn test_calculate_d_constants() {
        let aodp = 1.040117522759639;
        let tsi = 35.85740444884659;
        let c1sq = 0.0000000000000005446643077867345;
        let s4 = 1.01222928;
        let c1 = 2.3338044215116538e-8;

        let d_constants = calculate_d_constants(aodp, tsi, c1sq, s4, c1);

        assert_eq!(d_constants.d2, 8.12550142270866e-14);
        assert_eq!(d_constants.d3, 4.2372075736327043e-19);
        assert_eq!(d_constants.d4, 2.5770097992217537e-24);
    }

    #[test]
    fn test_update_for_secular_gravity_and_atmospheric_drag() {
        let xnodp = 0.07010615558630984;
        let aodp = 1.040117522759639;
        let eta = 0.3234711976798404;
        let xmdot = 0.07006729335201786;
        let omgdot = -0.00002971792465285666;
        let xnodot = -0.00003096311254169127;
        let xnodcf = -0.000000000002535821899421168;
        let omgcof = 0.00000016348304905484922;
        let xmcof = -0.000049353388663657485;
        let delmo = 0.6963086765241224;
        let sinmo = 0.9362350466329594;
        let t2cof = 0.00000003500706632267481;
        let t3cof = 0.00000000000008234434284266006;
        let t4cof = 0.00000000000000000032351134586589164;
        let t5cof = 0.0000000000000000000000015781283156000947;

        let c_constants = types::CConstants {
            c1: 2.3338044215116538e-8,
            c2: 0.0003492882575298811,
            c3: 0.004037532255765166,
            c4: 0.000377201121554739,
            c5: 0.012334919304344908,
        };
        let d_constants = types::DConstants {
            d2: 8.12550142270866e-14,
            d3: 4.2372075736327043e-19,
            d4: 2.5770097992217537e-24,
        };

        let sgaaduo = update_for_secular_gravity_and_atmospheric_drag(
            TSINCE,
            XMO,
            xmdot,
            OMEGAO,
            omgdot,
            XNODEO,
            xnodot,
            xnodcf,
            BSTAR,
            omgcof,
            xmcof,
            eta,
            aodp,
            xnodp,
            EO,
            delmo,
            sinmo,
            t2cof,
            t3cof,
            t4cof,
            t5cof,
            c_constants,
            d_constants,
        );

        assert_eq!(sgaaduo.e, 0.0086731);
        assert_eq!(sgaaduo.a, 1.040117522759639);
        assert_eq!(sgaaduo.xl, 4.873641689736749);
        assert_eq!(sgaaduo.beta, 0.9999623879608622);
        assert_eq!(sgaaduo.xn, 0.07010615556528188);
    }

    #[test]
    fn test_long_period_periodics() {
        let sgaaduo = types::SecularGravityAndAtmosphericDragUpdateOutput {
            e: 0.0086731,
            a: 1.040117522759639,
            xl: 0.07010615558630984,
            beta: 0.9999623879608622,
            xn: 0.07010615556528188,
            xnode: 0.0,
            omega: 1.6161106125158646,
        };

        let xlcof = 0.001935745758076399;
        let aycof = 0.0011203600994653647;

        let (xlt, ayn, axn) = long_period_periodics(&sgaaduo, xlcof, aycof);

        assert_eq!(xlt, 0.07010542434717408);
        assert_eq!(ayn, 0.009741425556977382);
        assert_eq!(axn, -0.0003928808433640341);
    }

    #[test]
    fn test_keplers_equation() {
        let xlt = 0.07011593807103583;
        let xnode = 0.0;
        let axn = 0.005255942497390392;
        let ayn = 0.007976339600468509;

        let keplers_equation_output = keplers_equation(xlt, xnode, axn, ayn, E6A);

        assert_eq!(keplers_equation_output.temp2, 0.06248313642895434);
        assert_eq!(keplers_equation_output.temp3, 0.0003281941220562471);
        assert_eq!(keplers_equation_output.temp4, 0.007960774282990234);
        assert_eq!(keplers_equation_output.temp5, 0.0052456858611741215);
        assert_eq!(keplers_equation_output.temp6, 0.0004980624833125527);
    }

    #[test]
    fn test_short_periodics() {
        let (
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
            xnode,
            cosio,
            sinio,
            xincl,
            rdot,
            xn,
            rfdot,
        ) = (
            1.0430516048741472,
            0.0005004479094222182,
            0.9999543754967234,
            -0.73895561738563,
            0.0005204771435457296,
            0.91298520579521,
            0.8538389429983284,
            2.867852612718983,
            -0.3908964405664701,
            -0.5205372795674639,
            2.02403913260325,
            0.29498270153483575,
            0.9555025932959105,
            1.27135891222375,
            0.0006636054871099481,
            0.07010615556528188,
            0.07271020474319065,
        );

        assert_eq!(
            types::ShortPeriodicsOutput {
                rk: 6657.7080462027025,
                uk: 2.867827155413039,
                xnodek: 2.0239238673191204,
                xinck: 1.2715395691080913,
                rdotk: 0.07238614054582597,
                rfdotk: 7.727982650878059
            },
            short_periodics(
                r, temp2, betal, x3thm1, temp1, x1mth2, cos2u, u, x7thm1, sin2u, xnode, cosio,
                sinio, xincl, rdot, xn, rfdot
            )
        )
    }

    #[test]
    fn test_calculate_orientation_vectors() {
        let uk = 2.867827155413039;
        let xinck = 1.2715395691080913;
        let xnodek = 2.0239238673191204;

        assert_eq!(
            types::OrientationVectors {
                ux: 0.3498156817129497,
                uy: -0.9004932019514819,
                uz: 0.25834276081762364,
                vx: 0.3735451525121787,
                vy: -0.11881912112969674,
                vz: -0.9199706709937114,
            },
            calculate_orientation_vectors(uk, xinck, xnodek)
        )
    }

    #[test]
    fn test_calculate_position_and_velocity() {
        let ov = types::OrientationVectors {
            ux: 0.3498156817129497,
            uy: -0.9004932019514819,
            uz: 0.25834276081762364,
            vx: 0.3735451525121787,
            vy: -0.11881912112969674,
            vz: -0.9199706709937114,
        };
        let spo = types::ShortPeriodicsOutput {
            rk: 6657.7080462027025,
            uk: 2.867827155413039,
            xnodek: 2.0239238673191204,
            xinck: 1.2715395691080913,
            rdotk: 0.07238614054582597,
            rfdotk: 7.727982650878059,
        };
        let pav = calculate_position_and_velocity(ov, spo);

        assert_abs_diff_eq!(POSITION_AND_VELOCITY_0.x, pav.x, epsilon = TOLERANCE);
        assert_abs_diff_eq!(POSITION_AND_VELOCITY_0.y, pav.y, epsilon = TOLERANCE);
        assert_abs_diff_eq!(POSITION_AND_VELOCITY_0.z, pav.z, epsilon = TOLERANCE);
        assert_abs_diff_eq!(POSITION_AND_VELOCITY_0.xdot, pav.xdot, epsilon = TOLERANCE);
        assert_abs_diff_eq!(POSITION_AND_VELOCITY_0.ydot, pav.ydot, epsilon = TOLERANCE);
        assert_abs_diff_eq!(POSITION_AND_VELOCITY_0.zdot, pav.zdot, epsilon = TOLERANCE);
    }

    #[test]
    fn test_sgp4() {
        let ouput = sgp4(
            TSINCE,
            EO,
            BSTAR,
            XINCL,
            OMEGAO,
            constants::CK4,
            XMO,
            XNODEO,
            E6A,
            XNO,
        );

        assert_abs_diff_eq!(POSITION_AND_VELOCITY_0.x, ouput.x, epsilon = TOLERANCE);
        assert_abs_diff_eq!(POSITION_AND_VELOCITY_0.y, ouput.y, epsilon = TOLERANCE);
        assert_abs_diff_eq!(POSITION_AND_VELOCITY_0.z, ouput.z, epsilon = TOLERANCE);
        assert_abs_diff_eq!(
            POSITION_AND_VELOCITY_0.xdot,
            ouput.xdot,
            epsilon = TOLERANCE
        );
        assert_abs_diff_eq!(
            POSITION_AND_VELOCITY_0.ydot,
            ouput.ydot,
            epsilon = TOLERANCE
        );
        assert_abs_diff_eq!(
            POSITION_AND_VELOCITY_0.zdot,
            ouput.zdot,
            epsilon = TOLERANCE
        );
    }
}
