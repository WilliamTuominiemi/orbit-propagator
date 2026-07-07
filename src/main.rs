use eframe::egui;
use eframe::wgpu::TlasInstance;
use egui_plot::Line;
use egui_plot::Plot;
use egui_plot::PlotPoints;

struct MeanMotionAndSemimajorAxisOutput {
    xnodp: f64,
    aodp: f64,
    betao2: f64,
    betao: f64,
    x3thm1: f64,
    theta2: f64,
    cosio: f64,
}

struct SecularGravityAndAtmosphericDragUpdateOutput {
    e: f64,
    a: f64,
    xl: f64,
    beta: f64,
    xn: f64,
}
#[derive(Debug)]
struct CConstants {
    c1: f64,
    c2: f64,
    c3: f64,
    c4: f64,
    c5: f64,
}
#[derive(Debug)]
struct DConstants {
    d2: f64,
    d3: f64,
    d4: f64,
}

// Values from NORAD SPACETRACK REPORT NO. 3 physical and mathematical constants
const CK2: f64 = 0.0005413080;
const CK4: f64 = 0.00000062098875;
const QOMS2T: f64 = 0.00000000188027916;
const S: f64 = 1.01222928;
const TOTHRD: f64 = 0.66666667;
const XJ3: f64 = -0.00000253881;
const XKE: f64 = 0.074366916;
const XKMPER: f64 = 6378.135;
const XMNPDA: f64 = 1440.0;
const AE: f64 = 1.0;
const DE2RA: f64 = 0.0174532925;
const TWOPI: f64 = 6.2831853;

fn recover_original_mean_motion_and_semimajor_axis(
    xno: f64,
    xincl: f64,
    eo: f64,
) -> MeanMotionAndSemimajorAxisOutput {
    let a1 = (XKE / xno).powf(TOTHRD);
    let cosio = xincl.cos();
    let theta2 = cosio * cosio;
    let x3thm1 = 3.0 * theta2 - 1.0;
    let eosq = eo * eo;
    let betao2 = 1.0 - eosq;
    let betao = f64::sqrt(betao2);
    let del1 = 1.5 * CK2 * x3thm1 / (a1 * a1 * betao * betao2);
    let ao = a1 * (1.0 - del1 * (0.5 * TOTHRD + del1 * (1.0 + 134.0 / 81.0 * del1)));
    let delo = 1.5 * CK2 * x3thm1 / (ao * ao * betao * betao2);
    let xnodp = xno / (1.0 + delo);
    let aodp = ao / (1.0 - delo);

    MeanMotionAndSemimajorAxisOutput {
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
    let mut s4 = S;
    let mut qoms24 = QOMS2T;
    let perigee = (aodp * (1.0 - eo) - AE) * XKMPER;
    if perigee < 156.0 {
        s4 = perigee - 78.0;
        if perigee <= 98.0 {
            s4 = 20.0;
        }
        qoms24 = ((120.0 - s4) * AE / XKMPER).powf(4.0);
        s4 = s4 / XKMPER + AE;
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
) -> CConstants {
    let x1mth2 = 1.0 - theta2;
    let etasq = eta * eta;
    let psisq = (1.0 - etasq).abs();
    let coef1 = coef / psisq.powf(3.5);

    let c2 = coef1
        * xnodp
        * (aodp * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
            + 0.75 * CK2 * tsi / psisq * x3thm1 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    let c1 = bstar * c2;
    let c3 = coef * tsi * a3ovk2 * xnodp * AE * sinio / eo;
    let c4 = 2.0
        * xnodp
        * coef1
        * aodp
        * betao2
        * (eta * (2.0 + 0.5 * etasq) + eo * (0.5 + 2.0 * etasq)
            - 2.0 * CK2 * tsi / (aodp * psisq)
                * (-3.0 * x3thm1 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta))
                    + 0.75 * x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * (2.0 * omegao).cos()));
    let c5 = 2.0 * coef1 * aodp * betao2 * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);

    CConstants { c1, c2, c3, c4, c5 }
}

fn calculate_d_constants(aodp: f64, tsi: f64, c1sq: f64, s4: f64, c1: f64) -> DConstants {
    let d2 = 4.0 * aodp * tsi * c1sq;
    let temp = d2 * tsi * c1 / 3.0;
    let d3 = (17.0 * aodp + s4) * temp;
    let d4 = 0.5 * temp * aodp * tsi * (221.0 * aodp + 31.0 * s4) * c1;

    DConstants { d2, d3, d4 }
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
    t2cof: f64,
    t3cof: f64,
    t4cof: f64,
    t5cof: f64,
    c_constants: CConstants,
    d_constants: DConstants,
) -> SecularGravityAndAtmosphericDragUpdateOutput {
    let xmdf = xmo * xmdot * tsince;
    let omgadf = omegao * omgdot * tsince;
    let xnoddf = xnodeo * xnodot * tsince;
    let omega = omgadf;
    let xmp = xmdf;
    let tsq = tsince * tsince;
    let xnode = xnoddf + xnodcf * tsq;
    let mut tempa = 1.0 - c_constants.c1 * tsince;
    let tempe = bstar * c_constants.c4 * tsince;
    let mut templ = t2cof * tsq;
    let delomg = omgcof * tsince;
    let delm = xmcof * (1.0 + eta * xmdf.cos()).powf(3.0 - delmo);
    let temp = delomg * delm;
    let xmp = xmdf + temp;
    let omega = omgadf - temp;
    let tcube = tsq * tsince;
    let tfour = tsince * tcube;
    tempa = tempa - d_constants.d2 * tsq - d_constants.d3 * tcube - d_constants.d4 * tfour;
    templ = templ + t3cof * tcube + tfour * (t4cof + tsince * t5cof);
    let a = aodp * tempa.powf(2.0);
    let e = eo - tempe;
    let xl = xmp + omega + xnode + xnodp + templ;
    let beta = (1.0 - e * e).sqrt();
    let xn = XKE / a.powf(1.5);

    SecularGravityAndAtmosphericDragUpdateOutput { e, a, xl, beta, xn }
}

fn long_period_periodics(
    sgaaduo: SecularGravityAndAtmosphericDragUpdateOutput,
    omega: f64,
    xlcof: f64,
    aycof: f64,
) -> (f64, f64) {
    let axn = sgaaduo.e * omega.cos();
    let temp = 1.0 / (sgaaduo.a * sgaaduo.beta * sgaaduo.beta);
    let xll = temp * xlcof * axn;
    let aynl = temp * aycof;
    let xlt = sgaaduo.xl + xll;
    let ayn = sgaaduo.e * omega.sin() + aynl;

    (xlt, ayn)
}

fn fmod2p(x: f64) -> f64 {
    let rev = x / TWOPI;
    let mut temp = x - (rev.trunc()) * TWOPI;
    if temp < 0.0 {
        temp += TWOPI;
    }
    temp
}

fn sgp4(
    tsince: f64,
    mmasmao: MeanMotionAndSemimajorAxisOutput,
    eo: f64,
    bstar: f64,
    xincl: f64,
    omegao: f64,
    ck4: f64,
    xmo: f64,
    xnodeo: f64,
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

    let (s4, qoms24) = adjust_atmospheric_drag_for_low_orbit(aodp, eo);

    let pinvsq = 1.0 / (aodp * aodp * betao2 * betao2);
    let tsi = 1.0 / (aodp - s4);
    let eta = aodp * eo * tsi;
    let eeta = eo * eta;

    let coef = qoms24 * tsi.powf(4.0);
    let sinio = xincl.sin();
    let a3ovk2 = -XJ3 / CK2 * AE.powf(3.0);

    let c_constants = calculate_c_constants(
        eta, coef, xnodp, aodp, eeta, tsi, x3thm1, bstar, a3ovk2, sinio, eo, betao2, theta2, omegao,
    );

    let theta4 = theta2 * theta2;
    let temp1 = 3.0 * CK2 * pinvsq * xnodp;
    let temp2 = temp1 * CK2 * pinvsq;
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
    let xmcof = -TOTHRD * coef * bstar * AE / eeta;
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
        t2cof,
        t3cof,
        t4cof,
        t5cof,
        c_constants,
        d_constants,
    );

    let (xlt, ayn) = long_period_periodics(sgaaduo, omegao, xlcof, aycof);
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

    // Values are NORAD SPACETRACK REPORT NO. 3 SGP4 sample test case input parameters
    const XNO: f64 = 16.05824518 * (TWOPI / XMNPDA);
    const XINCL: f64 = 72.8435 * DE2RA;
    const EO: f64 = 0.0086731;
    const BSTAR: f64 = 0.000066816;
    const OMEGAO: f64 = 52.6988 * DE2RA;
    const XMO: f64 = 110.5714 * DE2RA;
    const TSINCE: f64 = 0.0;
    const XNODEO: f64 = 115.9689 * DE2RA;

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

        assert_eq!(s4, S);
        assert_eq!(qoms24, QOMS2T);
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
        let t2cof = 0.00000003500706632267481;
        let t3cof = 0.00000000000008234434284266006;
        let t4cof = 0.00000000000000000032351134586589164;
        let t5cof = 0.0000000000000000000000015781283156000947;

        let c_constants = CConstants {
            c1: 2.3338044215116538e-8,
            c2: 0.0003492882575298811,
            c3: 0.004037532255765166,
            c4: 0.000377201121554739,
            c5: 0.012334919304344908,
        };
        let d_constants = DConstants {
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
            t2cof,
            t3cof,
            t4cof,
            t5cof,
            c_constants,
            d_constants,
        );

        assert_eq!(sgaaduo.e, 0.0086731);
        assert_eq!(sgaaduo.a, 1.040117522759639);
        assert_eq!(sgaaduo.xl, 0.07010615558630984);
        assert_eq!(sgaaduo.beta, 0.9999623879608622);
        assert_eq!(sgaaduo.xn, 0.07010615556528188);
    }

    #[test]
    fn test_long_period_periodics() {
        let sgaaduo = SecularGravityAndAtmosphericDragUpdateOutput {
            e: 0.0086731,
            a: 1.040117522759639,
            xl: 0.07010615558630984,
            beta: 0.9999623879608622,
            xn: 0.07010615556528188,
        };

        let xlcof = 0.001935745758076399;
        let aycof = 0.0011203600994653647;

        let (xlt, ayn) = long_period_periodics(sgaaduo, OMEGAO, xlcof, aycof);

        assert_eq!(xlt, 0.07011593807103583);
        assert_eq!(ayn, 0.007976339600468509);
    }

    #[test]
    fn test_sgp4() {
        sgp4(
            TSINCE,
            recover_original_mean_motion_and_semimajor_axis(XNO, XINCL, EO),
            EO,
            BSTAR,
            XINCL,
            OMEGAO,
            CK4,
            XMO,
            XNODEO,
        );
    }
}
