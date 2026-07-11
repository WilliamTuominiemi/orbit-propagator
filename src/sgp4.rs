use crate::constants;
use crate::helpers;
use crate::test_constants;
use crate::types;

pub struct Sgp4 {
    xmdot: f64,
    omgdot: f64,
    omgcof: f64,
    xmcof: f64,
    xnodcf: f64,
    xlcof: f64,
    aycof: f64,
    delmo: f64,
    sinmo: f64,
    x7thm1: f64,
    c_constants: types::CConstants,
    d_constants: types::DConstants,
    t2cof: f64,
    t3cof: f64,
    t4cof: f64,
    t5cof: f64,
}

impl Sgp4 {
    pub fn new(eo: f64, bstar: f64, xincl: f64, omegao: f64, ck4: f64, xmo: f64, xno: f64) -> Self {
        let mmasmao = helpers::recover_original_mean_motion_and_semimajor_axis(xno, xincl, eo);
        let (s4, qoms24) = helpers::adjust_atmospheric_drag_for_low_orbit(mmasmao.aodp, eo);

        let pinvsq = 1.0 / (mmasmao.aodp * mmasmao.aodp * mmasmao.betao2 * mmasmao.betao2);
        let tsi = 1.0 / (mmasmao.aodp - s4);
        let eta = mmasmao.aodp * eo * tsi;
        let eeta = eo * eta;

        let coef = qoms24 * tsi.powf(4.0);
        let sinio = xincl.sin();
        let a3ovk2 = -constants::XJ3 / constants::CK2 * constants::AE.powf(3.0);
        let x1mth2 = 1.0 - mmasmao.theta2;

        let c_constants = helpers::calculate_c_constants(
            eta,
            coef,
            mmasmao.xnodp,
            mmasmao.aodp,
            eeta,
            tsi,
            mmasmao.x3thm1,
            bstar,
            a3ovk2,
            sinio,
            eo,
            mmasmao.betao2,
            mmasmao.theta2,
            omegao,
        );

        let theta4 = mmasmao.theta2 * mmasmao.theta2;
        let temp1 = 3.0 * constants::CK2 * pinvsq * mmasmao.xnodp;
        let temp2 = temp1 * constants::CK2 * pinvsq;
        let temp3 = 1.25 * ck4 * pinvsq * pinvsq * mmasmao.xnodp;

        let xmdot = mmasmao.xnodp
            + 0.5 * temp1 * mmasmao.betao * mmasmao.x3thm1
            + 0.0625 * temp2 * mmasmao.betao * (13.0 - 78.0 * mmasmao.theta2 + 137.0 * theta4);
        let x1m5th = 1.0 - 5.0 * mmasmao.theta2;
        let omgdot = -0.5 * temp1 * x1m5th
            + 0.0625 * temp2 * (7.0 - 114.0 * mmasmao.theta2 + 395.0 * theta4)
            + temp3 * (3.0 - 36.0 * mmasmao.theta2 + 49.0 * theta4);
        let xhdot1 = -temp1 * mmasmao.cosio;
        let xnodot = xhdot1
            + (0.5 * temp2 * (4.0 - 19.0 * mmasmao.theta2)
                + 2.0 * temp3 * (3.0 - 7.0 * mmasmao.theta2))
                * mmasmao.cosio;
        let omgcof = bstar * c_constants.c3 * omegao.cos();
        let xmcof = -constants::TOTHRD * coef * bstar * constants::AE / eeta;
        let xnodcf = 3.5 * mmasmao.betao2 * xhdot1 * c_constants.c1;
        let t2cof = 1.5 * c_constants.c1;
        let xlcof = 0.125 * a3ovk2 * sinio * (3.0 + 5.0 * mmasmao.cosio) / (1.0 + mmasmao.cosio);
        let aycof = 0.25 * a3ovk2 * sinio;
        let delmo = (1.0 + eta * xmo.cos()).powf(3.0);
        let sinmo = xmo.sin();
        let x7thm1 = 7.0 * mmasmao.theta2 - 1.0;
        let c1sq = c_constants.c1 * c_constants.c1;

        let d_constants =
            helpers::calculate_d_constants(mmasmao.aodp, tsi, c1sq, s4, c_constants.c1);

        let t3cof = d_constants.d2 + 2.0 * c1sq;
        let t4cof =
            0.25 * (3.0 * d_constants.d3 + c_constants.c1 * (12.0 * d_constants.d2 + 10.0 * c1sq));
        let t5cof = 0.2
            * (3.0 * d_constants.d4
                + 12.0 * c_constants.c1 * d_constants.d3
                + 6.0 * d_constants.d2 * d_constants.d2
                + 15.0 * c1sq * (2.0 * d_constants.d2 + c1sq));

        Sgp4 {
            xmdot,
            omgdot,
            omgcof,
            xmcof,
            xnodcf,
            xlcof,
            aycof,
            delmo,
            sinmo,
            x7thm1,
            c_constants,
            d_constants,
            t2cof,
            t3cof,
            t4cof,
            t5cof,
        }
    }

    pub fn propagate(&self, tsince: f64) -> types::PositionAndVelocity {
        types::PositionAndVelocity {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            xdot: 0.0,
            ydot: 0.0,
            zdot: 0.0,
        }
    }
}
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let sgp4 = Sgp4::new(
            test_constants::EO,
            test_constants::BSTAR,
            test_constants::XINCL,
            test_constants::OMEGAO,
            constants::CK4,
            test_constants::XMO,
            test_constants::XNO,
        );

        assert_eq!(sgp4.xmdot, 0.07006729335201786);
        assert_eq!(sgp4.omgdot, -0.00002971792465285666);
        assert_eq!(sgp4.omgcof, 0.00000016348304905484922);
        assert_eq!(sgp4.xmcof, -0.000049353388663657485);
        assert_eq!(sgp4.xnodcf, -0.000000000002535821899421168);
        assert_eq!(sgp4.xlcof, 0.001935745758076399);
        assert_eq!(sgp4.aycof, 0.0011203600994653647);
        assert_eq!(sgp4.delmo, 0.6963086765241224);
        assert_eq!(sgp4.sinmo, 0.9362350466329594);
        assert_eq!(sgp4.x7thm1, -0.3908964405664701);
        assert_eq!(sgp4.t2cof, 0.00000003500706632267481);
        assert_eq!(sgp4.t3cof, 0.00000000000008234434284266006);
        assert_eq!(sgp4.t4cof, 0.00000000000000000032351134586589164);
        assert_eq!(sgp4.t5cof, 0.0000000000000000000000015781283156000947);

        assert_eq!(sgp4.c_constants.c1, 2.3338044215116538e-8);
        assert_eq!(sgp4.c_constants.c2, 0.0003492882575298811);
        assert_eq!(sgp4.c_constants.c3, 0.004037532255765166);
        assert_eq!(sgp4.c_constants.c4, 0.000377201121554739);
        assert_eq!(sgp4.c_constants.c5, 0.012334919304344908);

        assert_eq!(sgp4.d_constants.d2, 8.12550142270866e-14);
        assert_eq!(sgp4.d_constants.d3, 4.2372075736327043e-19);
        assert_eq!(sgp4.d_constants.d4, 2.5770097992217537e-24);
    }
}
