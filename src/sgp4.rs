use crate::constants;
use crate::helpers;
use crate::test_constants;
use crate::types;

pub struct Sgp4 {
    eo: f64,
    bstar: f64,
    xincl: f64,
    omegao: f64,
    xmo: f64,
    xnodeo: f64,
    e6a: f64,
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
    x1mth2: f64,
    xnodot: f64,
    eta: f64,
    sinio: f64,
    c_constants: types::CConstants,
    d_constants: types::DConstants,
    t2cof: f64,
    t3cof: f64,
    t4cof: f64,
    t5cof: f64,
    mmasmao: types::MeanMotionAndSemimajorAxisOutput,
}

impl Sgp4 {
    pub fn new(
        eo: f64,
        bstar: f64,
        xincl: f64,
        omegao: f64,
        xmo: f64,
        xno: f64,
        xnodeo: f64,
        e6a: f64,
    ) -> Self {
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
        let temp3 = 1.25 * constants::CK4 * pinvsq * pinvsq * mmasmao.xnodp;

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
            eo,
            bstar,
            xincl,
            omegao,
            xmo,
            xnodeo,
            e6a,
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
            x1mth2,
            xnodot,
            eta,
            sinio,
            c_constants,
            d_constants,
            t2cof,
            t3cof,
            t4cof,
            t5cof,
            mmasmao,
        }
    }

    pub fn propagate(&self, tsince: f64) -> types::PositionAndVelocity {
        let sgaaduo = self.update_for_secular_gravity_and_atmospheric_drag(tsince);
        let lpo = self.long_period_periodics(&sgaaduo);
        let keo = self.keplers_equation(lpo.xlt, sgaaduo.xnode, lpo.axn, lpo.ayn);
        let sppq = self.short_period_prelimenary_quantities(keo, lpo.axn, lpo.ayn, sgaaduo.a);

        let spo = self.short_periodics(sppq, sgaaduo.xnode, sgaaduo.xn);
        let ov = self.calculate_orientation_vectors(spo.uk, spo.xinck, spo.xnodek);

        self.calculate_position_and_velocity(ov, spo)
    }

    fn update_for_secular_gravity_and_atmospheric_drag(
        &self,
        tsince: f64,
    ) -> types::SecularGravityAndAtmosphericDragUpdateOutput {
        let xmdf = self.xmo + self.xmdot * tsince;
        let omgadf = self.omegao + self.omgdot * tsince;
        let xnoddf = self.xnodeo + self.xnodot * tsince;
        let tsq = tsince * tsince;
        let xnode = xnoddf + self.xnodcf * tsq;
        let mut tempa = 1.0 - self.c_constants.c1 * tsince;
        let mut tempe = self.bstar * self.c_constants.c4 * tsince;
        let mut templ = self.t2cof * tsq;
        let delomg = self.omgcof * tsince;
        let delm = self.xmcof * ((1.0 + self.eta * xmdf.cos()).powf(3.0) - self.delmo);
        let temp = delomg + delm;
        let xmp = xmdf + temp;
        let omega = omgadf - temp;
        let tcube = tsq * tsince;
        let tfour = tsince * tcube;
        tempa = tempa
            - self.d_constants.d2 * tsq
            - self.d_constants.d3 * tcube
            - self.d_constants.d4 * tfour;
        tempe = tempe + self.bstar * self.c_constants.c5 * (xmp.sin() - self.sinmo);
        templ = templ + self.t3cof * tcube + tfour * (self.t4cof + tsince * self.t5cof);
        let a = self.mmasmao.aodp * tempa.powf(2.0);
        let e = self.eo - tempe;
        let xl = xmp + omega + xnode + self.mmasmao.xnodp * templ;
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
        &self,
        sgaaduo: &types::SecularGravityAndAtmosphericDragUpdateOutput,
    ) -> types::LongPeriodicsOutput {
        let axn = sgaaduo.e * sgaaduo.omega.cos();
        let temp = 1.0 / (sgaaduo.a * sgaaduo.beta * sgaaduo.beta);
        let xll = temp * self.xlcof * axn;
        let aynl = temp * self.aycof;
        let xlt = sgaaduo.xl + xll;
        let ayn = sgaaduo.e * sgaaduo.omega.sin() + aynl;

        types::LongPeriodicsOutput { xlt, ayn, axn }
    }

    fn short_periodics(
        &self,
        sppq: types::ShortPeriodPrelimenaryQuantities,
        xnode: f64,
        xn: f64,
    ) -> types::ShortPeriodicsOutput {
        let rk = (sppq.r * (1.0 - 1.5 * sppq.temp2 * sppq.betal * self.mmasmao.x3thm1)
            + 0.5 * sppq.temp1 * self.x1mth2 * sppq.cos2u)
            * constants::XKMPER;
        let uk = sppq.u - 0.25 * sppq.temp2 * self.x7thm1 * sppq.sin2u;
        let xnodek = xnode + 1.5 * sppq.temp2 * self.mmasmao.cosio * sppq.sin2u;
        let xinck = self.xincl + 1.5 * sppq.temp2 * self.mmasmao.cosio * self.sinio * sppq.cos2u;
        let rdotk =
            (sppq.rdot - xn * sppq.temp1 * self.x1mth2 * sppq.sin2u) * constants::XKMPER / 60.0;
        let rfdotk = (sppq.rfdot
            + xn * sppq.temp1 * (self.x1mth2 * sppq.cos2u + 1.5 * self.mmasmao.x3thm1))
            * constants::XKMPER
            / 60.0;

        types::ShortPeriodicsOutput {
            rk,
            uk,
            xnodek,
            xinck,
            rdotk,
            rfdotk,
        }
    }

    fn keplers_equation(
        &self,
        xlt: f64,
        xnode: f64,
        axn: f64,
        ayn: f64,
    ) -> types::KeplersEquationOutput {
        let capu = helpers::fmod2p(xlt - xnode);
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

            if (epw - temp2).abs() <= self.e6a {
                break;
            }

            temp2 = epw;
        }

        types::KeplersEquationOutput {
            temp3,
            temp4,
            temp5,
            temp6,
            sinepw,
            cosepw,
        }
    }

    fn short_period_prelimenary_quantities(
        &self,
        keo: types::KeplersEquationOutput,
        axn: f64,
        ayn: f64,
        a: f64,
    ) -> types::ShortPeriodPrelimenaryQuantities {
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
        let u = helpers::actan(sinu, cosu);
        let sin2u = 2.0 * sinu * cosu;
        let cos2u = 2.0 * cosu * cosu - 1.0;
        let temp = 1.0 / pl;
        let temp1 = constants::CK2 * temp;
        temp2 = temp1 * temp;

        types::ShortPeriodPrelimenaryQuantities {
            r,
            rdot,
            rfdot,
            temp2,
            betal,
            temp1,
            cos2u,
            u,
            sin2u,
        }
    }

    fn calculate_orientation_vectors(
        &self,
        uk: f64,
        xinck: f64,
        xnodek: f64,
    ) -> types::OrientationVectors {
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
        &self,
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
}

fn sut() -> Sgp4 {
    Sgp4::new(
        test_constants::EO,
        test_constants::BSTAR,
        test_constants::XINCL,
        test_constants::OMEGAO,
        test_constants::XMO,
        test_constants::XNO,
        test_constants::XNODEO,
        test_constants::E6A,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_new() {
        let sgp4 = sut();

        assert_abs_diff_eq!(
            sgp4.xmdot,
            0.07006729335201786,
            epsilon = test_constants::SMALL_TOLERANCE
        );

        assert_abs_diff_eq!(
            sgp4.omgdot,
            -0.00002971792465285666,
            epsilon = test_constants::SMALL_TOLERANCE
        );
        assert_abs_diff_eq!(
            sgp4.omgcof,
            0.00000016348304905484922,
            epsilon = test_constants::SMALL_TOLERANCE
        );
        assert_abs_diff_eq!(
            sgp4.xmcof,
            -0.000049353388663657485,
            epsilon = test_constants::SMALL_TOLERANCE
        );
        assert_abs_diff_eq!(
            sgp4.xnodcf,
            -0.000000000002535821899421168,
            epsilon = test_constants::SMALL_TOLERANCE
        );
        assert_abs_diff_eq!(
            sgp4.xlcof,
            0.001935745758076399,
            epsilon = test_constants::SMALL_TOLERANCE
        );
        assert_abs_diff_eq!(
            sgp4.aycof,
            0.0011203600994653647,
            epsilon = test_constants::SMALL_TOLERANCE
        );
        assert_abs_diff_eq!(
            sgp4.delmo,
            0.6963086765241224,
            epsilon = test_constants::SMALL_TOLERANCE
        );
        assert_abs_diff_eq!(
            sgp4.sinmo,
            0.9362350466329594,
            epsilon = test_constants::SMALL_TOLERANCE
        );
        assert_abs_diff_eq!(
            sgp4.x7thm1,
            -0.3908964405664701,
            epsilon = test_constants::SMALL_TOLERANCE
        );
        assert_abs_diff_eq!(
            sgp4.t2cof,
            0.00000003500706632267481,
            epsilon = test_constants::SMALL_TOLERANCE
        );
        assert_abs_diff_eq!(
            sgp4.t3cof,
            0.00000000000008234434284266006,
            epsilon = test_constants::SMALL_TOLERANCE
        );
        assert_abs_diff_eq!(
            sgp4.t4cof,
            0.00000000000000000032351134586589164,
            epsilon = test_constants::SMALL_TOLERANCE
        );
        assert_abs_diff_eq!(
            sgp4.t5cof,
            0.0000000000000000000000015781283156000947,
            epsilon = test_constants::SMALL_TOLERANCE
        );

        let expected_c_constants = types::CConstants {
            c1: 2.3338044215116538e-8,
            c2: 0.0003492882575298811,
            c3: 0.004037532255765166,
            c4: 0.000377201121554739,
            c5: 0.012334919304344908,
        };

        assert_eq!(sgp4.c_constants, expected_c_constants);

        let expected_d_constants = types::DConstants {
            d2: 8.12550142270866e-14,
            d3: 4.2372075736327043e-19,
            d4: 2.5770097992217537e-24,
        };

        assert_eq!(sgp4.d_constants, expected_d_constants);
    }

    #[test]
    fn test_update_for_secular_gravity_and_atmospheric_drag() {
        let sgp4 = sut();

        let sgaaduo = sgp4.update_for_secular_gravity_and_atmospheric_drag(test_constants::TSINCE);

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

        let sgp4 = sut();

        let lpo = sgp4.long_period_periodics(&sgaaduo);

        assert_eq!(lpo.xlt, 0.07010542434717408);
        assert_eq!(lpo.ayn, 0.009741425556977382);
        assert_eq!(lpo.axn, -0.0003928808433640341);
    }

    #[test]
    fn test_keplers_equation() {
        let xlt = 0.07011593807103583;
        let xnode = 0.0;
        let axn = 0.005255942497390392;
        let ayn = 0.007976339600468509;

        let sgp4 = sut();

        let keplers_equation_output = sgp4.keplers_equation(xlt, xnode, axn, ayn);

        assert_eq!(keplers_equation_output.temp3, 0.0003281941220562471);
        assert_eq!(keplers_equation_output.temp4, 0.007960774282990234);
        assert_eq!(keplers_equation_output.temp5, 0.0052456858611741215);
        assert_eq!(keplers_equation_output.temp6, 0.0004980624833125527);
    }

    #[test]
    fn test_short_periodics() {
        let xnode = 2.02403913260325;
        let xn = 0.07010615556528188;

        let sppq = types::ShortPeriodPrelimenaryQuantities {
            r: 1.0430516048741472,
            rdot: 0.0006636054871099481,
            rfdot: 0.07271020474319065,
            temp2: 0.0005004479094222182,
            betal: 0.9999543754967234,
            temp1: 0.0005204771435457296,
            cos2u: 0.8538389429983284,
            u: 2.867852612718983,
            sin2u: -0.5205372795674639,
        };

        let sgp4 = sut();

        assert_eq!(
            types::ShortPeriodicsOutput {
                rk: 6657.7080462027025,
                uk: 2.867827155413039,
                xnodek: 2.0239238673191204,
                xinck: 1.2715395691080913,
                rdotk: 0.07238614054582597,
                rfdotk: 7.727982650878059
            },
            sgp4.short_periodics(sppq, xnode, xn,)
        )
    }

    #[test]
    fn test_calculate_orientation_vectors() {
        let uk = 2.867827155413039;
        let xinck = 1.2715395691080913;
        let xnodek = 2.0239238673191204;

        let sgp4 = sut();

        assert_eq!(
            types::OrientationVectors {
                ux: 0.3498156817129497,
                uy: -0.9004932019514819,
                uz: 0.25834276081762364,
                vx: 0.3735451525121787,
                vy: -0.11881912112969674,
                vz: -0.9199706709937114,
            },
            sgp4.calculate_orientation_vectors(uk, xinck, xnodek)
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

        let sgp4 = sut();

        let pav = sgp4.calculate_position_and_velocity(ov, spo);

        assert_abs_diff_eq!(
            test_constants::POSITION_AND_VELOCITY_0.x,
            pav.x,
            epsilon = test_constants::TOLERANCE
        );
        assert_abs_diff_eq!(
            test_constants::POSITION_AND_VELOCITY_0.y,
            pav.y,
            epsilon = test_constants::TOLERANCE
        );
        assert_abs_diff_eq!(
            test_constants::POSITION_AND_VELOCITY_0.z,
            pav.z,
            epsilon = test_constants::TOLERANCE
        );
        assert_abs_diff_eq!(
            test_constants::POSITION_AND_VELOCITY_0.xdot,
            pav.xdot,
            epsilon = test_constants::TOLERANCE
        );
        assert_abs_diff_eq!(
            test_constants::POSITION_AND_VELOCITY_0.ydot,
            pav.ydot,
            epsilon = test_constants::TOLERANCE
        );
        assert_abs_diff_eq!(
            test_constants::POSITION_AND_VELOCITY_0.zdot,
            pav.zdot,
            epsilon = test_constants::TOLERANCE
        );
    }
}
