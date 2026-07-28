use crate::constants;
use crate::test_constants;
use crate::types;

pub fn recover_original_mean_motion_and_semimajor_axis(
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

pub fn adjust_atmospheric_drag_for_low_orbit(aodp: f64, eo: f64) -> (f64, f64) {
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

pub fn calculate_c_constants(
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

pub fn calculate_d_constants(
    aodp: f64,
    tsi: f64,
    c1sq: f64,
    s4: f64,
    c1: f64,
) -> types::DConstants {
    let d2 = 4.0 * aodp * tsi * c1sq;
    let temp = d2 * tsi * c1 / 3.0;
    let d3 = (17.0 * aodp + s4) * temp;
    let d4 = 0.5 * temp * aodp * tsi * (221.0 * aodp + 31.0 * s4) * c1;

    types::DConstants { d2, d3, d4 }
}

pub fn fmod2p(x: f64) -> f64 {
    let rev = x / constants::TWOPI;
    let mut temp = x - (rev.trunc()) * constants::TWOPI;
    if temp < 0.0 {
        temp += constants::TWOPI;
    }
    temp
}

pub fn actan(sinx: f64, cosx: f64) -> f64 {
    let angle = sinx.atan2(cosx);
    if angle < 0.0 {
        angle + constants::TWOPI
    } else {
        angle
    }
}

fn format_tle_split_with_negative_exponent(input: &str) -> f64 {
    assert!(
        input.chars().any(|c| matches!(c, '-')),
        "value is in wrong format"
    );

    let exponent = input.chars().last().unwrap().to_digit(10).unwrap() as i32;
    let base = format!("0.{}", &input[..input.len() - 2])
        .parse::<f64>()
        .unwrap();

    base * 10.0_f64.powi(-exponent)
}

pub fn text_to_tle(input: &String) -> types::TLE {
    let mut row = 0;

    let mut name = "".to_string();

    let mut first_row: Vec<&str> = vec![];
    let mut second_row: Vec<&str> = vec![];

    for split in input.split_whitespace() {
        if row == 0 && split.len() == 1 && split == "1" {
            row = 1;
            continue;
        } else if row == 1 && split.len() == 1 && split == "2" {
            row = 2;
            continue;
        }

        if row == 0 {
            name += split;
        } else if row == 1 {
            first_row.push(split);
        } else if row == 2 {
            second_row.push(split);
        }
    }

    let fdomm = first_row[3];
    let formatted_first_derivative_of_mean_motion = format!("0{fdomm}");

    types::TLE {
        name,
        number: first_row[0].to_string(),
        international_designator: first_row[1].to_string(),
        epoch_year_julian_fraction: first_row[2].parse::<f64>().unwrap(),
        first_derivative_of_mean_motion: formatted_first_derivative_of_mean_motion
            .parse::<f64>()
            .unwrap(),
        second_derivative_of_mean_motion: format_tle_split_with_negative_exponent(first_row[4]),
        drag_term: format_tle_split_with_negative_exponent(first_row[5]),
        ephemeris_type: first_row[6].parse::<u32>().unwrap(),
        element_number_check_sum: first_row[7].parse::<u32>().unwrap(),
        inclination: second_row[1].parse::<f64>().unwrap(),
        right_ascension_of_ascending_node: second_row[2].parse::<f64>().unwrap(),
        eccentricity: format!("0.{}", second_row[3]).parse::<f64>().unwrap(),
        argument_of_perigee: second_row[4].parse::<f64>().unwrap(),
        mean_anomaly: second_row[5].parse::<f64>().unwrap(),
        mean_motion: second_row[6][..10].parse::<f64>().unwrap(),
    }
}

pub fn validate_tle(input: &String) -> types::ValidationResponse {
    let mut valid = true;
    let mut message = String::new();

    if input.is_empty() {
        return types::ValidationResponse {
            valid: false,
            message: "Input field is empty".to_string(),
        };
    }

    let mut row = 0;
    let mut first_row_columns = 0;
    let mut second_row_columns = 0;

    if input.chars().nth(0) == Some('1') {
        valid = false;
        message = "Satellite name is missing, needs to be added before first row".to_string();
    }

    for split in input.split_whitespace() {
        if row == 0 && split.len() == 1 && split == "1" {
            row = 1;
            continue;
        } else if row == 1 && split.len() == 1 && split == "2" {
            row = 2;
            continue;
        }

        if row == 1 {
            first_row_columns += 1;
        } else if row == 2 {
            second_row_columns += 1;
        }
    }

    if row != 2 {
        valid = false;
        message = "Row missing".to_string();
    } else if first_row_columns != 8 {
        valid = false;
        message = format!(
            "Wrong number of columns in row 1. Should be 8, is {}",
            first_row_columns
        );
    } else if second_row_columns != 7 {
        valid = false;
        message = format!(
            "Wrong number of columns in row 2. Should be 7, is {}",
            second_row_columns
        );
    }

    if valid {
        let previous_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {})); // suppress default panic output
        let result = std::panic::catch_unwind(|| text_to_tle(input));
        std::panic::set_hook(previous_hook);

        if result.is_err() {
            valid = false;
            message = "Failed to parse TLE data".to_string();
        }
    }

    types::ValidationResponse {
        valid,
        message: message.to_string(),
    }
}

use approx::assert_abs_diff_eq;

pub fn assert_approx(first: f64, second: f64) {
    assert_abs_diff_eq!(first, second, epsilon = test_constants::MID_TOLERANCE);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_recover_original_mean_motion_and_semimajor_axis() {
        let mmasmao = recover_original_mean_motion_and_semimajor_axis(
            test_constants::XNO,
            test_constants::XINCL,
            test_constants::EO,
        );

        assert_approx(mmasmao.xnodp, 0.07010615558630984);
        assert_approx(mmasmao.aodp, 1.040117522759639);
        assert_approx(mmasmao.betao2, 0.99992477733639);
        assert_approx(mmasmao.betao, 0.9999623879608622);
        assert_approx(mmasmao.x3thm1, -0.73895561738563);
        assert_approx(mmasmao.theta2, 0.08701479420478998);
        assert_approx(mmasmao.cosio, 0.29498270153483575);
    }

    #[test]
    fn test_adjust_atmospheric_drag_for_low_orbit() {
        let aodp = 1.040117522759639; // From previous test, same sample test case
        let (s4, qoms24) = adjust_atmospheric_drag_for_low_orbit(aodp, test_constants::EO);

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
            eta,
            coef,
            xnodp,
            aodp,
            eeta,
            tsi,
            x3thm1,
            test_constants::BSTAR,
            a3ovk2,
            sinio,
            test_constants::EO,
            betao2,
            theta2,
            test_constants::OMEGAO,
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
    fn test_format_tle_split_with_negative_exponent() {
        assert_eq!(
            format_tle_split_with_negative_exponent("13844-3"),
            0.00013844
        );
        assert_eq!(
            format_tle_split_with_negative_exponent("66816-4"),
            0.000066816
        );
    }

    #[test]
    fn test_text_to_tle() {
        let input = test_constants::SGP4TLE.to_string();

        let tle = text_to_tle(&input);

        assert_eq!(
            tle,
            types::TLE {
                name: "SGP4(SGP4)".to_string(),
                number: "88888U".to_string(),
                international_designator: "98067A".to_string(),
                epoch_year_julian_fraction: 80275.98708465,
                first_derivative_of_mean_motion: 0.00073094,
                second_derivative_of_mean_motion: 0.00013844,
                drag_term: 0.000066816,
                ephemeris_type: 0,
                element_number_check_sum: 8,
                inclination: 72.8435,
                right_ascension_of_ascending_node: 115.9689,
                eccentricity: 0.0086731,
                argument_of_perigee: 52.6988,
                mean_anomaly: 110.5714,
                mean_motion: 16.0582451,
            }
        )
    }

    #[test]
    fn test_text_to_tle_2() {
        let input = "ISS (ZARYA)
1 25544U 98067A   26209.15252568  .00010831  00000-0  20282-3 0  9993
2 25544  51.6320  97.3682 0007093 345.6120  14.4666 15.49220842578109"
            .to_string();

        let tle = text_to_tle(&input);

        assert_eq!(
            tle,
            types::TLE {
                name: "ISS(ZARYA)".to_string(),
                number: "25544U".to_string(),
                international_designator: "98067A".to_string(),
                epoch_year_julian_fraction: 26209.15252568,
                first_derivative_of_mean_motion: 0.00010831,
                second_derivative_of_mean_motion: 0.0,
                drag_term: 0.00020282,
                ephemeris_type: 0,
                element_number_check_sum: 9993,
                inclination: 51.6320,
                right_ascension_of_ascending_node: 97.3682,
                eccentricity: 0.0007093,
                argument_of_perigee: 345.6120,
                mean_anomaly: 14.4666,
                mean_motion: 15.4922084,
            }
        )
    }
}
