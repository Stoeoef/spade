use crate::{DelaunayTriangulation, Point2, Triangulation, TriangulationExt};

use alloc::{vec, vec::Vec};

fn fuzz_test(vertices: Vec<Point2<f64>>) {
    let mut clone = vertices.clone();
    clone.sort_by(|l, r| l.partial_cmp(r).unwrap());
    clone.dedup();
    let expected_size = clone.len();

    let triangulation = DelaunayTriangulation::<_>::bulk_load(vertices).unwrap();
    triangulation.sanity_check();
    assert_eq!(triangulation.num_vertices(), expected_size);
}

#[test]
fn test_bulk_load_fuzz_1() {
    fuzz_test(vec![
        Point2::new(10.0, 0.0),
        Point2::new(0.0, 2.0),
        Point2::new(0.0, 16.0),
        Point2::new(0.0, 0.0),
        Point2::new(0.0, 16.0),
        Point2::new(4.0, 4.0),
    ]);
}

#[test]
fn test_bulk_load_fuzz_2() {
    fuzz_test(vec![
        Point2::new(10.0, 1.0),
        Point2::new(1.0, 1.0),
        Point2::new(3.0, 2.0),
        Point2::new(0.0, 0.0),
        Point2::new(-19.0, 0.0),
        Point2::new(0.0, -1.0),
        Point2::new(-1.0, -1.0),
        Point2::new(-1.0, 0.0),
    ]);
}

#[test]
fn test_bulk_load_fuzz_3() {
    fuzz_test(vec![
        Point2::new(10.0, 3.0),
        Point2::new(0.0, -19.0),
        Point2::new(16.0, 0.0),
        Point2::new(0.0, 10.0),
        Point2::new(0.0, 8.0),
        Point2::new(0.0, 0.0),
        Point2::new(0.0, 1.0),
        Point2::new(-4.0, 0.0),
        Point2::new(-20.0, 0.0),
        Point2::new(8.0, 0.0),
        Point2::new(0.0, 0.0),
        Point2::new(-7.0, -1.0),
        Point2::new(0.0, 4.0),
    ]);
}

#[test]
fn test_bulk_load_fuzz_4() {
    fuzz_test(vec![
        Point2::new(10.0, 1.0),
        Point2::new(1.0, 3.0),
        Point2::new(0.0, 0.0),
        Point2::new(3.0, 0.0),
        Point2::new(-4.0, 0.0),
        Point2::new(0.0, 0.0),
        Point2::new(16.0, -19.0),
        Point2::new(-1.0, 10.0),
        Point2::new(8.0, 0.0),
        Point2::new(0.0, 0.0),
        Point2::new(0.0, -4.0),
    ]);
}

#[test]
fn test_bulk_load_fuzz_5() {
    fuzz_test(vec![
        Point2::new(10.0, 16.0),
        Point2::new(0.0, -4.0),
        Point2::new(7.0, -19.0),
        Point2::new(4.0, 0.0),
        Point2::new(16.0, 0.0),
        Point2::new(-20.0, 10.0),
        Point2::new(0.0, 0.0),
        Point2::new(0.0, -1.0),
        Point2::new(-1.0, 0.0),
    ]);
}

#[test]
fn test_bulk_load_fuzz_6() {
    fuzz_test(vec![
        Point2::new(10.0, 0.0),
        Point2::new(16.0, 16.0),
        Point2::new(0.0, 0.0),
        Point2::new(16.0, 0.0),
        Point2::new(0.0, 10.0),
        Point2::new(-4.0, 0.0),
        Point2::new(0.0, 0.0),
        Point2::new(-19.0, 0.0),
        Point2::new(0.0, -20.0),
        Point2::new(10.0, 0.0),
        Point2::new(8.0, 0.0),
        Point2::new(0.0, 1.0),
    ]);
}

#[test]
fn test_bulk_load_fuzz_7() {
    fuzz_test(vec![
        Point2::new(10.0, -1.0),
        Point2::new(-1.0, -1.0),
        Point2::new(-1.0, -1.0),
        Point2::new(-1.0, -1.0),
        Point2::new(-1.0, -1.0),
        Point2::new(-1.0, -1.0),
        Point2::new(-1.0, -1.0),
        Point2::new(-1.0, -1.0),
        Point2::new(-1.0, -1.0),
        Point2::new(-1.0, -1.0),
        Point2::new(-1.0, -1.0),
        Point2::new(-1.0, -2.0),
        Point2::new(-1.0, -1.0),
        Point2::new(-1.0, -1.0),
        Point2::new(0.0, -7.0),
        Point2::new(0.0, 2.0),
        Point2::new(0.0, 16.0),
        Point2::new(0.0, 0.0),
        Point2::new(3.0, 0.0),
        Point2::new(16.0, 0.0),
    ]);
}

#[test]
fn test_bulk_load_fuzz_8() {
    fuzz_test(vec![
        Point2::new(10.0, 3.0),
        Point2::new(2.0, 0.0),
        Point2::new(0.0, 4.0),
        Point2::new(0.0, 0.0),
        Point2::new(10.0, 2.0),
        Point2::new(8.0, 0.0),
        Point2::new(0.0, -19.0),
        Point2::new(0.0, 0.0),
        Point2::new(-20.0, 10.0),
        Point2::new(0.0, -1.0),
        Point2::new(-1.0, 0.0),
        Point2::new(1.0, 0.0),
    ]);
}

#[test]
fn test_bulk_load_fuzz_9() {
    fuzz_test(vec![
        Point2::new(10.0, 3.0),
        Point2::new(2.0, 0.0),
        Point2::new(16.0, 9.0),
        Point2::new(0.0, 2.0),
        Point2::new(0.0, -19.0),
        Point2::new(0.0, 0.0),
        Point2::new(-20.0, 10.0),
        Point2::new(0.0, 3.0),
        Point2::new(1.0, 0.0),
        Point2::new(0.0, -1.0),
        Point2::new(-1.0, -1.0),
    ]);
}

#[test]
fn test_bulk_load_fuzz_10() {
    fuzz_test(vec![
        Point2::new(0.0, 1.469368e-39),
        Point2::new(1.83671e-40, 0.0),
        Point2::new(0.03515625, 0.0),
        Point2::new(0.0, 0.0),
        Point2::new(0.0, 9.1835e-41),
        Point2::new(1.85823e-40, 0.0),
    ]);
}

#[test]
fn test_bulk_load_fuzz_11() {
    fuzz_test(vec![
        Point2::new(0.0, -4.391866e-18),
        Point2::new(0.0, 9.403955e-38),
        Point2::new(0.0, 0.0),
        Point2::new(9.408547e-38, 0.0),
    ]);
}

#[test]
fn test_bulk_load_fuzz_12() {
    fuzz_test(vec![
        Point2::new(12.078431, 12.078431),
        Point2::new(12.078431, 12.078431),
        Point2::new(12.078431, 12.111878),
        Point2::new(3.0196078, 3.0196078),
        Point2::new(12.078431, 12.078919),
        Point2::new(12.078431, 12.078431),
        Point2::new(12.124999, 12.078431),
        Point2::new(12.078431, 12.078431),
        Point2::new(12.098451, 12.078431),
        Point2::new(12.07843, 12.07843),
        Point2::new(12.078431, 12.078431),
        Point2::new(12.078431, 12.015931),
        Point2::new(12.015931, 12.078431),
    ]);
}

#[test]
fn test_bulk_load_fuzz_13() {
    fuzz_test(vec![
        Point2::new(-1.7936620343357662e-43, 1.7936620343357662e-43),
        Point2::new(-1.793662034335766e-43, 1.7936620343357662e-43),
        Point2::new(1.793662034335766e-43, -1.793662034335766e-43),
        Point2::new(1.793662034335766e-43, 1.7936620343357662e-43),
        Point2::new(1.7936620343357662e-43, -1.793662034335766e-43),
        Point2::new(-1.793662034335766e-43, -1.793662034335766e-43),
        Point2::new(-1.793662034335766e-43, -1.7936620343357662e-43),
        Point2::new(1.793662034335766e-43, -1.7936620343357662e-43),
        Point2::new(1.793662034335766e-43, 1.793662034335766e-43),
        Point2::new(-1.7936620343357662e-43, -1.7936620343357662e-43),
        Point2::new(-1.7936620343357662e-43, -1.793662034335766e-43),
        Point2::new(1.7936620343357662e-43, 1.793662034335766e-43),
    ]);
}

#[test]
fn test_bulk_load_fuzz_14() {
    fuzz_test(vec![
        Point2::new(-1.999_999_739_220_655_4, 4.000_466_406_403_405),
        Point2::new(-1.999_877_928_659_842_3, 4.062_966_407_334_729),
        Point2::new(-1.999_999_998_972_344, -1.999_999_998_981_323_3),
        Point2::new(-1.999_999_739_220_655_2, 4.062_966_406_403_405),
        Point2::new(-1.999_999_872_312_472_1, 4.062_966_407_334_727),
        Point2::new(-1.999_999_739_220_655_4, 4.062_966_406_403_405),
        Point2::new(-1.999_877_928_659_842_3, 4.062_966_407_334_727),
        Point2::new(-1.999_999_998_972_342_3, -1.999_999_996_711_224_5),
    ]);
}

#[test]
fn test_bulk_load_fuzz_15() {
    fuzz_test(vec![
        Point2::new(-1.9999999964929343, 1.9999999999999998),
        Point2::new(1.9999999999999998, 1.9999998807907102),
        Point2::new(1.9999999999999998, -2.0009765624999996),
        Point2::new(1.9999999999955662, 1.9999999973661031),
        Point2::new(1.9999999839783411, -1.9999999973660463),
        Point2::new(1.9999999973660552, 1.9999999999938607),
        Point2::new(1.9999999973661031, 1.9999999999999998),
        Point2::new(-1.750244140624953, 1.9999999999999998),
        Point2::new(1.9999990463547872, 1.9999999999999998),
        Point2::new(1.9999999973660465, 1.9999990463547872),
        Point2::new(1.9999999999999998, 1.999999997366057),
        Point2::new(1.9999999973661031, 1.999999217682898),
        Point2::new(-2.2499999999999996, -1.9399394989086431),
        Point2::new(1.999999999985449, 1.9999999999956228),
        Point2::new(1.9999999999956228, 2.2499999947322062),
        Point2::new(1.9999999850944616, 1.9999999899154657),
        Point2::new(1.9999999999999458, -1.9999999973661031),
        Point2::new(1.9843749999999998, 1.9999990464712025),
        Point2::new(1.9999999999999998, 1.9999999973661031),
        Point2::new(1.999999110747012, 1.9999999999999998),
        Point2::new(-1.999999999999953, 1.9999999999999998),
        Point2::new(1.9531240463547872, 1.9999999999999998),
        Point2::new(1.9999999973661022, 1.9999999973661031),
    ]);
}

#[test]
fn test_bulk_load_fuzz_16() {
    fuzz_test(vec![
        Point2::new(-1.9999999999999458, 1.9999999999999998),
        Point2::new(1.9999999999956228, 1.9843749973660463),
        Point2::new(-1.9999999999999998, 1.9999999999999998),
        Point2::new(1.9999998807907102, 1.9999999969149937),
    ]);
}

#[test]
fn test_bulk_load_fuzz_17() {
    fuzz_test(vec![
        Point2::new(2.0000000074796844, 2.0),
        Point2::new(2.0, 2.0000000298023224),
        Point2::new(2.001953125, 8.0),
        Point2::new(2.0, 2.000000029802323),
        Point2::new(2.0000000000000004, 2.0000000298023224),
    ]);
}

#[test]
fn test_bulk_load_fuzz_18() {
    fuzz_test(vec![
        Point2::new(2.001953125, -2.0),
        Point2::new(2.0000000298023224, 2.0001220703125),
        Point2::new(2.0019569396972656, -1.9999971389770508),
        Point2::new(2.0000000298023224, 2.0),
        Point2::new(2.0000000298023224, -2.001953125),
        Point2::new(2.124999999709075, 3.000000029802323),
        Point2::new(2.0000000000001137, 2.0),
        Point2::new(2.001953125029104, -2.0004883110523224),
        Point2::new(2.001953125, 2.0),
        Point2::new(2.0, -2.0004883110523224),
    ]);
}
