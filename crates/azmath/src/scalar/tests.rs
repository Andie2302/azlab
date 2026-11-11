//! Tests der Skalar-Traits. Bewusst generisch geschrieben: neue Zahlentypen
//! und Backends durchlaufen dieselben Prüfungen.

use super::*;

fn poly<T: Ring>(x: T) -> T {
    x * x + x + T::ONE
}

#[test]
fn ring_works_for_ints_and_floats() {
    assert_eq!(poly(3_i32), 13);
    assert_eq!(poly(3.0_f64), 13.0);
}

#[test]
fn powi_uses_repeated_squaring() {
    assert_eq!(Field::powi(2.0_f64, 10), 1024.0);
    assert_eq!(Field::powi(2.0_f32, -2), 0.25);
    assert_eq!(Field::powi(5.0_f64, 0), 1.0);
    assert_eq!(Field::powi(2.0_f64, i32::MIN), 0.0);
}

#[test]
fn real_basics() {
    assert_eq!(Real::abs(-2.5_f32), 2.5);
    assert_eq!(Real::max(1.0_f64, 2.0), 2.0);
    assert_eq!(<f32 as Field>::from_f64(0.5), 0.5);
    assert!(Real::is_nan(<f64 as Real>::NAN));
}

#[test]
fn conjugate_is_identity_for_reals() {
    assert_eq!(Conjugate::conj(-3.0_f64), -3.0);
    assert_eq!(Conjugate::norm_sqr(-3.0_f64), 9.0);
    assert_eq!(Conjugate::scale(2.0_f32, 4.0), 8.0);
}

#[cfg(any(feature = "std", feature = "libm"))]
mod backend {
    use super::*;

    fn close<T: Real>(a: T, b: T) -> bool {
        (a - b).abs() <= T::from_f64(8.0) * T::EPSILON * b.abs().max(T::ONE)
    }

    fn check<T: RealElementary>() {
        let two = T::from_f64(2.0);
        let three = T::from_f64(3.0);
        let four = T::from_f64(4.0);
        let five = T::from_f64(5.0);
        assert!(close((three * three + four * four).sqrt(), five));
        assert!(close(T::ONE.exp(), T::E));
        assert!(close(T::E.ln(), T::ONE));
        assert!(close(T::ONE.atan2(T::ONE), T::PI / four));
        assert!(close(two.powf(T::from_f64(0.5)), two.sqrt()));
        assert!(close(T::from_f64(27.0).cbrt(), three));
        assert!(close(three.hypot(four), five));
        assert!(close(two.mul_add(three, T::ONE), T::from_f64(7.0)));
        assert_eq!(T::from_f64(2.5).round(), three);
        assert_eq!(T::from_f64(-2.5).floor(), -three);
        let (s, c) = T::from_f64(0.3).sin_cos();
        assert!(close(s * s + c * c, T::ONE));
    }

    #[test]
    fn f32_backend() {
        check::<f32>();
    }

    #[test]
    fn f64_backend() {
        check::<f64>();
    }
}
