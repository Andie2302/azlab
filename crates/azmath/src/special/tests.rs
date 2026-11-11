//! Tests der speziellen Funktionen gegen bekannte Werte und Identitäten.

#![cfg(any(feature = "std", feature = "libm"))]

use super::{beta, erf, erfc, gamma, gamma_p, gamma_q, ln_gamma};
use crate::quad::GaussLegendre;
use crate::scalar::RealElementary;
use crate::testutil::close;

const ULPS: f64 = 2048.0;

fn check<T: RealElementary>() {
    let f = T::from_f64;
    let root_pi = T::PI.sqrt();

    // Γ(n) = (n − 1)!, Γ(1/2) = √π, Γ(x + 1) = x·Γ(x).
    let mut factorial = T::ONE;
    for n in 1..=10_u32 {
        let n = f(f64::from(n));
        assert!(close(gamma(n) / factorial, T::ONE, ULPS));
        assert!(close(ln_gamma(n), factorial.ln(), ULPS));
        factorial *= n;
    }
    assert!(close(gamma(f(0.5)), root_pi, ULPS));
    assert!(close(gamma(f(1.5)), root_pi * f(0.5), ULPS));
    assert!(close(gamma(f(-0.5)), -root_pi * f(2.0), ULPS));
    for x in [f(0.3), f(1.7), f(4.2), f(-1.3)] {
        assert!(close(gamma(x + T::ONE) / (x * gamma(x)), T::ONE, ULPS));
    }
    // Spiegelungsformel und Logarithmus.
    assert!(close(gamma(f(0.3)) * gamma(f(0.7)), T::PI / (T::PI * f(0.3)).sin(), ULPS));
    assert!(close(ln_gamma(f(4.2)), gamma(f(4.2)).ln(), ULPS));
    assert!(close(ln_gamma(f(-1.3)), gamma(f(-1.3)).abs().ln(), ULPS));
    // Große Argumente: Γ läuft über, ln Γ bleibt endlich und erfüllt die
    // Funktionalgleichung.
    assert!(ln_gamma(f(500.0)).is_finite());
    assert!(close(ln_gamma(f(501.0)) - ln_gamma(f(500.0)), f(500.0).ln(), 1.0e6));
    assert!(!gamma(T::ZERO).is_finite());

    // Beta: B(2, 3) = 1/12, Symmetrie, B(1/2, 1/2) = π.
    assert!(close(beta(f(2.0), f(3.0)), T::ONE / f(12.0), ULPS));
    assert!(close(beta(f(1.3), f(2.7)), beta(f(2.7), f(1.3)), ULPS));
    assert!(close(beta(f(0.5), f(0.5)), T::PI, ULPS));

    // Unvollständige Gammafunktion: P(1, x) = 1 − e^(−x), P + Q = 1.
    for x in [f(0.1), f(0.9), f(2.5), f(8.0)] {
        assert!(close(gamma_p(T::ONE, x), -(-x).exp_m1(), ULPS));
        assert!(close(gamma_q(T::ONE, x), (-x).exp(), ULPS * 4.0));
        for a in [f(0.5), f(2.3), f(6.0)] {
            assert!(close(gamma_p(a, x) + gamma_q(a, x), T::ONE, ULPS));
            assert!(gamma_p(a, x) > T::ZERO && gamma_p(a, x) < T::ONE);
        }
    }
    assert_eq!((gamma_p(f(2.0), T::ZERO), gamma_q(f(2.0), T::ZERO)), (T::ZERO, T::ONE));
    assert!(gamma_p(T::ZERO, T::ONE).is_nan() && gamma_q(T::ONE, -T::ONE).is_nan());
    assert!(gamma_p(T::NAN, T::ONE).is_nan());
    // Gegenprobe mit Quadratur: P(a, x) = (1/Γ(a))·∫₀ˣ t^(a−1)·e^(−t) dt.
    let rule = GaussLegendre::<T, 16>::new();
    let (a, x) = (f(2.5), f(1.8));
    let integral = rule.integrate_panels(|t| t.powf(a - T::ONE) * (-t).exp(), T::ZERO, x, 8);
    assert!((gamma_p(a, x) - integral / gamma(a)).abs() <= T::EPSILON.sqrt());

    // Fehlerfunktion: ungerade, erf + erfc = 1, Grenzwerte, Ableitung.
    assert_eq!(erf(T::ZERO), T::ZERO);
    assert_eq!(erfc(T::ZERO), T::ONE);
    for x in [f(0.2), f(0.8), f(1.5), f(3.0)] {
        assert!(close(erf(-x), -erf(x), ULPS));
        assert!(close(erf(x) + erfc(x), T::ONE, ULPS));
        assert!(close(erfc(-x), f(2.0) - erfc(x), ULPS));
        // erf(x) = (2/√π)·∫₀ˣ e^(−t²) dt
        let integral = rule.integrate(|t| (-t * t).exp(), T::ZERO, x);
        assert!((erf(x) - f(2.0) / root_pi * integral).abs() <= T::EPSILON.sqrt());
    }
    assert!(close(erf(f(6.0)), T::ONE, ULPS));
    // Für große x ist erfc winzig, aber nicht null und nicht ausgelöscht.
    let tail = erfc(f(5.0));
    assert!(tail > T::ZERO && tail < f(1.0e-10));
    // Bekannte Schranken: (1 − 1/(2x²)) < erfc(x)·x·√π·e^(x²) < 1.
    let ratio = tail * f(5.0) * root_pi * f(25.0).exp();
    assert!(ratio > f(0.98) && ratio < T::ONE);
}

#[test]
fn special_f32() {
    check::<f32>();
}

#[test]
fn special_f64() {
    check::<f64>();
}
