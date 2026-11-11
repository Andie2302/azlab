//! Tests der Quadratur gegen bekannte Integrale.

use super::{QuadError, QuadOptions, adaptive_simpson, simpson, trapezoid_rule};
use crate::scalar::Real;
use crate::testutil::close;

fn rules<T: Real>() {
    let f = T::from_f64;
    let third = T::ONE / f(3.0);
    // ∫₀¹ x² dx = 1/3: Simpson ist für Polynome bis Grad 3 exakt.
    assert!(close(simpson(|x: T| x * x, T::ZERO, T::ONE, 2), third, 8.0));
    assert!(close(simpson(|x: T| x * x * x, T::ZERO, f(2.0), 3), f(4.0), 8.0));
    // Trapezregel: exakt für Geraden, sonst wird sie mit mehr Intervallen besser.
    assert!(close(trapezoid_rule(|x: T| f(2.0) * x + T::ONE, T::ZERO, f(3.0), 1), f(12.0), 8.0));
    let coarse = (trapezoid_rule(|x: T| x * x, T::ZERO, T::ONE, 4) - third).abs();
    let fine = (trapezoid_rule(|x: T| x * x, T::ZERO, T::ONE, 8) - third).abs();
    assert!(fine < coarse * f(0.3));
    // Vertauschte Grenzen wechseln das Vorzeichen.
    assert!(close(simpson(|x: T| x * x, T::ONE, T::ZERO, 2), -third, 8.0));
    // Null Intervalle werden wie eines bzw. zwei behandelt.
    assert!(close(trapezoid_rule(|x: T| x, T::ZERO, T::ONE, 0), f(0.5), 8.0));
    assert!(close(simpson(|x: T| x, T::ZERO, T::ONE, 0), f(0.5), 8.0));

    // Adaptiv: 1/(1 + x²) über [0, 1] ohne Elementarfunktionen prüfbar über
    // den Vergleich mit einer feinen Simpsonsumme.
    let options = QuadOptions::<T>::default();
    let reference = simpson(|x: T| T::ONE / (T::ONE + x * x), T::ZERO, T::ONE, 2000);
    assert!(adaptive_simpson(|x: T| T::ONE / (T::ONE + x * x), T::ZERO, T::ONE, &options)
        .is_ok_and(|q| close(q.value, reference, 4096.0) && q.evaluations >= 5));
    // Polynom Grad 3: sofort exakt.
    assert!(adaptive_simpson(|x: T| x * x * x, T::ZERO, f(2.0), &options)
        .is_ok_and(|q| close(q.value, f(4.0), 8.0) && q.evaluations == 5));
    // Knick: |x − 0,3| über [0, 1] = (0,3² + 0,7²)/2.
    assert!(adaptive_simpson(|x: T| (x - f(0.3)).abs(), T::ZERO, T::ONE, &options)
        .is_ok_and(|q| close(q.value, (f(0.3) * f(0.3) + f(0.7) * f(0.7)) * f(0.5), 4096.0)));

    // Fehlerfälle.
    assert_eq!(
        adaptive_simpson(|x: T| x, T::ZERO, T::NAN, &options),
        Err(QuadError::InvalidInput)
    );
    assert_eq!(
        adaptive_simpson(|x: T| x, T::ZERO, T::ONE, &options.with_rel_tol(T::ZERO)),
        Err(QuadError::InvalidInput)
    );
    assert!(matches!(
        adaptive_simpson(
            |x: T| (x - f(0.3)).abs(),
            T::ZERO,
            T::ONE,
            &options.with_max_evaluations(7)
        ),
        Err(QuadError::NotConverged(q)) if q.evaluations <= 9
    ));
    assert!(adaptive_simpson(|x: T| x, T::ZERO, T::ONE, &options.with_abs_tol(f(0.1))).is_ok());
}

#[test]
fn rules_f32() {
    rules::<f32>();
}

#[test]
fn rules_f64() {
    rules::<f64>();
}

#[cfg(any(feature = "std", feature = "libm"))]
mod gauss {
    use super::*;
    use crate::quad::GaussLegendre;
    use crate::scalar::RealElementary;

    fn check<T: RealElementary>() {
        let f = T::from_f64;
        // Zwei Stützstellen: ±1/√3, Gewichte 1.
        let two = GaussLegendre::<T, 2>::new();
        let node = (T::ONE / f(3.0)).sqrt();
        assert!(close(two.nodes()[1], node, 8.0) && close(two.nodes()[0], -node, 8.0));
        assert!(close(two.weights()[0], T::ONE, 8.0) && close(two.weights()[1], T::ONE, 8.0));
        // Drei Stützstellen: 0 und ±√(3/5), Gewichte 8/9 und 5/9.
        let three = GaussLegendre::<T, 3>::default();
        assert!(close(three.nodes()[1], T::ZERO, 8.0));
        assert!(close(three.nodes()[2], (f(3.0) / f(5.0)).sqrt(), 8.0));
        assert!(close(three.weights()[1], f(8.0) / f(9.0), 8.0));
        assert!(close(three.weights()[0], f(5.0) / f(9.0), 8.0));

        // N Stützstellen integrieren Grad 2N − 1 exakt: ∫₀¹ x⁵ dx = 1/6.
        assert!(close(three.integrate(|x| x.powi(5), T::ZERO, T::ONE), T::ONE / f(6.0), 16.0));
        // Gewichte summieren sich zur Intervalllänge.
        let ten = GaussLegendre::<T, 10>::new();
        assert!(close(ten.weights().iter().fold(T::ZERO, |s, &w| s + w), f(2.0), 16.0));
        // ∫₀^π sin x dx = 2, ∫₀¹ eˣ dx = e − 1.
        assert!(close(ten.integrate(|x| x.sin(), T::ZERO, T::PI), f(2.0), 64.0));
        assert!(close(ten.integrate(|x| x.exp(), T::ZERO, T::ONE), T::E - T::ONE, 64.0));
        // Viele Teilintervalle gleichen die niedrige Ordnung aus.
        let panels = two.integrate_panels(|x| x.sin(), T::ZERO, T::PI, 64);
        assert!((panels - f(2.0)).abs() < f(1.0e-5));
        assert!(close(two.integrate_panels(|x| x * x, T::ZERO, T::ONE, 0), T::ONE / f(3.0), 8.0));
        // Entartete Regeln.
        assert!(close(GaussLegendre::<T, 1>::new().integrate(|x| x, T::ZERO, f(2.0)), f(2.0), 8.0));
        assert_eq!(GaussLegendre::<T, 0>::new().integrate(|x| x, T::ZERO, T::ONE), T::ZERO);
    }

    #[test]
    fn gauss_f32() {
        check::<f32>();
    }

    #[test]
    fn gauss_f64() {
        check::<f64>();
    }
}
