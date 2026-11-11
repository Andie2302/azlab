//! Tests für Polynome und ihre Nullstellen.

use super::{
    BracketKind, Polynomial, RootBracket, clenshaw, horner, horner_with_derivative,
    isolate_real_roots, refine_root,
};
use crate::interval::Interval;
use crate::number::{Complex, Dual};
use crate::scalar::Real;
use crate::solve::{SolveError, Tolerance};
use crate::testutil::close;

fn coefficients<T: Real, const N: usize>(c: [f64; N]) -> [T; N] {
    c.map(T::from_f64)
}

fn interval<T: Real>(lo: f64, hi: f64) -> Interval<T> {
    Interval::new(T::from_f64(lo), T::from_f64(hi)).unwrap_or(Interval::ENTIRE)
}

fn empty_bracket<T: Real>() -> RootBracket<T> {
    RootBracket { interval: Interval::point(T::ZERO), kind: BracketKind::Unresolved }
}

fn evaluation<T: Real>() {
    let f = T::from_f64;
    // 1 + 2x + 3x²
    let p = Polynomial(coefficients::<T, 3>([1.0, 2.0, 3.0]));
    assert_eq!(p.eval(f(2.0)), f(17.0));
    assert_eq!(p.eval_with_derivative(f(2.0)), (f(17.0), f(14.0)));
    assert_eq!(p.derivative(), Polynomial(coefficients([2.0, 6.0, 0.0])));
    assert_eq!(p.degree(), Some(2));
    assert_eq!(Polynomial(coefficients::<T, 3>([1.0, 2.0, 0.0])).degree(), Some(1));
    assert_eq!(Polynomial::<T, 3>::ZERO.degree(), None);
    assert_eq!(p.coefficients().len(), 3);

    let q = Polynomial(coefficients::<T, 3>([1.0, 1.0, 1.0]));
    assert_eq!(p + q, Polynomial(coefficients([2.0, 3.0, 4.0])));
    assert_eq!(p - q, Polynomial(coefficients([0.0, 1.0, 2.0])));
    assert_eq!(-q, Polynomial(coefficients([-1.0, -1.0, -1.0])));
    assert_eq!(q * f(2.0), Polynomial(coefficients([2.0, 2.0, 2.0])));
    let mut acc = p;
    acc += q;
    acc -= q;
    assert_eq!(acc, p);

    // Auswertung in anderen Zahlbereichen über `map`.
    let dual = p.map(Dual::constant).eval(Dual::variable(f(2.0)));
    assert_eq!((dual.re, dual.eps), (f(17.0), f(14.0)));
    // 1 + 2i + 3i² = −2 + 2i
    let complex = p.map(Complex::from_real).eval(Complex::I);
    assert_eq!(complex, Complex::new(f(-2.0), f(2.0)));
    // Wertebereich über [0, 1] ist [1, 6]; der Einschluss enthält ihn.
    let range = p.map(Interval::point).eval(interval(0.0, 1.0));
    assert!(range.contains(f(1.0)) && range.contains(f(6.0)));

    // Slices.
    assert_eq!(horner(&coefficients::<T, 3>([1.0, 2.0, 3.0]), f(2.0)), f(17.0));
    assert_eq!(horner::<T>(&[], f(2.0)), T::ZERO);
    assert_eq!(
        horner_with_derivative(&coefficients::<T, 3>([1.0, 2.0, 3.0]), f(2.0)),
        (f(17.0), f(14.0))
    );
    // T₂(x) = 2x² − 1
    assert_eq!(clenshaw(&coefficients::<T, 3>([0.0, 0.0, 1.0]), f(0.5)), f(-0.5));
    // 3·T₀ + 2·T₁ = 3 + 2x
    assert_eq!(clenshaw(&coefficients::<T, 2>([3.0, 2.0]), f(0.5)), f(4.0));
    assert_eq!(clenshaw::<T>(&[], f(0.5)), T::ZERO);
}

fn isolation<T: Real>() {
    let f = T::from_f64;
    let tol = Tolerance::<T>::default();

    // (x − 1)(x − 2)(x − 3) = −6 + 11x − 6x² + x³
    let cubic = coefficients::<T, 4>([-6.0, 11.0, -6.0, 1.0]);
    let mut out = [empty_bracket::<T>(); 16];
    let found = isolate_real_roots(&cubic, interval(0.0, 4.5), f(1.0e-4), &mut out);
    assert!(found.is_ok_and(|count| {
        let brackets = out.get(..count).unwrap_or(&[]);
        count <= 6
            && [1.0, 2.0, 3.0].into_iter().all(|root| {
                brackets
                    .iter()
                    .any(|b| b.kind == BracketKind::Unique && b.interval.contains(f(root)))
            })
    }));
    // Verfeinern der ersten eindeutigen Klammer liefert die Nullstelle 1.
    assert!(out
        .iter()
        .find(|b| b.kind == BracketKind::Unique)
        .is_some_and(|b| refine_root(&cubic, b, &tol).is_ok_and(|s| close(s.x, T::ONE, 1024.0))));

    // Keine reelle Nullstelle: 1 + x².
    let mut out = [empty_bracket::<T>(); 4];
    assert_eq!(
        isolate_real_roots(&coefficients::<T, 3>([1.0, 0.0, 1.0]), interval(-3.0, 3.0), f(1.0e-3), &mut out),
        Ok(0)
    );
    // Nullstellen außerhalb des Suchintervalls zählen nicht.
    assert_eq!(isolate_real_roots(&cubic, interval(3.5, 4.5), f(1.0e-3), &mut out), Ok(0));

    // Doppelte Nullstelle (x − 1)²: wird nicht übersehen, bleibt aber
    // unaufgelöst.
    let double = coefficients::<T, 3>([1.0, -2.0, 1.0]);
    let mut out = [empty_bracket::<T>(); 32];
    let found = isolate_real_roots(&double, interval(0.0, 3.0), f(1.0e-3), &mut out);
    assert!(found.is_ok_and(|count| {
        let brackets = out.get(..count).unwrap_or(&[]);
        count >= 1
            && brackets.iter().any(|b| b.interval.contains(T::ONE))
            && brackets.iter().all(|b| b.kind == BracketKind::Unresolved)
    }));

    // Grad 7: x·(x² − 1)·(x² − 4)·(x² − 9) = −36x + 49x³ − 14x⁵ + x⁷ mit den
    // Nullstellen 0, ±1, ±2, ±3. Keine fehlt.
    let septic = coefficients::<T, 8>([0.0, -36.0, 0.0, 49.0, 0.0, -14.0, 0.0, 1.0]);
    let mut out = [empty_bracket::<T>(); 64];
    let found = isolate_real_roots(&septic, interval(-3.5, 3.6), f(1.0e-3), &mut out);
    assert!(found.is_ok_and(|count| {
        let brackets = out.get(..count).unwrap_or(&[]);
        [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
            .into_iter()
            .all(|root| brackets.iter().any(|b| b.interval.contains(f(root))))
    }));

    // Fehlerfälle.
    let mut none: [RootBracket<T>; 0] = [];
    assert_eq!(
        isolate_real_roots(&cubic, interval(0.0, 4.5), f(1.0e-3), &mut none),
        Err(SolveError::CapacityExceeded)
    );
    let mut out = [empty_bracket::<T>(); 4];
    assert_eq!(
        isolate_real_roots(&cubic, interval(0.0, 4.5), T::ZERO, &mut out),
        Err(SolveError::InvalidInput)
    );
    assert_eq!(
        isolate_real_roots(&coefficients::<T, 2>([0.0, 0.0]), interval(0.0, 1.0), f(0.1), &mut out),
        Err(SolveError::InvalidInput)
    );
    assert_eq!(
        isolate_real_roots(&cubic, Interval::ENTIRE, f(0.1), &mut out),
        Err(SolveError::InvalidInput)
    );
}

#[test]
fn evaluation_f32() {
    evaluation::<f32>();
}

#[test]
fn evaluation_f64() {
    evaluation::<f64>();
}

#[test]
fn isolation_f32() {
    isolation::<f32>();
}

#[test]
fn isolation_f64() {
    isolation::<f64>();
}

#[cfg(any(feature = "std", feature = "libm"))]
mod roots {
    use super::*;
    use crate::poly::{
        complex_roots, solve_cubic, solve_linear, solve_quadratic, solve_quadratic_complex,
        solve_quartic,
    };
    use crate::scalar::RealElementary;

    /// Gleich viele Nullstellen, jede relativ nahe am erwarteten Wert.
    fn roots_match<T: Real>(found: &[T], expected: &[f64], ulps: f64) -> bool {
        found.len() == expected.len()
            && found
                .iter()
                .zip(expected)
                .all(|(&x, &e)| close(x / T::from_f64(e), T::ONE, ulps))
    }

    fn closed_forms<T: RealElementary>() {
        let c = coefficients::<T, 2>;
        assert!(roots_match(solve_linear(c([-6.0, 3.0])).as_slice(), &[2.0], 8.0));
        assert!(solve_linear(c([1.0, 0.0])).is_empty());

        let c = coefficients::<T, 3>;
        // (x − 1)(x − 2), doppelte Nullstelle, keine reelle Nullstelle.
        assert!(roots_match(solve_quadratic(c([2.0, -3.0, 1.0])).as_slice(), &[1.0, 2.0], 8.0));
        assert!(roots_match(solve_quadratic(c([1.0, -2.0, 1.0])).as_slice(), &[1.0], 8.0));
        assert_eq!(solve_quadratic(c([1.0, 0.0, 1.0])).len(), 0);
        // Auslöschungsfall x² − 10⁸x + 1: Nullstellen 10⁻⁸ und 10⁸.
        assert!(roots_match(
            solve_quadratic(c([1.0, -1.0e8, 1.0])).as_slice(),
            &[1.0e-8, 1.0e8],
            16.0
        ));
        // Leitkoeffizient null, alles null.
        assert!(roots_match(solve_quadratic(c([-6.0, 3.0, 0.0])).as_slice(), &[2.0], 8.0));
        assert!(solve_quadratic(c([0.0, 0.0, 0.0])).is_empty());

        let c = coefficients::<T, 4>;
        assert!(roots_match(
            solve_cubic(c([-6.0, 11.0, -6.0, 1.0])).as_slice(),
            &[1.0, 2.0, 3.0],
            1024.0
        ));
        assert!(roots_match(solve_cubic(c([-1.0, 0.0, 0.0, 1.0])).as_slice(), &[1.0], 64.0));
        assert!(roots_match(solve_cubic(c([-8.0, 12.0, -6.0, 1.0])).as_slice(), &[2.0], 64.0));
        assert!(roots_match(solve_cubic(c([2.0, -3.0, 1.0, 0.0])).as_slice(), &[1.0, 2.0], 8.0));

        let c = coefficients::<T, 5>;
        // Biquadratischer Weg: Nullstellen 1, 2, 3, 4.
        assert!(roots_match(
            solve_quartic(c([24.0, -50.0, 35.0, -10.0, 1.0])).as_slice(),
            &[1.0, 2.0, 3.0, 4.0],
            1024.0
        ));
        // Weg über die Resolvente: Nullstellen 1, 2, 3, 5.
        assert!(roots_match(
            solve_quartic(c([30.0, -61.0, 41.0, -11.0, 1.0])).as_slice(),
            &[1.0, 2.0, 3.0, 5.0],
            1024.0
        ));
        assert!(solve_quartic(c([1.0, 0.0, 0.0, 0.0, 1.0])).is_empty());
        assert!(roots_match(
            solve_quartic(c([-6.0, 11.0, -6.0, 1.0, 0.0])).as_slice(),
            &[1.0, 2.0, 3.0],
            1024.0
        ));
    }

    /// Zu jedem erwarteten Wert gibt es eine berechnete Nullstelle in der
    /// Nähe.
    fn all_present<T: RealElementary>(found: &[Complex<T>], expected: &[(f64, f64)]) -> bool {
        let radius = T::EPSILON * T::from_f64(131_072.0);
        expected.iter().all(|&(re, im)| {
            let target = Complex::new(T::from_f64(re), T::from_f64(im));
            found.iter().any(|&z| (z - target).norm() <= radius)
        })
    }

    fn complex<T: RealElementary>() {
        let real = |v: f64| Complex::from_real(T::from_f64(v));

        // x² + 1 = 0 → ±i
        let pair = solve_quadratic_complex([real(1.0), real(0.0), real(1.0)]);
        assert!(all_present(&pair, &[(0.0, 1.0), (0.0, -1.0)]));
        // x² − 2x + 1 = 0 → doppelte Nullstelle 1
        let pair = solve_quadratic_complex([real(1.0), real(-2.0), real(1.0)]);
        assert!(all_present(&pair, &[(1.0, 0.0)]) && pair[0] == pair[1]);

        // Lockere Toleranz: Das Verfahren endet an der Rundungsgrenze.
        let loose = T::EPSILON * T::from_f64(4096.0);
        let tol = Tolerance::<T>::default().with_abs_x(loose).with_rel_x(loose);

        // x⁴ − 1 = 0 → ±1, ±i
        let mut roots = [Complex::new(T::ZERO, T::ZERO); 4];
        let quartic = [real(-1.0), real(0.0), real(0.0), real(0.0), real(1.0)];
        assert!(complex_roots(&quartic, &mut roots, &tol).is_ok());
        assert!(all_present(&roots, &[(1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)]));

        // (x − 1)(x − 2)(x − 3)(x − 4)(x − 5)
        let mut roots = [Complex::new(T::ZERO, T::ZERO); 5];
        let quintic = [
            real(-120.0),
            real(274.0),
            real(-225.0),
            real(85.0),
            real(-15.0),
            real(1.0),
        ];
        assert!(complex_roots(&quintic, &mut roots, &tol).is_ok());
        assert!(all_present(
            &roots,
            &[(1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0), (5.0, 0.0)]
        ));

        // Führende Nullen werden ignoriert: Grad 1, Nullstelle 2.
        let mut single = [Complex::new(T::ZERO, T::ZERO); 1];
        assert!(complex_roots(&[real(-2.0), real(1.0), real(0.0)], &mut single, &tol).is_ok());
        assert!(all_present(&single, &[(2.0, 0.0)]));

        // Fehlerfälle.
        let mut short = [Complex::new(T::ZERO, T::ZERO); 3];
        assert_eq!(
            complex_roots(&quartic, &mut short, &tol),
            Err(SolveError::CapacityExceeded)
        );
        let mut long = [Complex::new(T::ZERO, T::ZERO); 5];
        assert_eq!(complex_roots(&quartic, &mut long, &tol), Err(SolveError::InvalidInput));
        assert_eq!(
            complex_roots(&[real(0.0), real(0.0)], &mut short, &tol),
            Err(SolveError::InvalidInput)
        );
    }

    #[test]
    fn closed_forms_f32() {
        closed_forms::<f32>();
    }

    #[test]
    fn closed_forms_f64() {
        closed_forms::<f64>();
    }

    #[test]
    fn complex_f32() {
        complex::<f32>();
    }

    #[test]
    fn complex_f64() {
        complex::<f64>();
    }
}
