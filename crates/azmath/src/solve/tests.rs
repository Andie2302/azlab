//! Tests der Löser.

use super::{
    ScalarResult, Solution, SolveError, Tolerance, bisection, brent, newton, newton_auto,
    newton_bracketed, secant,
};
use crate::number::Dual;
use crate::scalar::Real;
use crate::testutil::close;

fn scalar<T: Real>() {
    let one = T::ONE;
    let two = T::from_f64(2.0);
    let three = T::from_f64(3.0);
    let tol = Tolerance::<T>::default();

    // x² − 2 = 0 → x = √2; geprüft über Rückeinsetzen.
    let f = |x: T| x * x - two;
    let df = |x: T| (x * x - two, x + x);
    let is_root = |result: ScalarResult<T>| result.is_ok_and(|s| close(s.x * s.x, two, 64.0));

    assert!(is_root(bisection(f, one, two, &tol)));
    assert!(is_root(brent(f, one, two, &tol)));
    assert!(is_root(brent(f, two, one, &tol)));
    assert!(is_root(secant(f, one, two, &tol)));
    assert!(is_root(newton(df, two, &tol)));
    assert!(is_root(newton_auto(|x: Dual<T>| x * x - Dual::constant(two), two, &tol)));
    assert!(is_root(newton_bracketed(df, one, two, &tol)));

    // x³ − x − 2 = 0 auf [1, 2].
    assert!(brent(|x: T| x * x * x - x - two, one, two, &tol)
        .is_ok_and(|s| close(s.x * s.x * s.x - s.x, two, 64.0)));

    // Nullstelle genau am Rand.
    assert_eq!(
        bisection(|x: T| x - one, one, two, &tol),
        Ok(Solution { x: one, residual: T::ZERO, iterations: 0 })
    );

    // Kein Vorzeichenwechsel.
    assert_eq!(bisection(f, two, three, &tol), Err(SolveError::NoBracket));
    assert_eq!(brent(f, two, three, &tol), Err(SolveError::NoBracket));
    assert_eq!(newton_bracketed(df, two, three, &tol), Err(SolveError::NoBracket));

    // Ungültige Eingabe.
    assert_eq!(bisection(f, T::NAN, two, &tol), Err(SolveError::InvalidInput));
    assert_eq!(brent(f, one, T::INFINITY, &tol), Err(SolveError::InvalidInput));

    // Nicht konvergiert: Der Fehler enthält den letzten Stand.
    let single = tol.with_max_iterations(1);
    assert_eq!(
        bisection(f, one, two, &single),
        Err(SolveError::NotConverged(Solution {
            x: T::from_f64(1.5),
            residual: T::from_f64(0.25),
            iterations: 1,
        }))
    );
    assert!(matches!(
        newton(df, T::from_f64(10.0), &single),
        Err(SolveError::NotConverged(s)) if s.iterations == 1 && s.x < T::from_f64(10.0)
    ));

    // Ableitung null am Startwert.
    assert!(matches!(
        newton(|x: T| (x * x + one, x + x), T::ZERO, &tol),
        Err(SolveError::Singular(s)) if s.x == T::ZERO
    ));

    // Lockere Toleranz für den Restfehler beendet früher.
    let loose = tol.with_abs_f(T::from_f64(0.5));
    assert!(bisection(f, one, two, &loose).is_ok_and(|s| s.iterations == 1));
}

/// Robustheit an den Grenzen des Zahlbereichs.
fn extremes<T: Real>(huge: f64) {
    let one = T::ONE;
    let tol = Tolerance::<T>::default().with_max_iterations(4000);
    let big = T::from_f64(huge);
    let f = |x: T| x - one;
    let df = |x: T| (x - one, one);

    // Klammern, deren Breite überläuft.
    let near_one = |result: ScalarResult<T>| result.is_ok_and(|s| close(s.x, one, 64.0));
    assert!(near_one(bisection(f, -big, big, &tol)));
    assert!(near_one(brent(f, -big, big, &tol)));
    assert!(near_one(newton_bracketed(df, -big, big, &tol)));

    // Unendliche Ableitung darf nicht wie Konvergenz aussehen.
    assert!(matches!(
        newton(|_x: T| (one, T::INFINITY), one, &tol),
        Err(SolveError::NonFinite(s)) if s.residual == one
    ));
    assert!(matches!(
        newton(|_x: T| (one, T::NAN), one, &tol),
        Err(SolveError::NonFinite(_))
    ));
    assert!(matches!(
        newton(|_x: T| (T::INFINITY, one), one, &tol),
        Err(SolveError::NonFinite(_))
    ));
    // Auch nach einem winzigen Schritt gilt ein unendlicher Restfehler
    // nicht als Erfolg.
    let mut calls = 0;
    assert!(matches!(
        newton(
            |_x: T| {
                calls += 1;
                if calls == 1 { (one, T::from_f64(1.0e20)) } else { (T::INFINITY, one) }
            },
            one,
            &tol
        ),
        Err(SolveError::NonFinite(s)) if s.iterations == 1
    ));
    let mut calls = 0;
    assert!(matches!(
        secant(
            |x: T| {
                calls += 1;
                if calls <= 2 { x - one } else { T::INFINITY }
            },
            T::ZERO,
            T::from_f64(0.5),
            &tol
        ),
        Err(SolveError::NonFinite(_))
    ));
    // Mit Klammer wird stattdessen halbiert.
    assert!(near_one(newton_bracketed(|x: T| (x - one, T::INFINITY), T::ZERO, T::from_f64(3.0), &tol)));
    // Sekante mit nicht endlichem Funktionswert.
    assert!(matches!(
        secant(|x: T| if x > one { T::INFINITY } else { x - one }, T::ZERO, T::from_f64(2.0), &tol),
        Err(SolveError::NonFinite(_))
    ));
}

#[test]
fn extremes_f32() {
    extremes::<f32>(3.0e38);
}

#[test]
fn extremes_f64() {
    extremes::<f64>(1.0e308);
}

#[test]
fn scalar_f32() {
    scalar::<f32>();
}

#[test]
fn scalar_f64() {
    scalar::<f64>();
}

#[cfg(any(feature = "std", feature = "libm"))]
mod systems {
    use super::*;
    use crate::linalg::{Matrix, Vector};
    use crate::number::DualN;
    use crate::scalar::RealElementary;
    use crate::solve::{SystemResult, newton_system, newton_system_with_jacobian};

    /// x² + y² = 4 und x·y = 1.
    fn circle_and_hyperbola<T: RealElementary>(
        p: Vector<DualN<T, 2>, 2>,
    ) -> Vector<DualN<T, 2>, 2> {
        let k = |value: f64| DualN::constant(T::from_f64(value));
        Vector([p[0] * p[0] + p[1] * p[1] - k(4.0), p[0] * p[1] - k(1.0)])
    }

    fn check<T: RealElementary>() {
        let f = T::from_f64;
        let tol = Tolerance::<T>::default();
        let start = Vector([f(2.0), f(0.5)]);

        let auto = newton_system(circle_and_hyperbola::<T>, start, &tol);
        let manual = newton_system_with_jacobian(
            |p: Vector<T, 2>| {
                (
                    Vector([p[0] * p[0] + p[1] * p[1] - f(4.0), p[0] * p[1] - f(1.0)]),
                    Matrix([[p[0] + p[0], p[1] + p[1]], [p[1], p[0]]]),
                )
            },
            start,
            &tol,
        );
        // Geprüft über Rückeinsetzen.
        let solves = |result: &SystemResult<T, 2>| {
            result.as_ref().is_ok_and(|s| {
                close(s.x[0] * s.x[0] + s.x[1] * s.x[1], f(4.0), 64.0)
                    && close(s.x[0] * s.x[1], T::ONE, 64.0)
            })
        };
        assert!(solves(&auto));
        assert!(solves(&manual));
        assert!(match (&auto, &manual) {
            (Ok(a), Ok(m)) => close(a.x[0], m.x[0], 64.0) && close(a.x[1], m.x[1], 64.0),
            _ => false,
        });

        // Lineares System 2x + y = 3, x + 3y = 5 → (0,8; 1,4), praktisch in
        // einem Schritt.
        let linear = newton_system(
            |p: Vector<DualN<T, 2>, 2>| {
                let k = |value: f64| DualN::constant(T::from_f64(value));
                Vector([k(2.0) * p[0] + p[1] - k(3.0), p[0] + k(3.0) * p[1] - k(5.0)])
            },
            Vector([T::ZERO, T::ZERO]),
            &tol,
        );
        assert!(linear.is_ok_and(|s| {
            close(s.x[0], f(0.8), 64.0) && close(s.x[1], f(1.4), 64.0) && s.iterations <= 3
        }));

        // Singuläre Jacobi-Matrix: Fehler mit dem Startwert als letztem Stand.
        let origin = Vector([T::ONE, T::ONE]);
        assert!(matches!(
            newton_system_with_jacobian(
                |p: Vector<T, 2>| (
                    Vector([p[0] + p[1], p[0] + p[1]]),
                    Matrix([[T::ONE, T::ONE], [T::ONE, T::ONE]]),
                ),
                origin,
                &tol,
            ),
            Err(SolveError::Singular(s)) if s.x == origin && s.iterations == 0
        ));

        // Zu wenige Iterationen: Fehler mit letztem Stand.
        let single = tol.with_max_iterations(1);
        assert!(matches!(
            newton_system(circle_and_hyperbola::<T>, Vector([f(10.0), f(1.0)]), &single),
            Err(SolveError::NotConverged(s)) if s.iterations == 1
        ));
    }

    #[test]
    fn systems_f32() {
        check::<f32>();
    }

    #[test]
    fn systems_f64() {
        check::<f64>();
    }
}
