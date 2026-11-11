//! Tests der Optimierer an Standardproblemen mit bekanntem Minimum.

#![cfg(any(feature = "std", feature = "libm"))]

use super::{
    bfgs, gradient_descent, lbfgs, levenberg_marquardt, levenberg_marquardt_auto, nelder_mead,
    simplex_around,
};
use crate::autodiff::{gradient, reverse_gradient};
use crate::linalg::{Matrix, Vector};
use crate::number::DualN;
use crate::scalar::{Algebra, Elementary, Field, RealElementary};
use crate::solve::{SolveError, SystemResult, Tolerance};

/// Rosenbrock-Funktion: schmales, gekrümmtes Tal, Minimum 0 bei (1, 1).
fn rosenbrock<S: Field>(p: Vector<S, 2>) -> S {
    let [x, y] = p.0;
    let (a, b) = (S::ONE - x, y - x * x);
    a * a + S::from_f64(100.0) * b * b
}

/// Quadratische Schale mit Minimum bei (1, −2, 3).
fn bowl<S: Field>(p: Vector<S, 3>) -> S {
    let [x, y, z] = p.0;
    let (a, b, c) = (x - S::ONE, y + S::from_f64(2.0), z - S::from_f64(3.0));
    a * a + S::from_f64(2.0) * b * b + S::from_f64(3.0) * c * c + a * b
}

fn near<T: RealElementary, const N: usize>(
    result: SystemResult<T, N>,
    expected: [f64; N],
    limit: T,
) -> bool {
    result.is_ok_and(|s| {
        s.x.iter().zip(expected).all(|(&got, want)| (got - T::from_f64(want)).abs() <= limit)
    })
}

fn gradient_based<T: RealElementary + 'static>() {
    let f = T::from_f64;
    let tol = Tolerance::<T>::default().with_max_iterations(20_000);
    let limit = T::EPSILON.sqrt() * f(100.0);
    let start = Vector([f(-1.2), f(1.0)]);
    let valley = |x: Vector<T, 2>| gradient(rosenbrock, x);

    assert!(near(bfgs(valley, start, &tol), [1.0, 1.0], limit));
    assert!(near(lbfgs::<T, _, 2, 5>(valley, start, &tol), [1.0, 1.0], limit));
    assert!(near(lbfgs::<T, _, 2, 1>(valley, start, &tol), [1.0, 1.0], limit));
    // Gradient aus dem Rückwärtsmodus liefert dasselbe.
    let reverse = |x: Vector<T, 2>| {
        reverse_gradient::<T, _, 2, 64>(|p| rosenbrock(p), x).unwrap_or((T::NAN, Vector::ZERO))
    };
    assert!(near(bfgs(reverse, start, &tol), [1.0, 1.0], limit));
    // BFGS braucht viel weniger Schritte als der Gradientenabstieg.
    let quick = bfgs(valley, start, &tol);
    let slow = gradient_descent(valley, start, &tol.with_abs_f(f(1.0e-3)));
    assert!(match (quick, slow) {
        (Ok(fast), Ok(steady)) => fast.iterations < 200 && steady.iterations > fast.iterations,
        _ => false,
    });

    // Gut konditioniert: Alle Verfahren finden das Minimum der Schale.
    let origin = Vector([T::ZERO; 3]);
    let shell = |x: Vector<T, 3>| gradient(bowl, x);
    // Minimum: 2a + b = 0 und 4b + a = 0 → a = b = 0, also (1, −2, 3).
    for result in [
        gradient_descent(shell, origin, &tol),
        bfgs(shell, origin, &tol),
        lbfgs::<T, _, 3, 4>(shell, origin, &tol),
    ] {
        assert!(near(result, [1.0, -2.0, 3.0], limit));
        assert!(result.is_ok_and(|s| s.residual <= limit));
    }
    // Start im Minimum: sofort fertig.
    let minimum = Vector([f(1.0), f(-2.0), f(3.0)]);
    assert!(bfgs(shell, minimum, &tol).is_ok_and(|s| s.iterations == 0 && s.residual == T::ZERO));

    // Fehlerfälle.
    assert!(matches!(
        bfgs(valley, start, &tol.with_max_iterations(2)),
        Err(SolveError::NotConverged(s)) if s.iterations == 2
    ));
    assert!(matches!(
        gradient_descent(|x: Vector<T, 2>| (T::NAN, x), start, &tol),
        Err(SolveError::NonFinite(_))
    ));
    assert_eq!(lbfgs::<T, _, 2, 0>(valley, start, &tol), Err(SolveError::InvalidInput));
}

fn derivative_free<T: RealElementary>() {
    let f = T::from_f64;
    let tol = Tolerance::<T>::default().with_max_iterations(5000);
    let limit = T::EPSILON.sqrt().sqrt() * f(10.0);

    let simplex = simplex_around::<T, 3, 4>(Vector([T::ZERO; 3]), T::ONE);
    assert!(simplex.is_some_and(|corners| {
        corners[0] == Vector::ZERO
            && corners[2] == Vector([T::ZERO, T::ONE, T::ZERO])
            && near(nelder_mead(bowl, corners, &tol), [1.0, -2.0, 3.0], limit)
    }));
    let simplex = simplex_around::<T, 2, 3>(Vector([f(-1.2), f(1.0)]), f(0.5));
    assert!(simplex.is_some_and(|corners| near(nelder_mead(rosenbrock, corners, &tol), [1.0, 1.0], limit)));
    // Nicht glatt: |x − 1| + |y + 2|.
    let kinked = |p: Vector<T, 2>| (p[0] - T::ONE).abs() + (p[1] + f(2.0)).abs();
    let simplex = simplex_around::<T, 2, 3>(Vector([T::ZERO; 2]), T::ONE);
    assert!(simplex.is_some_and(|corners| near(nelder_mead(kinked, corners, &tol), [1.0, -2.0], limit)));

    // Fehlerfälle.
    assert!(simplex_around::<T, 2, 2>(Vector([T::ZERO; 2]), T::ONE).is_none());
    assert_eq!(
        nelder_mead(rosenbrock, [Vector([T::ZERO; 2]); 2], &tol),
        Err(SolveError::InvalidInput)
    );
    let simplex = simplex_around::<T, 2, 3>(Vector([f(-1.2), f(1.0)]), f(0.5));
    assert!(simplex.is_some_and(|corners| matches!(
        nelder_mead(rosenbrock, corners, &tol.with_max_iterations(3)),
        Err(SolveError::NotConverged(s)) if s.iterations == 3
    )));
    assert!(simplex.is_some_and(|corners| matches!(
        nelder_mead(|_p: Vector<T, 2>| T::NAN, corners, &tol),
        Err(SolveError::NonFinite(_))
    )));
}

fn least_squares<T: RealElementary>() {
    let f = T::from_f64;
    let tol = Tolerance::<T>::default().with_max_iterations(500);
    let limit = T::EPSILON.sqrt() * f(100.0);

    // Messwerte des Modells y = a·e^(b·t) mit a = 2, b = −0,5.
    let times: [T; 6] = core::array::from_fn(|i| f(i as f64) * f(0.5));
    let data = times.map(|t| f(2.0) * (f(-0.5) * t).exp());
    let start = Vector([T::ONE, T::ZERO]);

    // Jacobi-Matrix automatisch …
    let auto = levenberg_marquardt_auto(
        |p: Vector<DualN<T, 2>, 2>| {
            Vector(core::array::from_fn::<_, 6, _>(|i| {
                p[0] * (p[1] * DualN::constant(times[i])).exp() - DualN::constant(data[i])
            }))
        },
        start,
        &tol,
    );
    assert!(near(auto, [2.0, -0.5], limit));
    assert!(auto.is_ok_and(|s| s.residual <= limit));
    // … und von Hand.
    let manual = levenberg_marquardt(
        |p: Vector<T, 2>| {
            let residuals = Vector(core::array::from_fn::<_, 6, _>(|i| {
                p[0] * (p[1] * times[i]).exp() - data[i]
            }));
            let jacobi = Matrix(core::array::from_fn::<_, 6, _>(|i| {
                let e = (p[1] * times[i]).exp();
                [e, p[0] * times[i] * e]
            }));
            (residuals, jacobi)
        },
        start,
        &tol,
    );
    assert!(near(manual, [2.0, -0.5], limit));

    // Lineares Problem (Ausgleichsgerade): Lösung der Normalgleichungen.
    // Punkte (0, 1), (1, 3), (2, 5), (3, 7) liegen auf y = 2t + 1.
    let line = levenberg_marquardt_auto(
        |p: Vector<DualN<T, 2>, 2>| {
            Vector(core::array::from_fn::<_, 4, _>(|i| {
                let t = DualN::constant(T::from_f64(i as f64));
                p[0] * t + p[1] - (DualN::constant(T::from_f64(2.0)) * t + DualN::ONE)
            }))
        },
        Vector([T::ZERO; 2]),
        &tol,
    );
    assert!(near(line, [2.0, 1.0], limit));

    // Fehlerfälle.
    assert!(matches!(
        levenberg_marquardt(
            |p: Vector<T, 1>| (Vector([T::NAN]), Matrix([[p[0]]])),
            Vector([T::ONE]),
            &tol
        ),
        Err(SolveError::NonFinite(_))
    ));
    assert!(matches!(
        levenberg_marquardt_auto(
            |p: Vector<DualN<T, 2>, 2>| {
                Vector(core::array::from_fn::<_, 6, _>(|i| {
                    p[0] * (p[1] * DualN::constant(times[i])).exp() - DualN::constant(data[i])
                }))
            },
            start,
            &tol.with_max_iterations(1)
        ),
        Err(SolveError::NotConverged(_))
    ));
}

#[test]
fn optim_f32() {
    gradient_based::<f32>();
    derivative_free::<f32>();
    least_squares::<f32>();
}

#[test]
fn optim_f64() {
    gradient_based::<f64>();
    derivative_free::<f64>();
    least_squares::<f64>();
}
