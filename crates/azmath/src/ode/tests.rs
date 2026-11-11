//! Tests der Differentialgleichungslöser.
//!
//! Verglichen wird mit der bekannten Lösung, die zur Laufzeit über die
//! Elementarfunktionen berechnet wird, und es werden Eigenschaften der
//! Verfahren geprüft. Keine vorab ausgerechneten Erwartungswerte.

#![cfg(any(feature = "std", feature = "libm"))]

use super::{
    OdeError, OdeOptions, dormand_prince, dormand_prince_observed, euler, implicit_euler, rk4,
    trapezoid,
};
use crate::linalg::Vector;
use crate::number::DualN;
use crate::scalar::RealElementary;
use crate::solve::Tolerance;

/// y' = −y mit der Lösung e^(−t).
fn decay<T: RealElementary>(_t: T, y: Vector<T, 1>) -> Vector<T, 1> {
    -y
}

/// y'' = −y als System: (y, v)' = (v, −y); Lösung (cos t, −sin t).
fn oscillator<T: RealElementary>(_t: T, y: Vector<T, 2>) -> Vector<T, 2> {
    Vector([y[1], -y[0]])
}

fn explicit<T: RealElementary>() {
    let f = T::from_f64;
    let one = Vector([T::ONE]);
    let end = f(1.0);
    let exact = (-end).exp();
    let error = |steps: u32, use_rk4: bool| {
        let result = if use_rk4 {
            rk4(decay::<T>, T::ZERO, one, end, steps)
        } else {
            euler(decay::<T>, T::ZERO, one, end, steps)
        };
        result.map(|state| (state.y[0] - exact).abs())
    };

    // Euler: 1. Ordnung, mehr Schritte verkleinern den Fehler deutlich.
    assert!(match (error(50, false), error(100, false)) {
        (Ok(coarse), Ok(fine)) => fine < coarse * f(0.75),
        _ => false,
    });
    // RK4 ist bei gleicher Schrittzahl viel genauer als Euler.
    assert!(match (error(50, false), error(50, true)) {
        (Ok(first_order), Ok(fourth_order)) => {
            fourth_order < first_order * f(0.01) && fourth_order <= T::EPSILON.sqrt()
        }
        _ => false,
    });
    // Zähler und Endzeit.
    assert!(rk4(decay::<T>, T::ZERO, one, end, 7)
        .is_ok_and(|state| state.steps == 7 && state.rejected == 0 && state.t == end));

    // Schwingung gegen cos und sin.
    let start = Vector([T::ONE, T::ZERO]);
    let limit = T::EPSILON.sqrt();
    assert!(rk4(oscillator::<T>, T::ZERO, start, f(2.0), 200).is_ok_and(|state| {
        (state.y[0] - f(2.0).cos()).abs() <= limit && (state.y[1] + f(2.0).sin()).abs() <= limit
    }));
    // Rückwärts integrieren führt zum Startwert zurück.
    assert!(rk4(oscillator::<T>, T::ZERO, start, f(2.0), 200)
        .and_then(|forward| rk4(oscillator::<T>, forward.t, forward.y, T::ZERO, 200))
        .is_ok_and(|back| {
            (back.y[0] - start[0]).abs() <= limit && (back.y[1] - start[1]).abs() <= limit
        }));

    // Fehlerfälle.
    assert_eq!(euler(decay::<T>, T::ZERO, one, end, 0), Err(OdeError::InvalidInput));
    assert_eq!(rk4(decay::<T>, T::NAN, one, end, 10), Err(OdeError::InvalidInput));
    assert_eq!(rk4(decay::<T>, T::ZERO, one, T::INFINITY, 10), Err(OdeError::InvalidInput));
    // Lösung läuft ins Unendliche: y' = y² mit Pol bei t = 1.
    assert!(matches!(
        euler(|_t: T, y: Vector<T, 1>| Vector([y[0] * y[0]]), T::ZERO, Vector([f(1.0e30)]), f(1.0e30), 4),
        Err(OdeError::NonFinite(_))
    ));
}

fn adaptive<T: RealElementary>() {
    let f = T::from_f64;
    let options = OdeOptions::<T>::default();
    let one = Vector([T::ONE]);
    let end = f(3.0);
    let limit = options.rel_tol * f(100.0);

    assert!(dormand_prince(decay::<T>, T::ZERO, one, end, &options).is_ok_and(|state| {
        state.t == end && state.steps > 0 && (state.y[0] - (-end).exp()).abs() <= limit
    }));

    let start = Vector([T::ONE, T::ZERO]);
    let forward = dormand_prince(oscillator::<T>, T::ZERO, start, f(2.0), &options);
    assert!(forward.is_ok_and(|state| {
        (state.y[0] - f(2.0).cos()).abs() <= limit && (state.y[1] + f(2.0).sin()).abs() <= limit
    }));
    // Rückwärts zum Startwert.
    assert!(forward
        .and_then(|state| dormand_prince(oscillator::<T>, state.t, state.y, T::ZERO, &options))
        .is_ok_and(|back| {
            back.t == T::ZERO
                && (back.y[0] - start[0]).abs() <= limit
                && (back.y[1] - start[1]).abs() <= limit
        }));

    // Beobachter: einmal je angenommenem Schritt, Zeiten wachsen, Ende ist das Ziel.
    let mut calls = 0_u32;
    let mut previous = T::ZERO;
    let mut ordered = true;
    let observed = dormand_prince_observed(decay::<T>, T::ZERO, one, end, &options, |t, _| {
        calls += 1;
        ordered = ordered && t > previous;
        previous = t;
    });
    assert!(observed.is_ok_and(|state| state.steps == calls && ordered && previous == end));

    // Strengere Toleranz braucht mehr Schritte.
    let loose = options.with_rel_tol(f(0.01)).with_abs_tol(f(0.01));
    assert!(
        match (
            dormand_prince(decay::<T>, T::ZERO, one, end, &loose),
            dormand_prince(decay::<T>, T::ZERO, one, end, &options),
        ) {
            (Ok(few), Ok(many)) => few.steps < many.steps,
            _ => false,
        }
    );

    // Leeres Intervall, erste Schrittweite, Fehlerfälle.
    assert!(dormand_prince(decay::<T>, end, one, end, &options)
        .is_ok_and(|state| state.steps == 0 && state.y == one));
    assert!(dormand_prince(decay::<T>, T::ZERO, one, end, &options.with_initial_step(f(0.5))).is_ok());
    assert_eq!(
        dormand_prince(decay::<T>, T::ZERO, one, T::NAN, &options),
        Err(OdeError::InvalidInput)
    );
    assert_eq!(
        dormand_prince(decay::<T>, T::ZERO, one, end, &options.with_rel_tol(-T::ONE)),
        Err(OdeError::InvalidInput)
    );
    assert_eq!(
        dormand_prince(decay::<T>, T::ZERO, one, end, &options.with_initial_step(T::ZERO)),
        Err(OdeError::InvalidInput)
    );
    // Zu wenige erlaubte Schritte: Der Fehler trägt den letzten Zustand.
    assert!(matches!(
        dormand_prince(decay::<T>, T::ZERO, one, end, &options.with_max_steps(1)),
        Err(OdeError::TooManySteps(state)) if state.t < end && state.steps + state.rejected == 1
    ));
    // Mindestschrittweite größer als das Intervall.
    assert!(matches!(
        dormand_prince(decay::<T>, T::ZERO, one, end, &options.with_min_step(f(10.0))),
        Err(OdeError::StepTooSmall(_))
    ));
}

fn implicit<T: RealElementary>() {
    let f = T::from_f64;
    let tol = Tolerance::<T>::default();
    let one = Vector([T::ONE]);
    let end = f(1.0);

    // Steif: y' = −1000·y mit nur zehn Schritten.
    let rate = f(1000.0);
    let stiff = |_t: T, y: Vector<DualN<T, 1>, 1>| y * DualN::constant(-rate);
    assert!(implicit_euler(stiff, T::ZERO, one, end, 10, &tol)
        .is_ok_and(|state| state.y[0] >= T::ZERO && state.y[0] < f(1.0e-6) && state.steps == 10));
    assert!(trapezoid(stiff, T::ZERO, one, end, 10, &tol)
        .is_ok_and(|state| state.y[0].abs() < T::ONE));
    // Das explizite Verfahren schaukelt sich bei derselben Schrittweite auf.
    assert!(match euler(|_t: T, y: Vector<T, 1>| y * (-rate), T::ZERO, one, end, 10) {
        Ok(state) => state.y[0].abs() > T::ONE,
        Err(OdeError::NonFinite(_)) => true,
        Err(_) => false,
    });

    // Genauigkeit an y' = −y: Die Trapezregel schlägt das implizite Euler-Verfahren.
    let slow = |_t: T, y: Vector<DualN<T, 1>, 1>| -y;
    let exact = (-end).exp();
    assert!(
        match (
            implicit_euler(slow, T::ZERO, one, end, 20, &tol),
            trapezoid(slow, T::ZERO, one, end, 20, &tol),
        ) {
            (Ok(first), Ok(second)) => {
                (second.y[0] - exact).abs() < (first.y[0] - exact).abs() * f(0.25)
            }
            _ => false,
        }
    );

    // Nichtlineares System: Schwingung, verglichen mit cos und sin.
    let wave = |_t: T, y: Vector<DualN<T, 2>, 2>| Vector([y[1], -y[0]]);
    assert!(trapezoid(wave, T::ZERO, Vector([T::ONE, T::ZERO]), end, 400, &tol).is_ok_and(|state| {
        (state.y[0] - end.cos()).abs() <= f(1.0e-3) && (state.y[1] + end.sin()).abs() <= f(1.0e-3)
    }));

    // Fehlerfälle.
    assert_eq!(implicit_euler(slow, T::ZERO, one, end, 0, &tol), Err(OdeError::InvalidInput));
    assert_eq!(trapezoid(slow, T::NAN, one, end, 5, &tol), Err(OdeError::InvalidInput));
}

#[test]
fn explicit_f32() {
    explicit::<f32>();
}

#[test]
fn explicit_f64() {
    explicit::<f64>();
}

#[test]
fn adaptive_f32() {
    adaptive::<f32>();
}

#[test]
fn adaptive_f64() {
    adaptive::<f64>();
}

#[test]
fn implicit_f32() {
    implicit::<f32>();
}

#[test]
fn implicit_f64() {
    implicit::<f64>();
}
