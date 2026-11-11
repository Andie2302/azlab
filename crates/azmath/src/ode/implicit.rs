//! Implizite Verfahren mit fester Schrittweite für steife Probleme.

use super::error::{OdeError, OdeResult};
use super::fixed::fixed_steps;
use crate::linalg::Vector;
use crate::number::DualN;
use crate::scalar::RealElementary;
use crate::solve::{Tolerance, newton_system};

/// Gemeinsame Umsetzung: Jeder Schritt löst
/// `z = y + h·((1 − theta)·f(t, y) + theta·f(t + h, z))` nach `z`.
/// `theta = 1` ist das implizite Euler-Verfahren, `theta = 1/2` die
/// Trapezregel.
fn theta_method<T, F, const N: usize>(
    mut f: F,
    t0: T,
    y0: Vector<T, N>,
    t1: T,
    steps: u32,
    tol: &Tolerance<T>,
    theta: T,
) -> OdeResult<T, N>
where
    T: RealElementary,
    F: FnMut(T, Vector<DualN<T, N>, N>) -> Vector<DualN<T, N>, N>,
{
    let explicit_weight = T::ONE - theta;
    fixed_steps(t0, y0, t1, steps, |t, y, h, state| {
        let lifted = y.map(DualN::constant);
        // Bekannter Anteil der rechten Seite.
        let known = if explicit_weight == T::ZERO {
            lifted
        } else {
            lifted + f(t, lifted) * DualN::constant(h * explicit_weight)
        };
        let implicit_step = DualN::constant(h * theta);
        let end = t + h;
        match newton_system(|z| z - known - f(end, z) * implicit_step, y, tol) {
            Ok(solution) => Ok(solution.x),
            Err(_) => Err(OdeError::SolverFailed(state)),
        }
    })
}

/// Implizites Euler-Verfahren von `t0` nach `t1` in `steps` Schritten.
///
/// 1. Ordnung, aber für jede Schrittweite stabil: das robuste Verfahren für
/// steife Probleme. `f` wird mit Dualzahlen aufgerufen, damit die
/// Jacobi-Matrix für das Newton-Verfahren automatisch entsteht; `tol`
/// steuert das Newton-Verfahren.
pub fn implicit_euler<T, F, const N: usize>(
    f: F,
    t0: T,
    y0: Vector<T, N>,
    t1: T,
    steps: u32,
    tol: &Tolerance<T>,
) -> OdeResult<T, N>
where
    T: RealElementary,
    F: FnMut(T, Vector<DualN<T, N>, N>) -> Vector<DualN<T, N>, N>,
{
    theta_method(f, t0, y0, t1, steps, tol, T::ONE)
}

/// Trapezregel (Crank-Nicolson) von `t0` nach `t1` in `steps` Schritten.
///
/// 2. Ordnung und stabil für steife Probleme, dämpft sehr schnelle Anteile
/// aber weniger stark als [`implicit_euler`]. Aufruf wie dort.
pub fn trapezoid<T, F, const N: usize>(
    f: F,
    t0: T,
    y0: Vector<T, N>,
    t1: T,
    steps: u32,
    tol: &Tolerance<T>,
) -> OdeResult<T, N>
where
    T: RealElementary,
    F: FnMut(T, Vector<DualN<T, N>, N>) -> Vector<DualN<T, N>, N>,
{
    theta_method(f, t0, y0, t1, steps, tol, T::from_f64(0.5))
}
