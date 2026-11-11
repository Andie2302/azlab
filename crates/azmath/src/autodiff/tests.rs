//! Tests der Autodiff-Hilfen. Alle Ableitungen sind von Hand bekannt.

use super::{derivative, gradient, hessian, jacobian, second_derivative, variables};
use crate::linalg::{Matrix, Vector};
use crate::scalar::Real;

fn check<T: Real>() {
    let f = T::from_f64;
    let v = |data: [f64; 2]| Vector(data.map(T::from_f64));

    // x³ bei x = 2: Wert 8, Ableitung 12, zweite Ableitung 12.
    assert_eq!(derivative(|x| x * x * x, f(2.0)), (f(8.0), f(12.0)));
    assert_eq!(second_derivative(|x| x * x * x, f(2.0)), (f(8.0), f(12.0), f(12.0)));

    // x²·y bei (3, 2): Wert 18, Gradient (2xy, x²) = (12, 9),
    // Hesse-Matrix [[2y, 2x], [2x, 0]] = [[4, 6], [6, 0]].
    assert_eq!(
        gradient(|p| p[0] * p[0] * p[1], v([3.0, 2.0])),
        (f(18.0), v([12.0, 9.0]))
    );
    assert_eq!(
        hessian(|p| p[0] * p[0] * p[1], v([3.0, 2.0])),
        (
            f(18.0),
            v([12.0, 9.0]),
            Matrix([[f(4.0), f(6.0)], [f(6.0), f(0.0)]])
        )
    );

    // (x·y, x + y) bei (3, 2): Wert (6, 5), Jacobi-Matrix [[y, x], [1, 1]].
    assert_eq!(
        jacobian(|p| Vector([p[0] * p[1], p[0] + p[1]]), v([3.0, 2.0])),
        (v([6.0, 5.0]), Matrix([[f(2.0), f(3.0)], [f(1.0), f(1.0)]]))
    );

    // Die Variablen tragen Einheitsrichtungen.
    let seeded = variables(v([3.0, 2.0]));
    assert_eq!(seeded[0].eps, [T::ONE, T::ZERO]);
    assert_eq!(seeded[1].eps, [T::ZERO, T::ONE]);
}

#[test]
fn autodiff_f32() {
    check::<f32>();
}

#[test]
fn autodiff_f64() {
    check::<f64>();
}
