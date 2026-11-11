//! Tests der Lie-Gruppen und Lie-Algebren.
//!
//! Geprüft werden ausschließlich mathematische Identitäten mit frei
//! gewählten Eingaben; es gibt keine vorab ausgerechneten Erwartungswerte.

#![cfg(any(feature = "std", feature = "libm"))]

use super::{
    LieAlgebra, LieGroup, Se2, Se2Tangent, Se3, Se3Tangent, So2, So2Tangent, So3, So3Tangent,
};
use crate::linalg::{Matrix, Vector};
use crate::number::Quaternion;
use crate::scalar::{Algebra, RealElementary};
use crate::testutil::close;

const ULPS: f64 = 512.0;

fn v<T: RealElementary, const N: usize>(data: [f64; N]) -> Vector<T, N> {
    Vector(data.map(T::from_f64))
}

fn close_v<T: RealElementary, const N: usize>(a: Vector<T, N>, b: Vector<T, N>) -> bool {
    a.0.into_iter().zip(b.0).all(|(x, y)| close(x, y, ULPS))
}

fn close_m<T: RealElementary, const R: usize, const C: usize>(
    a: Matrix<T, R, C>,
    b: Matrix<T, R, C>,
) -> bool {
    a.0.into_iter().zip(b.0).all(|(x, y)| close_v(Vector(x), Vector(y)))
}

/// Linke Seite der Jacobi-Identität; muss null ergeben.
fn jacobi<A: LieAlgebra>(x: A, y: A, z: A) -> A {
    x.bracket(y.bracket(z)) + y.bracket(z.bracket(x)) + z.bracket(x.bracket(y))
}

fn same_so3<T: RealElementary>(a: So3<T>, b: So3<T>) -> bool {
    let p = v::<T, 3>([1.0, -2.0, 0.5]);
    let q = v::<T, 3>([0.3, 0.7, -1.1]);
    close_v(a.rotate(p), b.rotate(p)) && close_v(a.rotate(q), b.rotate(q))
}

fn same_se2<T: RealElementary>(a: Se2<T>, b: Se2<T>) -> bool {
    [[0.0, 0.0], [1.5, -2.0], [-0.7, 0.4]]
        .into_iter()
        .all(|p| close_v(a.act(v(p)), b.act(v(p))))
}

fn same_se3<T: RealElementary>(a: Se3<T>, b: Se3<T>) -> bool {
    [[0.0, 0.0, 0.0], [1.0, -2.0, 0.5], [0.3, 0.7, -1.1], [-0.6, 0.2, 0.9]]
        .into_iter()
        .all(|p| close_v(a.act(v(p)), b.act(v(p))))
}

fn so2<T: RealElementary>() {
    let f = T::from_f64;
    let xi = So2Tangent { angle: f(0.7) };
    let g = So2::exp(xi);
    let h = So2::from_angle(f(-0.4));
    let k = So2::from_angle(f(1.9));
    let p = v::<T, 2>([1.5, -2.0]);

    assert!(close(g.log().angle, xi.angle, ULPS));
    assert!(close_v(So2::exp(g.log()).rotate(p), g.rotate(p)));
    assert!(close_v(g.compose(g.inverse()).rotate(p), p));
    assert!(close_v(So2::<T>::identity().rotate(p), p));
    assert!(close_v(g.compose(h).rotate(p), g.rotate(h.rotate(p))));
    assert!(close_v(
        g.compose(h).compose(k).rotate(p),
        g.compose(h.compose(k)).rotate(p)
    ));
    assert!(close(g.rotate(p).norm(), p.norm(), ULPS));
    assert!(close_v(g.matrix() * p, g.rotate(p)));
    assert_eq!(g.adjoint(xi), xi);
    assert_eq!(xi.bracket(h.log()), So2Tangent::ZERO);
    assert_eq!(xi + xi, xi.scale(f(2.0)));
    assert_eq!(xi - So2Tangent { angle: f(0.7) }, So2Tangent::ZERO);
    assert_eq!(-(-xi), xi);
}

fn so3<T: RealElementary>() {
    let f = T::from_f64;
    let xi = So3Tangent { phi: v::<T, 3>([0.3, -0.4, 0.5]) };
    let eta = So3Tangent { phi: v::<T, 3>([-0.2, 0.6, 0.1]) };
    let zeta = So3Tangent { phi: v::<T, 3>([0.7, 0.1, -0.3]) };
    let g = So3::exp(xi);
    let h = So3::exp(eta);
    let k = So3::exp(zeta);
    let p = v::<T, 3>([1.0, -2.0, 0.5]);

    assert!(close_v(g.log().phi, xi.phi));
    assert!(same_so3(So3::exp(g.log()), g));
    assert!(same_so3(g.compose(g.inverse()), So3::identity()));
    assert!(close_v(g.compose(h).rotate(p), g.rotate(h.rotate(p))));
    assert!(same_so3(g.compose(h).compose(k), g.compose(h.compose(k))));
    assert!(close(g.rotate(p).norm(), p.norm(), ULPS));
    assert!(close_v(g.matrix() * p, g.rotate(p)));
    assert!(same_so3(
        g.compose(So3::exp(eta)).compose(g.inverse()),
        So3::exp(g.adjoint(eta))
    ));

    // Lie-Klammer.
    assert!(close_v(xi.bracket(eta).phi, -eta.bracket(xi).phi));
    assert!(close_v(jacobi(xi, eta, zeta).phi, Vector::ZERO));
    assert!(close_v(xi.hat() * p, xi.phi.cross(p)));
    assert!(close_v((xi + eta - eta).phi, xi.phi));
    assert!(close_v(xi.scale(f(2.0)).phi, (xi + xi).phi));

    // Erzeugen.
    assert!(So3::from_axis_angle(xi.phi, xi.phi.norm()).is_some_and(|r| same_so3(r, g)));
    assert!(So3::from_axis_angle(Vector::<T, 3>::ZERO, f(1.0)).is_none());
    assert!(So3::from_quaternion(g.quaternion() * f(3.0)).is_some_and(|r| same_so3(r, g)));
    assert!(So3::from_quaternion(Quaternion::<T>::ZERO).is_none());
    // q und −q sind dieselbe Drehung.
    assert!(So3::from_quaternion(-g.quaternion())
        .is_some_and(|r| same_so3(r, g) && close_v(r.log().phi, xi.phi)));
}

fn se2<T: RealElementary>() {
    let f = T::from_f64;
    let xi = Se2Tangent { rho: v::<T, 2>([0.5, -1.0]), theta: f(0.8) };
    let eta = Se2Tangent { rho: v::<T, 2>([-0.3, 0.2]), theta: f(-0.5) };
    let zeta = Se2Tangent { rho: v::<T, 2>([0.9, 0.4]), theta: f(0.3) };
    let g = Se2::exp(xi);
    let h = Se2::exp(eta);
    let k = Se2::exp(zeta);
    let p = v::<T, 2>([1.5, -2.0]);

    let back = g.log();
    assert!(close_v(back.rho, xi.rho) && close(back.theta, xi.theta, ULPS));
    assert!(same_se2(Se2::exp(g.log()), g));
    assert!(same_se2(g.compose(g.inverse()), Se2::identity()));
    assert!(close_v(g.compose(h).act(p), g.act(h.act(p))));
    assert!(same_se2(g.compose(h).compose(k), g.compose(h.compose(k))));
    assert!(same_se2(
        g.compose(Se2::exp(eta)).compose(g.inverse()),
        Se2::exp(g.adjoint(eta))
    ));

    // Reine Translation: Grenzfall kleiner Winkel.
    let shift = Se2Tangent { rho: xi.rho, theta: T::ZERO };
    assert!(close_v(Se2::exp(shift).act(Vector::ZERO), xi.rho));
    assert!(close_v(Se2::exp(shift).log().rho, xi.rho));
    assert!(same_se2(Se2::new(So2::identity(), xi.rho), Se2::exp(shift)));

    // Lie-Klammer.
    let forward = xi.bracket(eta);
    let backward = -eta.bracket(xi);
    assert!(close_v(forward.rho, backward.rho) && close(forward.theta, backward.theta, ULPS));
    let cycle = jacobi(xi, eta, zeta);
    assert!(close_v(cycle.rho, Vector::ZERO) && close(cycle.theta, T::ZERO, ULPS));
}

fn se3<T: RealElementary>() {
    let xi = Se3Tangent {
        rho: v::<T, 3>([0.5, -1.0, 0.25]),
        phi: v::<T, 3>([0.3, -0.4, 0.5]),
    };
    let eta = Se3Tangent {
        rho: v::<T, 3>([-0.3, 0.2, 0.6]),
        phi: v::<T, 3>([-0.2, 0.6, 0.1]),
    };
    let zeta = Se3Tangent {
        rho: v::<T, 3>([0.9, 0.4, -0.7]),
        phi: v::<T, 3>([0.7, 0.1, -0.3]),
    };
    let g = Se3::exp(xi);
    let h = Se3::exp(eta);
    let k = Se3::exp(zeta);
    let p = v::<T, 3>([1.0, -2.0, 0.5]);

    let back = g.log();
    assert!(close_v(back.rho, xi.rho) && close_v(back.phi, xi.phi));
    assert!(same_se3(Se3::exp(g.log()), g));
    assert!(same_se3(g.compose(g.inverse()), Se3::identity()));
    assert!(close_v(g.compose(h).act(p), g.act(h.act(p))));
    assert!(same_se3(g.compose(h).compose(k), g.compose(h.compose(k))));
    assert!(same_se3(
        g.compose(Se3::exp(eta)).compose(g.inverse()),
        Se3::exp(g.adjoint(eta))
    ));

    // Reine Translation: Grenzfall kleiner Winkel.
    let shift = Se3Tangent { rho: xi.rho, phi: Vector::ZERO };
    assert!(close_v(Se3::exp(shift).act(Vector::ZERO), xi.rho));
    assert!(close_v(Se3::exp(shift).log().rho, xi.rho));
    assert!(same_se3(Se3::new(So3::identity(), xi.rho), Se3::exp(shift)));

    // Lie-Klammer.
    let forward = xi.bracket(eta);
    let backward = -eta.bracket(xi);
    assert!(close_v(forward.rho, backward.rho) && close_v(forward.phi, backward.phi));
    let cycle = jacobi(xi, eta, zeta);
    assert!(close_v(cycle.rho, Vector::ZERO) && close_v(cycle.phi, Vector::ZERO));
}

fn matrices<T: RealElementary>() {
    let xi = So3Tangent { phi: v::<T, 3>([0.3, -0.4, 0.5]) };
    let eta = So3Tangent { phi: v::<T, 3>([-0.2, 0.6, 0.1]) };
    let a = xi.hat();
    let b = eta.hat();
    let identity = Matrix::<T, 3, 3>::identity();

    assert!(close_m(Matrix::<T, 3, 3>::ZERO.exp(), identity));
    assert!(close_m(a.exp() * (-a).exp(), identity));
    // Das Matrix-Exponential der schiefsymmetrischen Matrix ist die Drehmatrix.
    assert!(close_m(a.exp(), So3::exp(xi).matrix()));
    // Große Norm: Skalieren und Quadrieren.
    let large = a * T::from_f64(8.0);
    assert!(close_m(large.exp() * (-large).exp(), identity));

    // Der Kommutator der Matrizen entspricht der Klammer der Drehvektoren.
    assert!(close_m(a.commutator(b), -b.commutator(a)));
    assert!(close_m(a.commutator(b), xi.bracket(eta).hat()));
    assert!(close_m(a.commutator(a), <Matrix<T, 3, 3> as Algebra>::ZERO));

    // Allgemeine, nicht schiefsymmetrische Matrix.
    let general = Matrix([
        [T::from_f64(0.2), T::from_f64(0.5)],
        [T::from_f64(-0.3), T::from_f64(0.1)],
    ]);
    assert!(close_m(general.exp() * (-general).exp(), Matrix::identity()));
}

#[test]
fn groups_f32() {
    so2::<f32>();
    so3::<f32>();
    se2::<f32>();
    se3::<f32>();
    matrices::<f32>();
}

#[test]
fn groups_f64() {
    so2::<f64>();
    so3::<f64>();
    se2::<f64>();
    se3::<f64>();
    matrices::<f64>();
}
