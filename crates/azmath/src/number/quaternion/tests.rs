//! Tests für [`Quaternion`].

use super::Quaternion;
use crate::number::{CayleyDickson, Complex};
use crate::scalar::{Algebra, Conjugate, Field, Real};
use crate::testutil::close;

fn q<T: Real>(w: f64, x: f64, y: f64, z: f64) -> Quaternion<T> {
    Quaternion::new(T::from_f64(w), T::from_f64(x), T::from_f64(y), T::from_f64(z))
}

fn close_q<T: Real>(a: Quaternion<T>, b: Quaternion<T>, ulps: f64) -> bool {
    close(a.w, b.w, ulps) && close(a.x, b.x, ulps) && close(a.y, b.y, ulps) && close(a.z, b.z, ulps)
}

fn algebra<T: Real + Conjugate<Real = T>>() {
    let one = Quaternion::<T>::ONE;
    let (i, j, k) = (Quaternion::<T>::I, Quaternion::<T>::J, Quaternion::<T>::K);
    assert_eq!(i * i, -one);
    assert_eq!(j * j, -one);
    assert_eq!(k * k, -one);
    assert_eq!(i * j * k, -one);
    assert_eq!(i * j, k);
    assert_eq!(j * k, i);
    assert_eq!(k * i, j);
    assert_eq!(j * i, -k);

    let a = q::<T>(1.0, 2.0, 3.0, 4.0);
    let b = q::<T>(5.0, 6.0, 7.0, 8.0);
    assert_eq!(a * b, q(-60.0, 12.0, 30.0, 24.0));
    assert_eq!(a + b, q(6.0, 8.0, 10.0, 12.0));
    assert_eq!(a - b, q(-4.0, -4.0, -4.0, -4.0));
    assert_eq!(a.conj(), q(1.0, -2.0, -3.0, -4.0));
    assert_eq!(a.norm_sqr(), T::from_f64(30.0));
    assert_eq!(a.dot(b), T::from_f64(70.0));
    assert_eq!(a * T::from_f64(2.0), q(2.0, 4.0, 6.0, 8.0));
    assert_eq!(a / T::from_f64(2.0), q(0.5, 1.0, 1.5, 2.0));
    assert_eq!(a.scalar(), T::ONE);
    assert_eq!(Quaternion::from_parts(a.scalar(), a.vector()), a);
    assert_eq!(Quaternion::<T>::from_f64(2.5), q(2.5, 0.0, 0.0, 0.0));
    assert_eq!(Quaternion::from(T::ONE), one);
    assert_eq!(Quaternion::from(Complex::new(T::ONE, T::ONE)), q(1.0, 1.0, 0.0, 0.0));

    assert!(close_q(a * a.recip(), one, 8.0));
    assert!(close_q(a.recip() * a, one, 8.0));
    assert!(close_q((a / b) * b, a, 16.0));
    assert_eq!(a.powi(2), a * a);

    let mut acc = a;
    acc += b;
    acc -= b;
    acc *= one;
    acc /= one;
    assert!(close_q(acc, a, 8.0));

    // Übereinstimmung mit der Cayley-Dickson-Konstruktion.
    let (ca, cb) = (
        CayleyDickson::<CayleyDickson<T>>::from(a),
        CayleyDickson::<CayleyDickson<T>>::from(b),
    );
    assert_eq!(Quaternion::from(ca * cb), a * b);
    assert_eq!(Quaternion::from(ca.conj()), a.conj());
    assert_eq!(Quaternion::from(ca), a);
    let (za, zb) = a.to_complex_pair();
    assert_eq!((za, zb), (Complex::new(T::ONE, T::from_f64(2.0)), Complex::new(T::from_f64(3.0), T::from_f64(4.0))));
    assert_eq!(Quaternion::from_complex_pair(za, zb), a);
}

#[test]
fn algebra_f32() {
    algebra::<f32>();
}

#[test]
fn algebra_f64() {
    algebra::<f64>();
}

#[test]
fn left_scalar() {
    assert_eq!(2.0_f64 * Quaternion::new(1.0, 2.0, 3.0, 4.0), Quaternion::new(2.0, 4.0, 6.0, 8.0));
    assert_eq!(2.0_f32 * Quaternion::new(1.0, 2.0, 3.0, 4.0), Quaternion::new(2.0, 4.0, 6.0, 8.0));
}

#[cfg(any(feature = "std", feature = "libm"))]
mod rotation {
    use super::*;
    use crate::scalar::RealElementary;

    const X: [f64; 3] = [1.0, 0.0, 0.0];
    const Z: [f64; 3] = [0.0, 0.0, 1.0];

    fn vec3<T: Real>(v: [f64; 3]) -> [T; 3] {
        v.map(T::from_f64)
    }

    fn close_v<T: Real>(a: [T; 3], b: [T; 3], ulps: f64) -> bool {
        a.into_iter().zip(b).all(|(p, r)| close(p, r, ulps))
    }

    /// Drehung um `axis`; die Achsen in diesen Tests sind nie null.
    fn rot<T: RealElementary>(axis: [f64; 3], angle: T) -> Quaternion<T> {
        Quaternion::from_axis_angle(vec3(axis), angle).unwrap_or(Quaternion::IDENTITY)
    }

    fn check<T: RealElementary>() {
        let one = Quaternion::<T>::IDENTITY;
        let quarter = T::PI * T::from_f64(0.5);
        let eighth = T::PI * T::from_f64(0.25);
        let rz90 = rot::<T>(Z, quarter);
        let rz45 = rot::<T>(Z, eighth);
        let rx90 = rot::<T>(X, quarter);

        // Norm und Normieren.
        assert!(close(q::<T>(1.0, 1.0, 1.0, 1.0).norm(), T::from_f64(2.0), 8.0));
        assert!(Quaternion::<T>::ZERO.normalize().is_none());
        assert!(
            q::<T>(0.0, 3.0, 0.0, 4.0)
                .normalize()
                .is_some_and(|n| close_q(n, q(0.0, 0.6, 0.0, 0.8), 8.0))
        );

        // Achse und Winkel.
        assert!(Quaternion::from_axis_angle(vec3::<T>([0.0, 0.0, 2.0]), quarter)
            .is_some_and(|r| close_q(r, rz90, 8.0)));
        assert!(Quaternion::from_axis_angle(vec3::<T>([0.0; 3]), quarter).is_none());
        let (axis, angle) = rz90.to_axis_angle();
        assert!(close_v(axis, vec3(Z), 8.0) && close(angle, quarter, 8.0));
        let (axis, angle) = one.to_axis_angle();
        assert!(close_v(axis, vec3(X), 8.0) && angle == T::ZERO);

        // Drehen: 90° um z bildet x auf y ab.
        assert!(close_v(rz90.rotate(vec3(X)), vec3([0.0, 1.0, 0.0]), 8.0));
        assert!(close_v(rx90.rotate(vec3([0.0, 1.0, 0.0])), vec3(Z), 8.0));
        let v = vec3::<T>([1.0, 2.0, 3.0]);
        // Hintereinanderausführung = Produkt.
        assert!(close_v((rz90 * rx90).rotate(v), rz90.rotate(rx90.rotate(v)), 16.0));
        // Übereinstimmung mit q·v·q*.
        let sandwich = rz45 * Quaternion::from_parts(T::ZERO, v) * rz45.conj();
        assert!(close_v(rz45.rotate(v), sandwich.vector(), 16.0));

        // Interpolation.
        assert!(close_q(one.slerp(rz90, T::ZERO), one, 8.0));
        assert!(close_q(one.slerp(rz90, T::ONE), rz90, 8.0));
        assert!(close_q(one.slerp(rz90, T::from_f64(0.5)), rz45, 8.0));
        assert!(close_q(one.slerp(-rz90, T::ONE), rz90, 8.0));
        assert!(close_q(rz45.slerp(rz45, T::from_f64(0.3)), rz45, 8.0));
        assert!(close_q(one.nlerp(rz90, T::from_f64(0.5)), rz45, 8.0));
        assert!(close_q(one.nlerp(-one, T::from_f64(0.5)), one, 8.0));

        // exp und ln.
        let small = q::<T>(0.5, 0.1, 0.2, 0.3);
        assert!(close_q(small.exp().ln(), small, 32.0));
        let p = q::<T>(1.0, 2.0, 3.0, 4.0);
        assert!(close_q(p.ln().exp(), p, 32.0));
        assert!(close_q(one.exp(), Quaternion::from_real(T::E), 8.0));
        assert!(close_q((-one).ln(), q(0.0, 1.0, 0.0, 0.0) * T::PI, 8.0));
        assert!(close_q((Quaternion::<T>::I * T::PI).exp(), -one, 8.0));
        assert!(close_q(Quaternion::from_real(T::E).ln(), one, 8.0));

        // Potenzen.
        assert!(close_q(rz90.powf(T::from_f64(0.5)), rz45, 16.0));
        assert!(close_q(rz90.sqrt(), rz45, 16.0));
        assert!(close_q(p.powf(T::from_f64(2.0)), p * p, 64.0));
        assert_eq!(Quaternion::<T>::ZERO.sqrt(), Quaternion::ZERO);
    }

    #[test]
    fn rotation_f32() {
        check::<f32>();
    }

    #[test]
    fn rotation_f64() {
        check::<f64>();
    }
}
