//! Tests für [`CayleyDickson`].

use super::{CayleyDickson, Octonion};
use crate::number::Complex;
use crate::scalar::{Algebra, Commutative, Conjugate, Field, Real, Ring};
use crate::testutil::close;

/// Reeller Grundtyp, über dem verdoppelt wird.
trait Base: Real + Conjugate<Real = Self> {}

impl<T: Real + Conjugate<Real = T>> Base for T {}

type Cd1<T> = CayleyDickson<T>;
type Cd2<T> = CayleyDickson<Cd1<T>>;

fn cd1<T: Base>(a: f64, b: f64) -> Cd1<T> {
    Cd1::new(T::from_f64(a), T::from_f64(b))
}

fn cd2<T: Base>(c: [f64; 4]) -> Cd2<T> {
    Cd2::new(cd1(c[0], c[1]), cd1(c[2], c[3]))
}

fn oct<T: Base>(c: [f64; 8]) -> Octonion<T> {
    Octonion::new(
        cd2([c[0], c[1], c[2], c[3]]),
        cd2([c[4], c[5], c[6], c[7]]),
    )
}

fn unit<T: Base>(index: usize) -> Octonion<T> {
    let mut c = [0.0; 8];
    if let Some(slot) = c.get_mut(index) {
        *slot = 1.0;
    }
    oct(c)
}

fn flat<T: Base>(o: Octonion<T>) -> [T; 8] {
    [
        o.a.a.a, o.a.a.b, o.a.b.a, o.a.b.b, o.b.a.a, o.b.a.b, o.b.b.a, o.b.b.b,
    ]
}

fn close_oct<T: Base>(x: Octonion<T>, y: Octonion<T>, ulps: f64) -> bool {
    flat(x).into_iter().zip(flat(y)).all(|(p, q)| close(p, q, ulps))
}

// Prüfungen zur Übersetzungszeit: welche Stufe sagt was zu.
fn is_field<A: Field>() {}
fn is_commutative<A: Commutative>() {}
fn is_ring<A: Ring>() {}
fn is_algebra<A: Algebra + Conjugate>() {}

#[test]
fn trait_levels() {
    is_field::<Cd1<f64>>();
    is_commutative::<Cd1<f64>>();
    is_field::<Cd2<f64>>();
    is_ring::<Cd2<f32>>();
    is_algebra::<Octonion<f64>>();
    is_algebra::<CayleyDickson<Octonion<f32>>>();
}

fn level_one_matches_complex<T: Base>() {
    let (z1, z2) = (cd1::<T>(1.0, 2.0), cd1::<T>(3.0, 4.0));
    let (c1, c2) = (Complex::from(z1), Complex::from(z2));
    assert_eq!(z1 * z2, cd1(-5.0, 10.0));
    assert_eq!(Complex::from(z1 * z2), c1 * c2);
    assert_eq!(Complex::from(z1 + z2), c1 + c2);
    assert_eq!(Complex::from(z1.conj()), c1.conj());
    assert_eq!(z1.norm_sqr(), c1.norm_sqr());
    assert_eq!(CayleyDickson::from(c1), z1);
    let quotient = Complex::from(cd1::<T>(-5.0, 10.0) / z2);
    assert!(close(quotient.re, T::ONE, 8.0) && close(quotient.im, T::from_f64(2.0), 8.0));
}

fn level_two_is_quaternion<T: Base>() {
    let one = Cd2::<T>::ONE;
    let i = cd2::<T>([0.0, 1.0, 0.0, 0.0]);
    let j = cd2::<T>([0.0, 0.0, 1.0, 0.0]);
    let k = cd2::<T>([0.0, 0.0, 0.0, 1.0]);
    assert_eq!(i * i, -one);
    assert_eq!(j * j, -one);
    assert_eq!(k * k, -one);
    assert_eq!(i * j, k);
    assert_eq!(j * k, i);
    assert_eq!(k * i, j);
    assert_eq!(j * i, -k);
    assert_eq!(i * j * k, -one);

    let q = cd2::<T>([1.0, 2.0, 3.0, 4.0]);
    assert_eq!(q.norm_sqr(), T::from_f64(30.0));
    assert_eq!(q.conj(), cd2([1.0, -2.0, -3.0, -4.0]));
    let unity = q * q.recip();
    assert!(close(unity.a.a, T::ONE, 8.0));
    assert!(close(unity.a.b, T::ZERO, 8.0));
    assert!(close(unity.b.a, T::ZERO, 8.0));
    assert!(close(unity.b.b, T::ZERO, 8.0));
    assert_eq!(Cd2::<T>::from_f64(2.5), cd2([2.5, 0.0, 0.0, 0.0]));
    assert_eq!(q.powi(2), q * q);
}

fn octonions<T: Base>() {
    let one = Octonion::<T>::ONE;
    for index in 1..8 {
        let e = unit::<T>(index);
        assert_eq!(e * e, -one);
    }

    // Nicht assoziativ: unabhängige Einheiten anti-assoziieren.
    let (e1, e2, e4) = (unit::<T>(1), unit::<T>(2), unit::<T>(4));
    let left = (e1 * e2) * e4;
    let right = e1 * (e2 * e4);
    assert_ne!(left, right);
    assert_eq!(left, -right);

    let x = oct::<T>([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);
    let y = oct::<T>([2.0, -1.0, 0.5, 3.0, -2.0, 1.0, 4.0, -3.0]);

    // Trotzdem: Norm multiplikativ, Inverses vorhanden, alternativ.
    assert!(close((x * y).norm_sqr(), x.norm_sqr() * y.norm_sqr(), 8.0));
    assert!(close_oct(x * x.inv(), one, 16.0));
    assert!(close_oct((x / y) * y, x, 64.0));
    assert!(close_oct((x * x) * y, x * (x * y), 16.0));

    let mut acc = x;
    acc += y;
    acc -= y;
    acc *= one;
    acc /= one;
    assert!(close_oct(acc, x, 8.0));
}

#[test]
fn all_levels_f32() {
    level_one_matches_complex::<f32>();
    level_two_is_quaternion::<f32>();
    octonions::<f32>();
}

#[test]
fn all_levels_f64() {
    level_one_matches_complex::<f64>();
    level_two_is_quaternion::<f64>();
    octonions::<f64>();
}

#[cfg(any(feature = "std", feature = "libm"))]
#[test]
fn norm() {
    assert!(close(cd1::<f64>(3.0, 4.0).norm(), 5.0, 8.0));
    assert!(close(cd2::<f32>([1.0, 1.0, 1.0, 1.0]).norm(), 2.0, 8.0));
}
