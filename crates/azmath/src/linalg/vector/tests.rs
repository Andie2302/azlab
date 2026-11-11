//! Tests für [`Vector`].

use super::Vector;
use crate::number::Complex;
use crate::scalar::{Algebra, Real};

fn v<T: Real, const N: usize>(data: [f64; N]) -> Vector<T, N> {
    Vector(data.map(T::from_f64))
}

fn algebra<T: Real>() {
    let a = v::<T, 3>([1.0, 2.0, 3.0]);
    let b = v::<T, 3>([4.0, 5.0, 6.0]);
    let two = T::from_f64(2.0);

    assert_eq!(a + b, v([5.0, 7.0, 9.0]));
    assert_eq!(b - a, v([3.0, 3.0, 3.0]));
    assert_eq!(-a, v([-1.0, -2.0, -3.0]));
    assert_eq!(a * two, v([2.0, 4.0, 6.0]));
    assert_eq!(b / two, v([2.0, 2.5, 3.0]));
    assert_eq!(a.dot(b), T::from_f64(32.0));
    assert_eq!(a.inner(b), T::from_f64(32.0));
    assert_eq!(a.cross(b), v([-3.0, 6.0, -3.0]));
    assert_eq!(v::<T, 3>([1.0, 0.0, 0.0]).cross(v([0.0, 1.0, 0.0])), v([0.0, 0.0, 1.0]));
    assert_eq!(a.component_mul(b), v([4.0, 10.0, 18.0]));
    assert_eq!(a.sum(), T::from_f64(6.0));
    assert_eq!(a.norm_sqr(), T::from_f64(14.0));

    let mut acc = a;
    acc += b;
    acc -= a;
    acc *= two;
    acc /= two;
    assert_eq!(acc, b);

    // Zugriff.
    assert_eq!(a[1], two);
    assert_eq!(a.get(2), Some(&T::from_f64(3.0)));
    assert_eq!(a.get(3), None);
    let mut c = a;
    c[0] = two;
    if let Some(slot) = c.get_mut(1) {
        *slot = T::ZERO;
    }
    assert_eq!(c, v([2.0, 0.0, 3.0]));
    assert_eq!(a.iter().count(), 3);
    assert_eq!((&a).into_iter().copied().fold(T::ZERO, |s, x| s + x), a.sum());
    assert_eq!(a.as_slice().len(), 3);
    assert_eq!(<[T; 3]>::from(a), *a.as_array());
    assert_eq!(Vector::from([T::ONE, T::ZERO]), v([1.0, 0.0]));

    // Konstruktoren.
    assert_eq!(Vector::<T, 3>::ZERO, v([0.0, 0.0, 0.0]));
    assert_eq!(Vector::<T, 2>::splat(two), v([2.0, 2.0]));
    assert_eq!(Vector::<T, 3>::unit(1), Some(v([0.0, 1.0, 0.0])));
    assert_eq!(Vector::<T, 3>::unit(3), None);
    assert_eq!(a.map(|x| x * x), v([1.0, 4.0, 9.0]));
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
fn concrete_types() {
    assert_eq!(2.0_f64 * Vector([1.0, 2.0]), Vector([2.0, 4.0]));
    assert_eq!(2.0_f32 * Vector([1.0, 2.0]), Vector([2.0, 4.0]));
    assert_eq!(Vector::<f64, 2>::default(), Vector([0.0, 0.0]));

    // Komplex: hermitesches und bilineares Produkt unterscheiden sich.
    let i = Complex::<f64>::I;
    let z = Vector([i, Complex::ONE]);
    assert_eq!(z.dot(z), Complex::ZERO);
    assert_eq!(z.inner(z), Complex::new(2.0, 0.0));
    assert_eq!(z.norm_sqr(), 2.0);
}

#[cfg(any(feature = "std", feature = "libm"))]
mod norms {
    use super::*;
    use crate::scalar::RealElementary;
    use crate::testutil::close;

    fn check<T: RealElementary>() {
        let a = v::<T, 2>([3.0, 4.0]);
        assert!(close(a.norm(), T::from_f64(5.0), 8.0));
        assert!(close(a.distance(v([0.0, 0.0])), T::from_f64(5.0), 8.0));
        assert!(a.normalize().is_some_and(|n| close(n[0], T::from_f64(0.6), 8.0)
            && close(n[1], T::from_f64(0.8), 8.0)));
        assert!(Vector::<T, 2>::ZERO.normalize().is_none());

        let z = Vector([Complex::new(T::ZERO, T::from_f64(3.0)), Complex::from_real(T::from_f64(4.0))]);
        assert!(close(z.norm(), T::from_f64(5.0), 8.0));
        assert!(z.normalize().is_some_and(|n| close(n.norm(), T::ONE, 8.0)));
    }

    /// Normen, deren Quadratsumme über- oder unterläuft.
    fn extremes<T: RealElementary>(tiny: f64, huge: f64) {
        use crate::linalg::Matrix;
        use crate::number::Quaternion;
        let root2 = T::from_f64(2.0).sqrt();
        for value in [T::from_f64(tiny), T::from_f64(huge)] {
            assert!(close(Vector([value, T::ZERO]).norm() / value, T::ONE, 8.0));
            assert!(close(Vector([value, -value]).norm() / value, root2, 8.0));
            assert!(Vector([value, value])
                .normalize()
                .is_some_and(|n| close(n[0] * root2, T::ONE, 8.0)));
            assert!(close(Matrix([[value, value], [value, value]]).norm() / value, T::from_f64(2.0), 8.0));
            let z = Vector([Complex::new(value, value)]);
            assert!(close(z.norm() / value, root2, 8.0));
            assert!(close(Quaternion::new(value, value, value, value).norm() / value, T::from_f64(2.0), 8.0));
            assert!(Quaternion::new(T::ZERO, value, T::ZERO, T::ZERO)
                .normalize()
                .is_some_and(|q| close(q.x, T::ONE, 8.0)));
        }
        // Sonderwerte.
        assert_eq!(Vector::<T, 3>::ZERO.norm(), T::ZERO);
        assert_eq!(Vector::<T, 0>::ZERO.norm(), T::ZERO);
        assert_eq!(Vector([T::INFINITY, T::ONE]).norm(), T::INFINITY);
        assert!(Vector([T::NAN, T::ONE]).norm().is_nan());
        assert!(Vector([T::NAN, T::ZERO]).norm().is_nan());
        assert!(Vector([T::INFINITY, T::ONE]).normalize().is_none());
    }

    #[test]
    fn extremes_f32() {
        extremes::<f32>(1.0e-30, 1.0e30);
    }

    #[test]
    fn extremes_f64() {
        extremes::<f64>(1.0e-200, 1.0e200);
    }

    #[test]
    fn norms_f32() {
        check::<f32>();
    }

    #[test]
    fn norms_f64() {
        check::<f64>();
    }
}
