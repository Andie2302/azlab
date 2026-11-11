//! Tests für [`Matrix`].

use super::Matrix;
use crate::linalg::Vector;
use crate::number::Complex;
use crate::scalar::{Algebra, Real, Ring};

fn m<T: Real, const R: usize, const C: usize>(rows: [[f64; C]; R]) -> Matrix<T, R, C> {
    Matrix(rows.map(|row| row.map(T::from_f64)))
}

fn v<T: Real, const N: usize>(data: [f64; N]) -> Vector<T, N> {
    Vector(data.map(T::from_f64))
}

fn cube<A: Ring>(x: A) -> A {
    x * x * x
}

fn algebra<T: Real>() {
    let a = m::<T, 2, 2>([[1.0, 2.0], [3.0, 4.0]]);
    let b = m::<T, 2, 2>([[5.0, 6.0], [7.0, 8.0]]);
    let two = T::from_f64(2.0);

    assert_eq!(a * b, m([[19.0, 22.0], [43.0, 50.0]]));
    assert_eq!(b * a, m([[23.0, 34.0], [31.0, 46.0]]));
    assert_eq!(a + b, m([[6.0, 8.0], [10.0, 12.0]]));
    assert_eq!(b - a, m([[4.0, 4.0], [4.0, 4.0]]));
    assert_eq!(-a, m([[-1.0, -2.0], [-3.0, -4.0]]));
    assert_eq!(a * two, m([[2.0, 4.0], [6.0, 8.0]]));
    assert_eq!(a / two, m([[0.5, 1.0], [1.5, 2.0]]));
    assert_eq!(a * v([1.0, 1.0]), v([3.0, 7.0]));
    assert_eq!(a.transpose(), m([[1.0, 3.0], [2.0, 4.0]]));
    assert_eq!(a.adjoint(), a.transpose());
    assert_eq!(a.trace(), T::from_f64(5.0));
    assert_eq!(a.diagonal(), v([1.0, 4.0]));
    assert_eq!(a.norm_sqr(), T::from_f64(30.0));

    // Rechteckig: 2×3 mal 3×2 ergibt 2×2.
    let r = m::<T, 2, 3>([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
    let s = m::<T, 3, 2>([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]);
    assert_eq!(r * s, m([[4.0, 5.0], [10.0, 11.0]]));
    assert_eq!(r.transpose(), m([[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]]));
    assert_eq!(r * v([1.0, 0.0, 1.0]), v([4.0, 10.0]));
    assert_eq!((r.rows(), r.cols()), (2, 3));

    // Einheitsmatrix und Matrizen als Ring.
    let id = Matrix::<T, 2, 2>::identity();
    assert_eq!(id, m([[1.0, 0.0], [0.0, 1.0]]));
    assert_eq!(<Matrix<T, 2, 2> as Algebra>::ONE, id);
    assert_eq!(a * id, a);
    assert_eq!(id * a, a);
    assert_eq!(cube(a), a * a * a);
    assert_eq!(Matrix::from_diagonal(v::<T, 2>([2.0, 3.0])), m([[2.0, 0.0], [0.0, 3.0]]));
    assert_eq!(Matrix::<T, 2, 3>::ZERO, m([[0.0; 3]; 2]));

    let mut acc = a;
    acc += b;
    acc -= a;
    acc *= id;
    acc *= two;
    acc /= two;
    assert_eq!(acc, b);

    // Zugriff.
    assert_eq!(a[(1, 0)], T::from_f64(3.0));
    assert_eq!(a.get(0, 1), Some(&two));
    assert_eq!(a.get(2, 0), None);
    assert_eq!(a.get(0, 2), None);
    let mut c = a;
    c[(0, 0)] = two;
    if let Some(slot) = c.get_mut(1, 1) {
        *slot = T::ZERO;
    }
    assert_eq!(c, m([[2.0, 2.0], [3.0, 0.0]]));
    assert_eq!(r.row(1), Some(v([4.0, 5.0, 6.0])));
    assert_eq!(r.row(2), None);
    assert_eq!(r.col(2), Some(v([3.0, 6.0])));
    assert_eq!(r.col(3), None);
    assert_eq!(Matrix::from_fn(|i, j| T::from_f64((i * 2 + j + 1) as f64)), a);
    assert_eq!(Matrix::from(*a.as_rows()), a);
    assert_eq!(a.map(|x| x * x), m([[1.0, 4.0], [9.0, 16.0]]));

    // Äußeres Produkt.
    assert_eq!(
        v::<T, 2>([1.0, 2.0]).outer(v::<T, 3>([3.0, 4.0, 5.0])),
        m([[3.0, 4.0, 5.0], [6.0, 8.0, 10.0]])
    );
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
    assert_eq!(2.0_f64 * Matrix([[1.0, 2.0]]), Matrix([[2.0, 4.0]]));
    assert_eq!(2.0_f32 * Matrix([[1.0], [2.0]]), Matrix([[2.0], [4.0]]));
    assert_eq!(Matrix::<f64, 1, 2>::default(), Matrix([[0.0, 0.0]]));

    // Komplex: Adjungierte konjugiert.
    let i = Complex::<f64>::I;
    let z = Matrix([[i, Complex::ONE]]);
    assert_eq!(z.adjoint(), Matrix([[-i], [Complex::ONE]]));
    assert_eq!(z.norm_sqr(), 2.0);
    assert_eq!(z * z.adjoint(), Matrix([[Complex::new(2.0, 0.0)]]));
}

#[cfg(any(feature = "std", feature = "libm"))]
#[test]
fn frobenius_norm() {
    use crate::testutil::close;
    assert!(close(m::<f64, 2, 2>([[1.0, 1.0], [1.0, 1.0]]).norm(), 2.0, 8.0));
    assert!(close(m::<f32, 1, 2>([[3.0, 4.0]]).norm(), 5.0, 8.0));
}
