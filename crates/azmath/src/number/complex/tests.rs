//! Tests für [`Complex`]. Erwartungswerte sind exakt von Hand herleitbar
//! oder mathematische Identitäten.

use super::Complex;
use crate::scalar::{Algebra, Conjugate, Field, Real};
use crate::testutil::close;

fn c<T: Real>(re: f64, im: f64) -> Complex<T> {
    Complex::new(T::from_f64(re), T::from_f64(im))
}

fn close_c<T: Real>(a: Complex<T>, b: Complex<T>, ulps: f64) -> bool {
    close(a.re, b.re, ulps) && close(a.im, b.im, ulps)
}

fn algebra<T: Real>() {
    let i = Complex::<T>::I;
    assert_eq!(i * i, -Complex::ONE);
    assert_eq!(c::<T>(1.0, 2.0) * c(3.0, 4.0), c(-5.0, 10.0));
    assert_eq!(c::<T>(1.0, 2.0) + c(3.0, 4.0), c(4.0, 6.0));
    assert_eq!(c::<T>(1.0, 2.0) - c(3.0, 4.0), c(-2.0, -2.0));
    assert!(close_c(c::<T>(-5.0, 10.0) / c(3.0, 4.0), c(1.0, 2.0), 8.0));
    assert_eq!(c::<T>(1.0, 2.0).conj(), c(1.0, -2.0));
    assert_eq!(c::<T>(3.0, 4.0).norm_sqr(), T::from_f64(25.0));
    assert_eq!(c::<T>(1.0, 2.0) * T::from_f64(2.0), c(2.0, 4.0));
    assert_eq!(c::<T>(1.0, 2.0).mul_i(), c::<T>(1.0, 2.0) * i);
    assert_eq!(Complex::<T>::from_f64(2.5), c(2.5, 0.0));

    let z = c::<T>(0.5, -1.5);
    assert!(close_c(z * z.recip(), Complex::ONE, 8.0));
    assert!(close_c(z.powi(3), z * z * z, 8.0));

    let mut acc = z;
    acc += z;
    acc *= c(0.0, 1.0);
    acc -= z;
    acc /= c(2.0, 0.0);
    assert!(close_c(acc, ((z + z) * i - z) / c(2.0, 0.0), 8.0));
}

fn division_avoids_overflow<T: Real>() {
    // Das Quadrat von 1e30 läuft in f32 über; die Division muss trotzdem
    // ein endliches Ergebnis liefern.
    let quotient = c::<T>(2.0e30, 2.0e30) / c(1.0e30, 1.0e30);
    assert!(close_c(quotient, c(2.0, 0.0), 8.0));
}

#[test]
fn algebra_f32() {
    algebra::<f32>();
    division_avoids_overflow::<f32>();
}

#[test]
fn algebra_f64() {
    algebra::<f64>();
    division_avoids_overflow::<f64>();
}

#[test]
fn left_scalar_and_conversions() {
    assert_eq!(2.0_f64 * Complex::new(1.0, 2.0), Complex::new(2.0, 4.0));
    assert_eq!(2.0_f32 * Complex::new(1.0, 2.0), Complex::new(2.0, 4.0));
    assert_eq!(Complex::from(3.0_f64), Complex::new(3.0, 0.0));
    assert_eq!(Complex::from((1.0_f32, 2.0)), Complex::new(1.0, 2.0));
    assert_eq!(<(f64, f64)>::from(Complex::new(1.0, 2.0)), (1.0, 2.0));
}

#[cfg(any(feature = "std", feature = "libm"))]
mod elementary {
    use super::*;
    use crate::scalar::{Elementary, RealElementary};

    fn check<T: RealElementary>() {
        let one = Complex::<T>::ONE;
        let i = Complex::<T>::I;

        // Eulersche Identität und Wurzel aus −1.
        assert!(close_c(Complex::new(T::ZERO, T::PI).exp(), -one, 8.0));
        assert_eq!(c::<T>(-1.0, 0.0).sqrt(), i);
        // Unterhalb des Schnitts (Imaginärteil −0) liegt der andere Zweig.
        assert_eq!(c::<T>(-1.0, -0.0).sqrt(), -i);
        assert!(close_c(c::<T>(3.0, 4.0).sqrt(), c(2.0, 1.0), 8.0));
        assert!(close_c(c::<T>(3.0, -4.0).sqrt(), c(2.0, -1.0), 8.0));
        assert!(close_c(c::<T>(-3.0, 4.0).sqrt(), c(1.0, 2.0), 8.0));
        assert_eq!(Complex::<T>::ZERO.sqrt(), Complex::ZERO);

        // Betrag, Winkel, Polarform.
        assert!(close(c::<T>(3.0, 4.0).norm(), T::from_f64(5.0), 8.0));
        assert!(close(i.arg(), T::PI * T::from_f64(0.5), 8.0));
        let (r, phi) = c::<T>(-1.0, 1.0).to_polar();
        assert!(close_c(Complex::from_polar(r, phi), c(-1.0, 1.0), 8.0));

        // Umkehrpaare.
        let z = c::<T>(0.5, 1.0);
        assert!(close_c(z.exp().ln(), z, 16.0));
        assert!(close_c(z.exp2().log2(), z, 16.0));
        assert!(close_c(z.ln().exp(), z, 16.0));
        assert!(close_c(c::<T>(100.0, 0.0).log10(), c(2.0, 0.0), 16.0));
        assert!(close_c(z.exp_m1(), z.exp() - one, 16.0));
        assert!(close_c(z.ln_1p(), (one + z).ln(), 16.0));
        assert!(close_c(z.powf(c(2.0, 0.0)), z * z, 16.0));
        assert_eq!(Complex::<T>::ZERO.powf(z), Complex::ZERO);
        assert_eq!(Complex::<T>::ZERO.powf(Complex::ZERO), one);

        // Trigonometrie und Hyperbelfunktionen.
        let w = c::<T>(0.3, 0.2);
        let (s, co) = (w.sin(), w.cos());
        assert!(close_c(s * s + co * co, one, 16.0));
        assert!(close_c(w.tan(), s / co, 16.0));
        assert!(close_c(w.tanh(), w.sinh() / w.cosh(), 16.0));
        assert!(close_c(w.cosh() * w.cosh() - w.sinh() * w.sinh(), one, 16.0));
        assert!(close_c(w.sin().asin(), w, 64.0));
        assert!(close_c(w.cos().acos(), w, 64.0));
        assert!(close_c(w.tan().atan(), w, 64.0));
        assert!(close_c(w.sinh().asinh(), w, 64.0));
        assert!(close_c(w.cosh().acosh(), w, 64.0));
        assert!(close_c(w.tanh().atanh(), w, 64.0));

        // Große Argumente laufen nicht in NaN.
        assert!(close_c(c::<T>(50.0, 1.0).tanh(), one, 8.0));
        assert!(close_c(c::<T>(1.0, 50.0).tan(), i, 8.0));
    }

    #[test]
    fn elementary_f32() {
        check::<f32>();
    }

    #[test]
    fn elementary_f64() {
        check::<f64>();
    }
}
