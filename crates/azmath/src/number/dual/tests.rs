//! Tests für [`Dual`] und [`DualN`]. Alle Ableitungen sind gegen die
//! bekannte analytische Formel geprüft.

use super::{Dual, DualN};
use crate::number::Complex;
use crate::scalar::{Algebra, Field, Real};

fn algebra<T: Real>() {
    let f = |v: f64| T::from_f64(v);
    let x = Dual::variable(f(2.0));
    let three = Dual::constant(f(3.0));

    // x³ → 3x², Produktregel, 1/x → −1/x².
    let cube = x * x * x;
    assert_eq!((cube.re, cube.eps), (f(8.0), f(12.0)));
    let poly = x * x + three * x;
    assert_eq!((poly.re, poly.eps), (f(10.0), f(7.0)));
    let inv = x.recip();
    assert_eq!((inv.re, inv.eps), (f(0.5), f(-0.25)));
    let quot = (x * x) / x;
    assert_eq!((quot.re, quot.eps), (f(2.0), f(1.0)));
    let pow = x.powi(3);
    assert_eq!((pow.re, pow.eps), (f(8.0), f(12.0)));
    let diff = -(x - three);
    assert_eq!((diff.re, diff.eps), (f(1.0), f(-1.0)));

    let mut acc = x;
    acc += x;
    acc *= x;
    acc -= three;
    acc /= x;
    // (2x·x − 3)/x = 2x − 3/x → 2 + 3/x²
    assert_eq!((acc.re, acc.eps), (f(2.5), f(2.75)));

    // Vergleiche sehen nur den Wert.
    assert_eq!(Dual::new(f(1.0), f(5.0)), Dual::new(f(1.0), f(-3.0)));
    assert!(Dual::new(f(1.0), f(5.0)) < Dual::new(f(2.0), f(-3.0)));

    // Real-Methoden.
    let neg = Dual::new(f(-2.0), f(1.0));
    assert_eq!((neg.abs().re, neg.abs().eps), (f(2.0), f(-1.0)));
    assert_eq!(neg.min(x).eps, f(1.0));
    assert_eq!(neg.max(three).eps, f(0.0));
    let flipped = x.copysign(neg);
    assert_eq!((flipped.re, flipped.eps), (f(-2.0), f(-1.0)));
    assert_eq!(Dual::<T>::from_f64(1.5).eps, T::ZERO);
    assert!(Dual::<T>::NAN.is_nan());
    assert!(x.is_finite());

    // Zweite Ableitung über Schachtelung: (x³)'' = 6x.
    let xx = Dual::variable(Dual::variable(f(2.0)));
    let nested = xx * xx * xx;
    assert_eq!(nested.re.re, f(8.0));
    assert_eq!(nested.re.eps, f(12.0));
    assert_eq!(nested.eps.re, f(12.0));
    assert_eq!(nested.eps.eps, f(12.0));

    // Komplexe Ableitung: (z²)' = 2z.
    let z = Dual::variable(Complex::new(f(1.0), f(1.0)));
    assert_eq!((z * z).eps, Complex::new(f(2.0), f(2.0)));
}

fn gradient<T: Real>() {
    let f = |v: f64| T::from_f64(v);
    // f(x, y) = x²·y bei (3, 2): Gradient (2xy, x²) = (12, 9).
    let x = DualN::new(f(3.0), [f(1.0), f(0.0)]);
    let y = DualN::new(f(2.0), [f(0.0), f(1.0)]);
    let value = x * x * y;
    assert_eq!(value.re, f(18.0));
    assert_eq!(value.eps, [f(12.0), f(9.0)]);

    let quot = value / y;
    assert_eq!(quot.re, f(9.0));
    assert_eq!(quot.eps, [f(6.0), f(0.0)]);

    assert_eq!(
        DualN::<T, 2>::variable(f(3.0), 1).map(|d| d.eps),
        Some([f(0.0), f(1.0)])
    );
    assert!(DualN::<T, 2>::variable(f(3.0), 2).is_none());
    assert_eq!(DualN::<T, 3>::constant(f(1.0)).eps, [T::ZERO; 3]);
    assert_eq!(DualN::<T, 2>::ONE.eps, [T::ZERO; 2]);
}

#[test]
fn default_is_zero() {
    assert_eq!(DualN::<f64, 2>::default().eps, [0.0; 2]);
    assert_eq!(Dual::<f32>::default().eps, 0.0);
}

#[test]
fn algebra_f32() {
    algebra::<f32>();
    gradient::<f32>();
}

#[test]
fn algebra_f64() {
    algebra::<f64>();
    gradient::<f64>();
}

#[cfg(any(feature = "std", feature = "libm"))]
mod elementary {
    use super::*;
    use crate::scalar::{Elementary, RealElementary};
    use crate::testutil::close;

    fn check<T: RealElementary>() {
        let f = |v: f64| T::from_f64(v);
        let x0 = f(0.5);
        let x = Dual::variable(x0);
        let tol = 32.0;
        let derivative_is = |d: Dual<T>, expected: T| close(d.eps, expected, tol);

        assert!(derivative_is(x.sqrt(), (x0.sqrt() * f(2.0)).recip()));
        assert!(derivative_is(x.exp(), x0.exp()));
        assert!(derivative_is(x.exp2(), x0.exp2() * f(core::f64::consts::LN_2)));
        assert!(derivative_is(x.exp_m1(), x0.exp()));
        assert!(derivative_is(x.ln(), f(2.0)));
        assert!(derivative_is(x.log2(), f(2.0 / core::f64::consts::LN_2)));
        assert!(derivative_is(x.log10(), f(2.0 / core::f64::consts::LN_10)));
        assert!(derivative_is(x.ln_1p(), f(1.0 / 1.5)));
        assert!(derivative_is(x.sin(), x0.cos()));
        assert!(derivative_is(x.cos(), -x0.sin()));
        assert!(derivative_is(x.tan(), T::ONE + x0.tan() * x0.tan()));
        assert!(derivative_is(x.asin(), f(0.75).sqrt().recip()));
        assert!(derivative_is(x.acos(), -f(0.75).sqrt().recip()));
        assert!(derivative_is(x.atan(), f(0.8)));
        assert!(derivative_is(x.sinh(), x0.cosh()));
        assert!(derivative_is(x.cosh(), x0.sinh()));
        assert!(derivative_is(x.tanh(), T::ONE - x0.tanh() * x0.tanh()));
        assert!(derivative_is(x.asinh(), f(1.25).sqrt().recip()));
        assert!(derivative_is(x.atanh(), f(1.0 / 0.75)));
        let two = Dual::variable(f(2.0));
        assert!(derivative_is(two.acosh(), f(3.0).sqrt().recip()));

        // Potenzen: konstanter Exponent, negative Basis, x^x.
        assert!(derivative_is(x.powf(Dual::constant(f(3.0))), f(0.75)));
        let neg = Dual::variable(f(-2.0));
        assert!(derivative_is(neg.powf(Dual::constant(f(2.0))), f(-4.0)));
        assert!(derivative_is(x.powf(x), x0.powf(x0) * (x0.ln() + T::ONE)));

        // Kettenregel: sin(x²)' = 2x·cos(x²).
        assert!(derivative_is((x * x).sin(), f(0.25).cos()));

        // Reelle Zusatzfunktionen.
        let eight = Dual::variable(f(8.0));
        assert!(derivative_is(eight.cbrt(), f(1.0 / 12.0)));
        assert!(derivative_is(x.atan2(Dual::constant(f(2.0))), f(2.0 / 4.25)));
        assert!(derivative_is(Dual::constant(f(2.0)).atan2(x), f(-2.0 / 4.25)));
        let three = Dual::variable(f(3.0));
        assert!(derivative_is(three.hypot(Dual::constant(f(4.0))), f(0.6)));
        assert_eq!(x.floor().eps, T::ZERO);
        assert_eq!(x.ceil().eps, T::ZERO);
        assert_eq!(x.round().eps, T::ZERO);
        assert_eq!(x.trunc().eps, T::ZERO);
        // x·x + x → 2x + 1
        let fused = x.mul_add(x, x);
        assert!(close(fused.re, f(0.75), tol) && derivative_is(fused, f(2.0)));

        // Gradient von sqrt(x² + y²) bei (3, 4): (0.6, 0.8).
        let gx = DualN::new(f(3.0), [f(1.0), f(0.0)]);
        let gy = DualN::new(f(4.0), [f(0.0), f(1.0)]);
        let norm = (gx * gx + gy * gy).sqrt();
        assert!(close(norm.re, f(5.0), tol));
        assert!(close(norm.eps[0], f(0.6), tol) && close(norm.eps[1], f(0.8), tol));
        let hyp = gx.hypot(gy);
        assert!(close(hyp.eps[0], f(0.6), tol) && close(hyp.eps[1], f(0.8), tol));

        // Zweite Ableitung: sin'' = −sin.
        let xx = Dual::variable(Dual::variable(x0));
        assert!(close(xx.sin().eps.eps, -x0.sin(), tol));
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
