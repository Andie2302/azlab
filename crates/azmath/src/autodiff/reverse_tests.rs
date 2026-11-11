//! Tests des Rückwärtsmodus: Vergleich mit dem Vorwärtsmodus und mit
//! Ableitungen von Hand.

use super::{Tape, TapeError, Var, gradient, reverse_gradient};
use crate::linalg::{Matrix, Vector};
use crate::scalar::{Algebra, Field, Real};
use crate::testutil::close;

/// Generisch geschrieben: läuft mit Dualzahlen und mit `Var`.
fn rational<S: Field>(p: Vector<S, 3>) -> S {
    let [x, y, z] = p.0;
    (x * y + z) / (x * x + S::from_f64(2.0)) - y * z.powi(3)
}

fn algebra<T: Real + 'static>() {
    let f = T::from_f64;
    let point = Vector([f(1.5), f(-0.5), f(2.0)]);

    // Gleicher Gradient wie im Vorwärtsmodus.
    let forward = gradient(rational, point);
    let reverse = reverse_gradient::<T, _, 3, 64>(|p| rational(p), point);
    assert!(reverse.is_ok_and(|(value, slope)| {
        close(value, forward.0, 16.0)
            && slope.iter().zip(forward.1.iter()).all(|(&a, &b)| close(a, b, 64.0))
    }));

    // Von Hand: x²·y bei (3, 2) → (12, 9).
    let tape = Tape::<T, 64>::new();
    let (x, y) = (Var::new(&tape, f(3.0)), Var::new(&tape, f(2.0)));
    let product = x * x * y;
    assert_eq!(product.value(), f(18.0));
    assert_eq!(product.gradient([x, y]), Ok([f(12.0), f(9.0)]));
    // Ableitung nach einem Zwischenergebnis und nach sich selbst.
    let square = x * x;
    assert_eq!((square * y).gradient([square]), Ok([f(2.0)]));
    assert_eq!(product.gradient([product]), Ok([T::ONE]));
    // Unbeteiligte Variable: Ableitung null.
    let other = Var::new(&tape, f(7.0));
    assert_eq!(product.gradient([other]), Ok([T::ZERO]));
    assert_eq!(tape.len(), 7);
    assert!(!tape.is_empty());

    // Konstanten belegen keinen Platz und haben keine Ableitung.
    let before = tape.len();
    let constant = Var::constant(f(4.0)) * Var::ONE + Var::from_f64(1.0);
    assert_eq!(tape.len(), before);
    assert_eq!(constant.value(), f(5.0));
    assert_eq!(constant.gradient([x]), Ok([T::ZERO]));
    assert_eq!((x * constant).gradient([x, constant]), Ok([f(5.0), T::ZERO]));
    assert_eq!(Var::<T, Tape<T, 4>>::constant(f(1.0)).gradient([]), Ok([]));

    // Vergleiche und Real-Methoden sehen nur den Wert.
    assert!(x > y && x == Var::constant(f(3.0)));
    assert_eq!((-x).abs().gradient([x]), Ok([T::ONE]));
    assert_eq!(x.min(y).gradient([x, y]), Ok([T::ZERO, T::ONE]));
    assert_eq!(x.max(y).gradient([x, y]), Ok([T::ONE, T::ZERO]));
    assert_eq!(x.copysign(-y).gradient([x]), Ok([-T::ONE]));
    assert_eq!((x / y).gradient([x, y]), Ok([f(0.5), f(-0.75)]));
    assert_eq!((x - y).gradient([x, y]), Ok([T::ONE, -T::ONE]));
    let mut acc = x;
    acc += y;
    acc *= x;
    acc -= y;
    acc /= y;
    // ((x + y)·x − y)/y → ∂/∂x = (2x + y)/y = 4
    assert!(acc.gradient([x]).is_ok_and(|[slope]| close(slope, f(4.0), 16.0)));

    // Band zu klein: Werte stimmen, die Ableitung meldet den Fehler.
    let small = Tape::<T, 3>::new();
    let (a, b) = (Var::new(&small, f(2.0)), Var::new(&small, f(5.0)));
    let overflow = a * b + a * a;
    assert_eq!(overflow.value(), f(14.0));
    assert_eq!(overflow.gradient([a, b]), Err(TapeError::CapacityExceeded));
    assert_eq!(
        reverse_gradient::<T, _, 3, 4>(|p| rational(p), point),
        Err(TapeError::CapacityExceeded)
    );
    // Nach dem Löschen ist das Band wieder benutzbar, alte Variablen nicht.
    small.clear();
    assert!(small.is_empty());
    let c = Var::new(&small, f(3.0));
    assert_eq!((c * c).gradient([c]), Ok([f(6.0)]));
    assert_eq!(c.gradient([overflow]), Err(TapeError::UnknownVariable));
    let stale = Var::new(&tape, f(1.0));
    tape.clear();
    assert_eq!(stale.gradient([stale]), Err(TapeError::UnknownVariable));

    // Zwei Bänder: Variablen dürfen nicht verwechselt werden, auch wenn
    // sie auf beiden Bändern dieselbe Nummer tragen.
    let (left, right) = (Tape::<T, 16>::new(), Tape::<T, 16>::new());
    let (u, w) = (Var::new(&left, f(2.0)), Var::new(&right, f(3.0)));
    assert_eq!(u.gradient([w]), Err(TapeError::UnknownVariable));
    assert_eq!(u.gradient([u]), Ok([T::ONE]));
    let mixed = u + w;
    assert_eq!(mixed.value(), f(5.0));
    assert_eq!(mixed.gradient([u]), Err(TapeError::UnknownVariable));
    assert_eq!(mixed.gradient([w]), Err(TapeError::UnknownVariable));
    assert_eq!((mixed * u).gradient([u]), Err(TapeError::UnknownVariable));
    assert_eq!((-mixed).gradient([u]), Err(TapeError::UnknownVariable));
    assert_eq!((u * u).gradient([u]), Ok([f(4.0)]));
    // Nach dem Löschen sind alte Variablen ungültig, auch wenn ihre Nummer
    // inzwischen neu vergeben ist.
    left.clear();
    let fresh = Var::new(&left, f(7.0));
    assert_eq!(u.gradient([u]), Err(TapeError::UnknownVariable));
    assert_eq!(fresh.gradient([u]), Err(TapeError::UnknownVariable));
    assert_eq!((u * fresh).gradient([fresh]), Err(TapeError::UnknownVariable));
    assert_eq!((u * fresh).value(), f(14.0));
    assert_eq!((fresh * fresh).gradient([fresh]), Ok([f(14.0)]));

    // Durch einen ganzen Algorithmus hindurch: Lösen von A(p)·x = b.
    // x₀ = 2/(2p − 1) → dx₀/dp = −4/(2p − 1)², bei p = 1: −4.
    let big = Tape::<T, 128>::new();
    let p = Var::new(&big, T::ONE);
    let k = |value: f64| Var::constant(T::from_f64(value));
    let system = Matrix([[p, k(1.0)], [k(1.0), k(2.0)]]);
    assert!(system.solve(Vector([k(1.0), k(0.0)])).is_ok_and(|solution| {
        close(solution[0].value(), f(2.0), 8.0)
            && solution[0].gradient([p]).is_ok_and(|[slope]| close(slope, f(-4.0), 8.0))
            && solution[1].gradient([p]).is_ok_and(|[slope]| close(slope, f(2.0), 8.0))
    }));
    assert_eq!(<Var<T, Tape<T, 4>> as Algebra>::ZERO.value(), T::ZERO);
}

#[test]
fn algebra_f32() {
    algebra::<f32>();
}

#[test]
fn algebra_f64() {
    algebra::<f64>();
}

#[cfg(any(feature = "std", feature = "libm"))]
mod elementary {
    use super::*;
    use crate::scalar::{Elementary, RealElementary};

    /// Benutzt jede Elementarfunktion; generisch für beide Modi.
    fn everything<S: RealElementary>(p: Vector<S, 2>) -> S {
        let [x, y] = p.0;
        x.sqrt() + x.exp() * y.sin() - x.ln() * y.cos() + (x * y).tan()
            + x.exp2() + x.exp_m1() + x.log2() + x.log10() + x.ln_1p()
            + x.powf(y) + x.powf(S::from_f64(3.0))
            + (x * S::from_f64(0.3)).asin() + (y * S::from_f64(0.3)).acos() + x.atan()
            + x.sinh() + y.cosh() + x.tanh()
            + x.asinh() + (x + S::ONE).acosh() + (y * S::from_f64(0.3)).atanh()
            + x.cbrt() + y.atan2(x) + x.hypot(y) + x.mul_add(y, x)
            + x.floor() + y.ceil() + x.round() + y.trunc()
    }

    fn check<T: RealElementary + 'static>() {
        let point = Vector([T::from_f64(0.7), T::from_f64(1.3)]);
        let forward = gradient(everything, point);
        let reverse = reverse_gradient::<T, _, 2, 128>(|p| everything(p), point);
        assert!(reverse.is_ok_and(|(value, slope)| {
            close(value, forward.0, 64.0)
                && slope.iter().zip(forward.1.iter()).all(|(&a, &b)| close(a, b, 256.0))
        }));
        // Konstanter Exponent an negativer Basis: (x²)' = 2x.
        let tape = Tape::<T, 8>::new();
        let x = Var::new(&tape, T::from_f64(-2.0));
        assert!(x
            .powf(Var::constant(T::from_f64(2.0)))
            .gradient([x])
            .is_ok_and(|[slope]| close(slope, T::from_f64(-4.0), 16.0)));
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

#[cfg(feature = "alloc")]
mod growing {
    use super::*;
    use crate::autodiff::GrowTape;

    fn check<T: Real>() {
        let f = T::from_f64;
        let tape = GrowTape::<T>::new();
        assert!(tape.is_empty());
        // Lange Rechnung: Summe von 1000 Quadraten, Ableitung 2·Σ = 2000·x.
        let x = Var::new(&tape, f(0.5));
        let mut sum = Var::ZERO;
        for _ in 0..1000 {
            sum += x * x;
        }
        assert!(sum.gradient([x]).is_ok_and(|[slope]| close(slope, f(1000.0), 64.0)));
        assert!(tape.len() > 2000);
        let point = Vector([f(1.5), f(-0.5), f(2.0)]);
        let variables = point.map(|value| Var::new(&tape, value));
        let forward = gradient(rational, point);
        assert!(rational(variables)
            .gradient_vector(variables)
            .is_ok_and(|slope| slope.iter().zip(forward.1.iter()).all(|(&a, &b)| close(a, b, 64.0))));
        tape.clear();
        assert_eq!(tape.len(), 0);
        assert_eq!(x.gradient([x]), Err(TapeError::UnknownVariable));
        assert!(GrowTape::<T>::default().is_empty());
    }

    #[test]
    fn growing_f32() {
        check::<f32>();
    }

    #[test]
    fn growing_f64() {
        check::<f64>();
    }
}
