//! Tests für [`Interval`].

use super::Interval;
use crate::scalar::{Algebra, Field, Real};
use crate::testutil::close;

fn iv<T: Real>(lo: f64, hi: f64) -> Interval<T> {
    Interval::new(T::from_f64(lo), T::from_f64(hi)).unwrap_or(Interval::ENTIRE)
}

/// `outer` enthält `[lo, hi]` und ist höchstens um Rundungsbreite größer.
fn encloses<T: Real>(outer: Interval<T>, lo: f64, hi: f64) -> bool {
    iv::<T>(lo, hi).is_subset_of(outer)
        && close(outer.lo(), T::from_f64(lo), 8.0)
        && close(outer.hi(), T::from_f64(hi), 8.0)
}

fn arithmetic<T: Real>() {
    let a = iv::<T>(1.0, 2.0);
    let b = iv::<T>(-3.0, 4.0);

    assert!(encloses(a + b, -2.0, 6.0));
    assert!(encloses(a - b, -3.0, 5.0));
    assert!(encloses(a * b, -6.0, 8.0));
    assert!(encloses(iv::<T>(-2.0, -1.0) * b, -8.0, 6.0));
    assert!(encloses(iv::<T>(6.0, 8.0) / iv(2.0, 4.0), 1.5, 4.0));
    assert!(encloses(iv::<T>(2.0, 4.0).recip(), 0.25, 0.5));
    assert_eq!(-b, iv(-4.0, 3.0));
    assert_eq!(a / b, Interval::ENTIRE);
    assert_eq!(a / iv(-1.0, 1.0), Interval::ENTIRE);

    // Rundung nach außen: 0,1 + 0,2 enthält 0,3.
    let sum = Interval::<T>::from_f64(0.1) + Interval::from_f64(0.2);
    assert!(sum.contains(T::from_f64(0.3)));
    assert!(sum.width() > T::ZERO);
    // Exakt darstellbare Konstanten bleiben Punkte.
    let half = Interval::<T>::from_f64(0.5);
    assert_eq!(half.lo(), half.hi());

    // Quadrat ist enger als das Produkt mit sich selbst.
    let c = iv::<T>(-2.0, 1.0);
    let square = c.sqr();
    assert_eq!(square.lo(), T::ZERO);
    assert!(square.hi() >= T::from_f64(4.0) && close(square.hi(), T::from_f64(4.0), 8.0));
    assert!((c * c).lo() < T::ZERO);
    assert!(square.is_subset_of(c.powi(2)));

    // x − x ist nicht null, enthält aber null.
    let difference = a - iv(1.0, 2.0);
    assert!(difference.contains_zero() && difference.width() > T::ZERO);

    let mut acc = a;
    acc += b;
    acc -= b;
    acc *= Interval::ONE;
    acc /= Interval::ONE;
    assert!(a.is_subset_of(acc));

    // Unendliche Grenzen und NaN.
    let product = Interval::<T>::ENTIRE * Interval::ZERO;
    assert!(product.contains_zero() && product.lo().is_finite() && product.hi().is_finite());
    assert_eq!(Interval::<T>::ENTIRE + a, Interval::ENTIRE);
    assert_eq!(Interval::<T>::point(T::NAN), Interval::ENTIRE);
    assert_eq!(Interval::<T>::hull_of(T::NAN, T::ONE), Interval::ENTIRE);
}

fn queries<T: Real>() {
    let a = iv::<T>(1.0, 3.0);
    let b = iv::<T>(2.0, 4.0);
    let one = T::ONE;
    let two = T::from_f64(2.0);

    assert!(Interval::new(two, one).is_none());
    assert!(Interval::new(T::NAN, one).is_none());
    assert_eq!(Interval::point(two), iv(2.0, 2.0));
    assert_eq!(Interval::hull_of(T::from_f64(3.0), one), a);

    assert_eq!((a.lo(), a.hi(), a.mid(), a.width()), (one, T::from_f64(3.0), two, two));
    assert_eq!(Interval::<T>::ENTIRE.mid(), T::ZERO);
    assert!(a.contains(one) && a.contains(two) && !a.contains(T::from_f64(3.5)));
    assert!(iv::<T>(-1.0, 1.0).contains_zero() && !a.contains_zero());
    assert!(iv::<T>(1.5, 2.5).is_subset_of(a) && !b.is_subset_of(a));

    assert_eq!(a.intersect(b), Some(iv(2.0, 3.0)));
    assert_eq!(a.intersect(iv(5.0, 6.0)), None);
    assert_eq!(a.hull(iv(5.0, 6.0)), iv(1.0, 6.0));

    assert!(a.certainly_lt(iv(4.0, 5.0)) && !a.certainly_lt(b));
    assert!(iv::<T>(4.0, 5.0).certainly_gt(a) && !b.certainly_gt(a));
    assert!(a.possibly_lt(b) && a.possibly_gt(b));
    assert!(!iv::<T>(4.0, 5.0).possibly_lt(a));

    assert_eq!(iv::<T>(-3.0, -1.0).abs(), a);
    assert_eq!(iv::<T>(-3.0, 2.0).abs(), iv(0.0, 3.0));
    assert_eq!(a.abs(), a);
}

#[test]
fn interval_f32() {
    arithmetic::<f32>();
    queries::<f32>();
}

#[test]
fn interval_f64() {
    arithmetic::<f64>();
    queries::<f64>();
}
