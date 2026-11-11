//! Zusammengesetzte Trapez- und Simpsonregel.

use crate::scalar::Real;

/// Zusammengesetzte Trapezregel mit `intervals` gleich breiten
/// Teilintervallen (mindestens eines). Fehler ∝ Breite².
pub fn trapezoid_rule<T, F>(mut f: F, a: T, b: T, intervals: u32) -> T
where
    T: Real,
    F: FnMut(T) -> T,
{
    let intervals = intervals.max(1);
    let h = (b - a) / T::from_f64(f64::from(intervals));
    let mut sum = (f(a) + f(b)) * T::from_f64(0.5);
    for index in 1..intervals {
        sum += f(a + h * T::from_f64(f64::from(index)));
    }
    sum * h
}

/// Zusammengesetzte Simpsonregel mit `intervals` Teilintervallen (wird auf
/// eine gerade Zahl, mindestens zwei, aufgerundet). Fehler ∝ Breite⁴.
pub fn simpson<T, F>(mut f: F, a: T, b: T, intervals: u32) -> T
where
    T: Real,
    F: FnMut(T) -> T,
{
    let intervals = intervals.max(2).saturating_add(intervals % 2);
    let h = (b - a) / T::from_f64(f64::from(intervals));
    let (two, four) = (T::from_f64(2.0), T::from_f64(4.0));
    let mut sum = f(a) + f(b);
    for index in 1..intervals {
        let weight = if index % 2 == 1 { four } else { two };
        sum += weight * f(a + h * T::from_f64(f64::from(index)));
    }
    sum * h / T::from_f64(3.0)
}
