//! Gammafunktion und Verwandte.
//!
//! `gamma` und `ln_gamma` nutzen die Lanczos-Näherung mit `g = 7` und neun
//! Koeffizienten (Quelle der Koeffizienten: Wikipedia, „Lanczos
//! approximation“; Genauigkeit etwa 15 Stellen). Die unvollständigen
//! Gammafunktionen nutzen Reihe und Kettenbruch nach Press u. a.,
//! „Numerical Recipes“, Abschnitt 6.2.

use crate::scalar::RealElementary;

/// Parameter `g` der Lanczos-Näherung.
const LANCZOS_G: f64 = 7.0;

/// Koeffizienten der Lanczos-Näherung für `g = 7`.
const LANCZOS: [f64; 9] = [
    0.999_999_999_999_809_9,
    676.520_368_121_885_1,
    -1_259.139_216_722_402_8,
    771.323_428_777_653_1,
    -176.615_029_162_140_6,
    12.507_343_278_686_905,
    -0.138_571_095_265_720_12,
    9.984_369_578_019_572e-6,
    1.505_632_735_149_311_6e-7,
];

/// Höchstzahl der Glieder von Reihe und Kettenbruch. Nahe `x ≈ a` wächst
/// der Bedarf etwa mit der Wurzel aus `a`; das reicht bis `a` um 10⁶.
const MAX_TERMS: usize = 10_000;

/// Die Lanczos-Summe `A(x)` und die Hilfsgröße `t = x + g − 1/2` für
/// `Γ(x) = √(2π)·t^(x − 1/2)·e^(−t)·A(x)`, gültig für `x >= 1/2`.
fn lanczos<T: RealElementary>(x: T) -> (T, T) {
    let shifted = x - T::ONE;
    let mut sum = T::from_f64(LANCZOS[0]);
    let mut denominator = shifted;
    for &coefficient in &LANCZOS[1..] {
        denominator += T::ONE;
        sum += T::from_f64(coefficient) / denominator;
    }
    (sum, shifted + T::from_f64(LANCZOS_G) + T::from_f64(0.5))
}

/// Gammafunktion `Γ(x)`; für natürliche Zahlen gilt `Γ(n) = (n − 1)!`.
///
/// An den Polstellen (null und negative ganze Zahlen) ist das Ergebnis
/// unendlich oder NaN.
pub fn gamma<T: RealElementary>(x: T) -> T {
    let half = T::from_f64(0.5);
    if x < half {
        // Spiegelung: Γ(x)·Γ(1 − x) = π / sin(πx)
        return T::PI / ((T::PI * x).sin() * gamma(T::ONE - x));
    }
    let (sum, t) = lanczos(x);
    // In zwei Faktoren potenzieren, damit t^(x − 1/2) nicht vorzeitig
    // überläuft.
    let power = t.powf((x - half) * half);
    (T::PI + T::PI).sqrt() * power * (-t).exp() * power * sum
}

/// Natürlicher Logarithmus von `|Γ(x)|`.
///
/// Bleibt auch für große `x` darstellbar, wo `Γ(x)` selbst überläuft.
pub fn ln_gamma<T: RealElementary>(x: T) -> T {
    let half = T::from_f64(0.5);
    if x < half {
        return (T::PI / (T::PI * x).sin().abs()).ln() - ln_gamma(T::ONE - x);
    }
    let (sum, t) = lanczos(x);
    (x - half) * t.ln() - t + ((T::PI + T::PI).sqrt() * sum).ln()
}

/// Betafunktion `B(a, b) = Γ(a)·Γ(b) / Γ(a + b)` für positive `a`, `b`.
pub fn beta<T: RealElementary>(a: T, b: T) -> T {
    (ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)).exp()
}

/// Gemeinsamer Vorfaktor `e^(−x)·x^a / Γ(a)`.
fn prefactor<T: RealElementary>(a: T, x: T) -> T {
    (a * x.ln() - x - ln_gamma(a)).exp()
}

/// `P(a, x)` über die Reihe; konvergiert schnell für `x < a + 1`. NaN, wenn
/// die Reihe nach [`MAX_TERMS`] Gliedern nicht konvergiert ist.
fn series<T: RealElementary>(a: T, x: T) -> T {
    let mut denominator = a;
    let mut term = a.recip();
    let mut sum = term;
    for _ in 0..MAX_TERMS {
        denominator += T::ONE;
        term *= x / denominator;
        sum += term;
        if term.abs() <= sum.abs() * T::EPSILON {
            return sum * prefactor(a, x);
        }
    }
    T::NAN
}

/// `Q(a, x)` über den Kettenbruch (Verfahren von Lentz); konvergiert schnell
/// für `x >= a + 1`. NaN, wenn er nach [`MAX_TERMS`] Gliedern nicht
/// konvergiert ist.
fn continued_fraction<T: RealElementary>(a: T, x: T) -> T {
    let tiny = T::EPSILON * T::EPSILON * T::EPSILON;
    let guard = |value: T| if value.abs() < tiny { tiny } else { value };
    let mut b = x + T::ONE - a;
    let mut c = tiny.recip();
    let mut d = guard(b).recip();
    let mut fraction = d;
    let mut index = T::ZERO;
    for _ in 0..MAX_TERMS {
        index += T::ONE;
        let numerator = -index * (index - a);
        b += T::from_f64(2.0);
        d = guard(numerator * d + b).recip();
        c = guard(b + numerator / c);
        let factor = d * c;
        fraction *= factor;
        if (factor - T::ONE).abs() <= T::EPSILON {
            return fraction * prefactor(a, x);
        }
    }
    T::NAN
}

/// Grenzwerte von `P(a, x)` für unendliche Argumente; `None`, wenn beide
/// endlich sind.
fn limit_p<T: RealElementary>(a: T, x: T) -> Option<T> {
    match (a.is_finite(), x.is_finite()) {
        (true, true) => None,
        // Das Integral bis unendlich ist die ganze Gammafunktion.
        (true, false) => Some(T::ONE),
        (false, true) => Some(T::ZERO),
        (false, false) => Some(T::NAN),
    }
}

/// `true`, wenn `a` so groß ist, dass `a + 1` nicht mehr von `a` zu
/// unterscheiden ist. Reihe und Kettenbruch zählen in Einerschritten und
/// liefern dann Unsinn.
fn unresolved<T: RealElementary>(a: T) -> bool {
    a + T::ONE == a
}

/// Regularisierte untere unvollständige Gammafunktion
/// `P(a, x) = γ(a, x) / Γ(a)` für `a > 0`, `x >= 0`.
///
/// Wächst von 0 bei `x = 0` auf 1 für `x → ∞`; für unendliches `x` ist das
/// Ergebnis 1. Außerhalb des Definitionsbereichs NaN.
///
/// Für sehr große `a` (über etwa 10⁶) mit `x` in der Nähe von `a`
/// konvergieren Reihe und Kettenbruch nicht in der vorgesehenen Zahl von
/// Gliedern; das Ergebnis ist dann NaN statt eines ungenauen Werts. Dasselbe
/// gilt, wenn `a` so groß ist, dass der Zahltyp `a + 1` nicht mehr von `a`
/// unterscheidet.
pub fn gamma_p<T: RealElementary>(a: T, x: T) -> T {
    if a <= T::ZERO || x < T::ZERO || a.is_nan() || x.is_nan() {
        return T::NAN;
    }
    if x == T::ZERO {
        return T::ZERO;
    }
    if let Some(limit) = limit_p(a, x) {
        return limit;
    }
    if unresolved(a) {
        return T::NAN;
    }
    if x < a + T::ONE { series(a, x) } else { T::ONE - continued_fraction(a, x) }
}

/// Regularisierte obere unvollständige Gammafunktion
/// `Q(a, x) = 1 − P(a, x)`, ohne Auslöschung für große `x`; für unendliches
/// `x` ist das Ergebnis 0. Grenzen wie bei [`gamma_p`].
pub fn gamma_q<T: RealElementary>(a: T, x: T) -> T {
    if a <= T::ZERO || x < T::ZERO || a.is_nan() || x.is_nan() {
        return T::NAN;
    }
    if x == T::ZERO {
        return T::ONE;
    }
    if let Some(limit) = limit_p(a, x) {
        return T::ONE - limit;
    }
    if unresolved(a) {
        return T::NAN;
    }
    if x < a + T::ONE { T::ONE - series(a, x) } else { continued_fraction(a, x) }
}
