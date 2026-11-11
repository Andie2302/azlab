//! Geschlossene Lösungsformeln für Polynome bis Grad 4.
//!
//! Alle Funktionen nehmen die Koeffizienten aufsteigend und liefern die
//! verschiedenen reellen Nullstellen aufsteigend sortiert. Ist der höchste
//! Koeffizient null, wird die Gleichung niedrigeren Grades gelöst; sind alle
//! Koeffizienten null, ist das Ergebnis leer.
//!
//! Mehrfache Nullstellen erscheinen einmal, sofern die Rechnung sie exakt
//! gleich liefert. Durch Rundung kann eine doppelte Nullstelle auch als zwei
//! sehr nahe Werte oder gar nicht erscheinen. Wo keine Nullstelle fehlen
//! darf, ist [`isolate_real_roots`](super::isolate_real_roots) das richtige
//! Werkzeug.

use super::eval::horner_with_derivative;
use super::real_roots::RealRoots;
use crate::number::Complex;
use crate::scalar::{Elementary, RealElementary};

/// Verbessert eine Näherung mit wenigen Newton-Schritten am ursprünglichen
/// Polynom.
fn polish<T: RealElementary>(coefficients: &[T], start: T) -> T {
    let mut x = start;
    for _ in 0..3 {
        let (value, slope) = horner_with_derivative(coefficients, x);
        if slope == T::ZERO {
            break;
        }
        let next = x - value / slope;
        if !next.is_finite() {
            break;
        }
        x = next;
    }
    x
}

/// Löst `c[0] + c[1]·x = 0`.
pub fn solve_linear<T: RealElementary>(c: [T; 2]) -> RealRoots<T, 1> {
    let [c0, c1] = c;
    let mut roots = RealRoots::empty();
    if c1 != T::ZERO {
        roots.insert(-c0 / c1);
    }
    roots
}

/// Löst `c[0] + c[1]·x + c[2]·x² = 0` in der auslöschungsfreien Form.
pub fn solve_quadratic<T: RealElementary>(c: [T; 3]) -> RealRoots<T, 2> {
    let [c0, c1, c2] = c;
    if c2 == T::ZERO {
        return RealRoots::from_slice(solve_linear([c0, c1]).as_slice());
    }
    let mut roots = RealRoots::empty();
    let discriminant = c1 * c1 - T::from_f64(4.0) * c2 * c0;
    if discriminant < T::ZERO || discriminant.is_nan() {
        return roots;
    }
    // q hat den größeren Betrag der beiden möglichen Zähler; so entsteht
    // keine Differenz fast gleicher Zahlen.
    let q = -(c1 + discriminant.sqrt().copysign(c1)) * T::from_f64(0.5);
    roots.insert(q / c2);
    if q != T::ZERO {
        roots.insert(c0 / q);
    }
    roots
}

/// Löst `c[0] + c[1]·x + c[2]·x² + c[3]·x³ = 0`.
///
/// Cardano bei einer reellen Nullstelle, trigonometrische Form bei dreien;
/// danach Newton-Politur.
pub fn solve_cubic<T: RealElementary>(c: [T; 4]) -> RealRoots<T, 3> {
    let [c0, c1, c2, c3] = c;
    if c3 == T::ZERO {
        return RealRoots::from_slice(solve_quadratic([c0, c1, c2]).as_slice());
    }
    let two = T::from_f64(2.0);
    let third = T::ONE / T::from_f64(3.0);
    let (b2, b1, b0) = (c2 / c3, c1 / c3, c0 / c3);

    // Substitution x = t − b2/3 ergibt t³ + p·t + q = 0.
    let shift = b2 * third;
    let p = b1 - b2 * shift;
    let q = two * shift * shift * shift - shift * b1 + b0;
    let half_q = q * T::from_f64(0.5);
    let p_third = p * third;
    let discriminant = half_q * half_q + p_third * p_third * p_third;

    let mut roots = RealRoots::empty();
    if discriminant > T::ZERO {
        // Eine reelle Nullstelle: t = u + v mit u·v = −p/3.
        let u = (-half_q - discriminant.sqrt().copysign(half_q)).cbrt();
        let v = if u == T::ZERO { T::ZERO } else { -p_third / u };
        roots.insert(polish(&c, u + v - shift));
    } else if p_third == T::ZERO {
        // Dreifache Nullstelle.
        roots.insert(polish(&c, -shift));
    } else if discriminant <= T::ZERO {
        // Drei reelle Nullstellen: t = 2r·cos(φ + 2πk/3).
        let r = (-p_third).sqrt();
        let cosine = (-half_q / (r * r * r)).max(-T::ONE).min(T::ONE);
        let step = two * T::PI * third;
        let mut angle = cosine.acos() * third;
        for _ in 0..3 {
            roots.insert(polish(&c, two * r * angle.cos() - shift));
            angle += step;
        }
    }
    roots
}

/// Löst `c[0] + c[1]·x + … + c[4]·x⁴ = 0`.
///
/// Verfahren von Ferrari über die kubische Resolvente, danach
/// Newton-Politur. Numerisch die heikelste der Formeln.
pub fn solve_quartic<T: RealElementary>(c: [T; 5]) -> RealRoots<T, 4> {
    let [c0, c1, c2, c3, c4] = c;
    if c4 == T::ZERO {
        return RealRoots::from_slice(solve_cubic([c0, c1, c2, c3]).as_slice());
    }
    let two = T::from_f64(2.0);
    let eight = T::from_f64(8.0);
    let (b3, b2, b1, b0) = (c3 / c4, c2 / c4, c1 / c4, c0 / c4);

    // Substitution x = y − b3/4 ergibt y⁴ + p·y² + q·y + r = 0.
    let shift = b3 * T::from_f64(0.25);
    let shift2 = shift * shift;
    let p = b2 - T::from_f64(6.0) * shift2;
    let q = b1 - two * b2 * shift + eight * shift2 * shift;
    let r = b0 - b1 * shift + b2 * shift2 - T::from_f64(3.0) * shift2 * shift2;

    let mut roots = RealRoots::empty();
    if q == T::ZERO {
        // Biquadratisch: quadratische Gleichung in z = y².
        for &z in solve_quadratic([r, p, T::ONE]).as_slice() {
            if z >= T::ZERO {
                let y = z.sqrt();
                roots.insert(polish(&c, y - shift));
                roots.insert(polish(&c, -y - shift));
            }
        }
        return roots;
    }

    // Resolvente 8m³ + 8p·m² + (2p² − 8r)·m − q² = 0 hat für q ≠ 0 eine
    // positive Nullstelle. Damit zerfällt die Gleichung in
    // (y² + s·y + p/2 + m − q/(2s))·(y² − s·y + p/2 + m + q/(2s)) = 0,
    // s = √(2m).
    let resolvent = solve_cubic([-(q * q), two * p * p - eight * r, eight * p, eight]);
    let Some(&m) = resolvent.as_slice().last() else {
        return roots;
    };
    if m > T::ZERO {
        let s = (m + m).sqrt();
        let base = p * T::from_f64(0.5) + m;
        let offset = q / (s + s);
        for (linear, constant) in [(s, base - offset), (-s, base + offset)] {
            for &y in solve_quadratic([constant, linear, T::ONE]).as_slice() {
                roots.insert(polish(&c, y - shift));
            }
        }
    }
    roots
}

/// Löst `c[0] + c[1]·x + c[2]·x² = 0` mit komplexen Koeffizienten.
///
/// Liefert immer beide Nullstellen (eine doppelte zweimal). Ist `c[2]` null,
/// steht die Nullstelle der linearen Gleichung an beiden Plätzen.
pub fn solve_quadratic_complex<T: RealElementary>(c: [Complex<T>; 3]) -> [Complex<T>; 2] {
    let [c0, c1, c2] = c;
    let zero = Complex::new(T::ZERO, T::ZERO);
    if c2 == zero {
        let root = -c0 / c1;
        return [root, root];
    }
    let root_of_discriminant = (c1 * c1 - c2 * c0 * T::from_f64(4.0)).sqrt();
    let plus = c1 + root_of_discriminant;
    let minus = c1 - root_of_discriminant;
    let larger = if plus.norm() >= minus.norm() { plus } else { minus };
    let q = -larger * T::from_f64(0.5);
    let first = q / c2;
    let second = if q == zero { first } else { c0 / q };
    [first, second]
}
