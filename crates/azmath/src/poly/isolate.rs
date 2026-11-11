//! Garantierte Einschließung aller reellen Nullstellen in einem Intervall.

use super::eval::horner;
use crate::interval::Interval;
use crate::scalar::{Algebra, Real};
use crate::solve::{ScalarResult, SolveError, Tolerance, brent};

/// Tiefe des Stapels für noch nicht untersuchte Teilintervalle.
const STACK: usize = 128;

/// Art einer gemeldeten Klammer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BracketKind {
    /// Genau eine Nullstelle: Das Polynom ist dort streng monoton und
    /// wechselt das Vorzeichen.
    Unique,
    /// Mögliche Nullstelle(n), nicht weiter trennbar – etwa mehrfache oder
    /// sehr dicht liegende Nullstellen, oder eine Nullstelle genau an einer
    /// Teilungsgrenze.
    Unresolved,
}

/// Intervall, das eine Nullstelle enthält oder enthalten kann.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RootBracket<T> {
    /// Das einschließende Intervall.
    pub interval: Interval<T>,
    /// Was über die Nullstellen darin bekannt ist.
    pub kind: BracketKind,
}

/// Einschluss von Wert und Ableitung des Polynoms über `x` (Horner mit
/// Intervallen).
fn enclose<T: Real>(coefficients: &[T], x: Interval<T>) -> (Interval<T>, Interval<T>) {
    let mut value = Interval::ZERO;
    let mut slope = Interval::ZERO;
    for &c in coefficients.iter().rev() {
        slope = slope * x + value;
        value = value * x + Interval::point(c);
    }
    (value, slope)
}

/// Wie [`enclose`], aber mit engerem Werte-Einschluss: Schnitt aus
/// Horner-Form und Mittelwertform `p(m) + p'(X)·(X − m)`.
fn enclose_tight<T: Real>(coefficients: &[T], x: Interval<T>) -> (Interval<T>, Interval<T>) {
    let (value, slope) = enclose(coefficients, x);
    let mid = Interval::point(x.mid());
    let (center, _) = enclose(coefficients, mid);
    let mean_value = center + slope * (x - mid);
    (value.intersect(mean_value).unwrap_or(value), slope)
}

/// Schließt alle reellen Nullstellen von `c[0] + c[1]·x + …` im Intervall
/// `search` ein.
///
/// Die Klammern werden aufsteigend nach `out` geschrieben; Rückgabe ist ihre
/// Anzahl.
///
/// **Garantie:** Jede reelle Nullstelle in `search` liegt in mindestens
/// einer gemeldeten Klammer – auch unter Rundungsfehlern, weil mit
/// Intervallarithmetik gerechnet wird. Es kann zu viel gemeldet werden
/// ([`BracketKind::Unresolved`], auch benachbarte Klammern zur selben
/// Nullstelle), nie zu wenig.
///
/// `min_width` ist die Breite, unter der nicht weiter geteilt wird.
///
/// Fehler: [`SolveError::InvalidInput`], wenn `search` nicht endlich,
/// `min_width` nicht positiv oder das Polynom identisch null ist;
/// [`SolveError::CapacityExceeded`], wenn `out` zu klein ist.
pub fn isolate_real_roots<T: Real>(
    coefficients: &[T],
    search: Interval<T>,
    min_width: T,
    out: &mut [RootBracket<T>],
) -> Result<usize, SolveError<(), T>> {
    if min_width <= T::ZERO
        || min_width.is_nan()
        || !search.lo().is_finite()
        || !search.hi().is_finite()
        || coefficients.iter().all(|c| *c == T::ZERO)
    {
        return Err(SolveError::InvalidInput);
    }

    let mut stack = [search; STACK];
    let mut depth = 1;
    let mut count = 0;
    while depth > 0 {
        depth -= 1;
        let Some(&current) = stack.get(depth) else {
            break;
        };
        let (value, slope) = enclose_tight(coefficients, current);
        if !value.contains_zero() {
            // Sicher keine Nullstelle.
            continue;
        }

        // Streng monoton: Die Vorzeichen an den Enden entscheiden, sofern
        // sie selbst sicher sind. Some(true) = genau eine Nullstelle,
        // Some(false) = keine, None = unentschieden.
        let mut verdict = None;
        if !slope.contains_zero() {
            let (at_lo, _) = enclose(coefficients, Interval::point(current.lo()));
            let (at_hi, _) = enclose(coefficients, Interval::point(current.hi()));
            if !at_lo.contains_zero() && !at_hi.contains_zero() {
                verdict = Some((at_lo.hi() < T::ZERO) != (at_hi.hi() < T::ZERO));
            }
        }

        let kind = match verdict {
            Some(false) => continue,
            Some(true) => BracketKind::Unique,
            None => {
                let mid = current.mid();
                let halves = if current.width() > min_width && depth + 2 <= STACK {
                    Interval::new(current.lo(), mid).zip(Interval::new(mid, current.hi()))
                } else {
                    None
                };
                match halves {
                    // Nur teilen, wenn beide Hälften echt kleiner werden.
                    Some((left, right)) if current.lo() < mid && mid < current.hi() => {
                        // Rechte Hälfte zuerst ablegen: So wird die linke
                        // zuerst untersucht und die Ausgabe ist aufsteigend.
                        if let Some(slot) = stack.get_mut(depth) {
                            *slot = right;
                        }
                        if let Some(slot) = stack.get_mut(depth + 1) {
                            *slot = left;
                        }
                        depth += 2;
                        continue;
                    }
                    _ => BracketKind::Unresolved,
                }
            }
        };

        let Some(slot) = out.get_mut(count) else {
            return Err(SolveError::CapacityExceeded);
        };
        *slot = RootBracket { interval: current, kind };
        count += 1;
    }
    Ok(count)
}

/// Verfeinert eine Klammer zu einer Nullstelle (Verfahren von Brent).
///
/// Für [`BracketKind::Unique`] gelingt das immer. Eine
/// [`BracketKind::Unresolved`]-Klammer hat nicht unbedingt einen
/// Vorzeichenwechsel; dann ist das Ergebnis [`SolveError::NoBracket`].
pub fn refine_root<T: Real>(
    coefficients: &[T],
    bracket: &RootBracket<T>,
    tol: &Tolerance<T>,
) -> ScalarResult<T> {
    brent(
        |x| horner(coefficients, x),
        bracket.interval.lo(),
        bracket.interval.hi(),
        tol,
    )
}
