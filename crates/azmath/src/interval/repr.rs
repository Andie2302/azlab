//! Der Typ [`Interval`].

use crate::scalar::Real;

/// Abgeschlossenes Intervall `[lo, hi]` mit `lo <= hi`.
///
/// Die Grenzen dürfen unendlich sein, aber nie NaN.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interval<T> {
    pub(super) lo: T,
    pub(super) hi: T,
}

impl<T: Real> Interval<T> {
    /// Die ganze Zahlengerade `[−∞, +∞]`.
    pub const ENTIRE: Self = Self { lo: T::NEG_INFINITY, hi: T::INFINITY };

    /// Erzeugt `[lo, hi]`. Gibt `None` zurück, wenn `lo > hi` oder eine
    /// Grenze NaN ist.
    #[inline]
    pub fn new(lo: T, hi: T) -> Option<Self> {
        if lo <= hi { Some(Self { lo, hi }) } else { None }
    }

    /// Das Intervall, das nur aus `x` besteht. Für NaN die ganze Gerade.
    #[inline]
    pub fn point(x: T) -> Self {
        if x.is_nan() { Self::ENTIRE } else { Self { lo: x, hi: x } }
    }

    /// Kleinstes Intervall, das `a` und `b` enthält, in beliebiger
    /// Reihenfolge. Für NaN die ganze Gerade.
    #[inline]
    pub fn hull_of(a: T, b: T) -> Self {
        if a.is_nan() || b.is_nan() {
            Self::ENTIRE
        } else {
            Self { lo: a.min(b), hi: a.max(b) }
        }
    }

    /// Intervall aus berechneten Grenzen, um eine Maschinenzahl nach außen
    /// erweitert. Fängt so den Rundungsfehler der Berechnung ab.
    #[inline]
    pub(super) fn outward(lo: T, hi: T) -> Self {
        if lo.is_nan() || hi.is_nan() {
            Self::ENTIRE
        } else {
            Self { lo: lo.next_down(), hi: hi.next_up() }
        }
    }

    /// Untere Grenze.
    #[inline]
    pub fn lo(self) -> T {
        self.lo
    }

    /// Obere Grenze.
    #[inline]
    pub fn hi(self) -> T {
        self.hi
    }

    /// Mittelpunkt. Für die ganze Gerade null.
    #[inline]
    pub fn mid(self) -> T {
        let half = T::from_f64(0.5);
        let mid = self.lo * half + self.hi * half;
        if mid.is_nan() { T::ZERO } else { mid }
    }

    /// Breite `hi − lo`.
    #[inline]
    pub fn width(self) -> T {
        self.hi - self.lo
    }

    /// `true`, wenn `x` im Intervall liegt.
    #[inline]
    pub fn contains(self, x: T) -> bool {
        self.lo <= x && x <= self.hi
    }

    /// `true`, wenn null im Intervall liegt.
    #[inline]
    pub fn contains_zero(self) -> bool {
        self.contains(T::ZERO)
    }

    /// `true`, wenn `self` ganz in `other` liegt.
    #[inline]
    pub fn is_subset_of(self, other: Self) -> bool {
        other.lo <= self.lo && self.hi <= other.hi
    }

    /// Durchschnitt; `None`, wenn sich die Intervalle nicht berühren.
    #[inline]
    pub fn intersect(self, other: Self) -> Option<Self> {
        Self::new(self.lo.max(other.lo), self.hi.min(other.hi))
    }

    /// Kleinstes Intervall, das beide enthält.
    #[inline]
    pub fn hull(self, other: Self) -> Self {
        Self { lo: self.lo.min(other.lo), hi: self.hi.max(other.hi) }
    }

    /// `true`, wenn jeder Punkt von `self` kleiner ist als jeder von `other`.
    #[inline]
    pub fn certainly_lt(self, other: Self) -> bool {
        self.hi < other.lo
    }

    /// `true`, wenn jeder Punkt von `self` größer ist als jeder von `other`.
    #[inline]
    pub fn certainly_gt(self, other: Self) -> bool {
        self.lo > other.hi
    }

    /// `true`, wenn irgendein Punkt von `self` kleiner ist als einer von
    /// `other`.
    #[inline]
    pub fn possibly_lt(self, other: Self) -> bool {
        self.lo < other.hi
    }

    /// `true`, wenn irgendein Punkt von `self` größer ist als einer von
    /// `other`.
    #[inline]
    pub fn possibly_gt(self, other: Self) -> bool {
        self.hi > other.lo
    }

    /// Beträge aller Punkte; exakt, ohne Erweiterung.
    #[inline]
    pub fn abs(self) -> Self {
        if self.lo >= T::ZERO {
            self
        } else if self.hi <= T::ZERO {
            Self { lo: -self.hi, hi: -self.lo }
        } else {
            Self { lo: T::ZERO, hi: (-self.lo).max(self.hi) }
        }
    }

    /// Quadrate aller Punkte. Enger als `self * self`, weil das Ergebnis nie
    /// negativ wird.
    #[inline]
    pub fn sqr(self) -> Self {
        let magnitude = self.abs();
        let widened = Self::outward(magnitude.lo * magnitude.lo, magnitude.hi * magnitude.hi);
        Self { lo: widened.lo.max(T::ZERO), hi: widened.hi }
    }
}
