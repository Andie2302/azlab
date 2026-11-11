//! Gemeinsamer Zugriff auf den Ableitungsanteil von `Dual` und `DualN`.

/// Ableitungsanteil: ein einzelner Wert oder ein Feld fester Länge.
pub(super) trait Tangent<T>: Copy {
    /// Wendet `f` auf jede Komponente an.
    fn map(self, f: impl Fn(T) -> T) -> Self;
    /// Verknüpft zwei Ableitungsanteile komponentenweise.
    fn zip(self, other: Self, f: impl Fn(T, T) -> T) -> Self;
}

impl<T: Copy> Tangent<T> for T {
    #[inline]
    fn map(self, f: impl Fn(T) -> T) -> Self {
        f(self)
    }

    #[inline]
    fn zip(self, other: Self, f: impl Fn(T, T) -> T) -> Self {
        f(self, other)
    }
}

impl<T: Copy, const N: usize> Tangent<T> for [T; N] {
    #[inline]
    fn map(self, f: impl Fn(T) -> T) -> Self {
        <[T; N]>::map(self, f)
    }

    #[inline]
    fn zip(self, other: Self, f: impl Fn(T, T) -> T) -> Self {
        let mut out = self;
        for (slot, rhs) in out.iter_mut().zip(other) {
            *slot = f(*slot, rhs);
        }
        out
    }
}
