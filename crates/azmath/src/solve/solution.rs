//! Der Typ [`Solution`].

/// Ergebnis eines Lösers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Solution<X, R> {
    /// Die Lösung bzw. die letzte Näherung.
    pub x: X,
    /// Restfehler: `|f(x)|` bzw. `‖F(x)‖`.
    pub residual: R,
    /// Anzahl der ausgeführten Iterationen.
    pub iterations: u32,
}
