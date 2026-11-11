//! Bequemer Einstieg in den Rückwärtsmodus.

use super::tape::{Tape, TapeError};
use super::var::Var;
use crate::linalg::Vector;
use crate::scalar::Algebra;

/// Wert und Gradient von `f` an der Stelle `x` im Rückwärtsmodus.
///
/// Legt ein Band mit `CAP` Plätzen auf dem Stack an; jede Variable und jede
/// Rechenoperation in `f` belegt einen. Reicht das nicht, ist das Ergebnis
/// [`TapeError::CapacityExceeded`].
///
/// `f` muss für jede Lebensdauer des Bands gelten; eine generische Funktion
/// übergibt man deshalb in einem Closure: `reverse_gradient(|p| f(p), x)`.
pub fn reverse_gradient<T, F, const N: usize, const CAP: usize>(
    f: F,
    x: Vector<T, N>,
) -> Result<(T, Vector<T, N>), TapeError>
where
    T: Algebra + 'static,
    F: for<'t> FnOnce(Vector<Var<'t, T, Tape<T, CAP>>, N>) -> Var<'t, T, Tape<T, CAP>>,
{
    let tape = Tape::<T, CAP>::new();
    let variables = x.map(|value| Var::new(&tape, value));
    let output = f(variables);
    let gradient = output.gradient_vector(variables)?;
    Ok((output.value(), gradient))
}
