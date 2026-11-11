//! Tensor-Operationen, die die Stufe ändern.
//!
//! Das Ergebnis wird in einen Puffer des Aufrufers geschrieben. Seine Form
//! ergibt sich aus den Formen der Eingaben und ist bei jeder Funktion
//! genannt.

use super::view::TensorView;
use crate::linalg::LinalgError;
use crate::scalar::Algebra;

/// Äußeres Produkt: `out[i…, j…] = a[i…] · b[j…]`.
///
/// Form des Ergebnisses: die Form von `a`, gefolgt von der Form von `b`.
/// Fehler, wenn `out` nicht `a.len() · b.len()` Elemente hat.
pub fn outer<T: Algebra, const RA: usize, const RB: usize>(
    a: TensorView<'_, T, RA>,
    b: TensorView<'_, T, RB>,
    out: &mut [T],
) -> Result<(), LinalgError> {
    let (left, right) = (a.data(), b.data());
    if left.len().checked_mul(right.len()) != Some(out.len()) {
        return Err(LinalgError::DimensionMismatch);
    }
    if right.is_empty() {
        return Ok(());
    }
    for (chunk, &x) in out.chunks_exact_mut(right.len()).zip(left) {
        for (slot, &y) in chunk.iter_mut().zip(right) {
            *slot = x * y;
        }
    }
    Ok(())
}

/// Verallgemeinertes Produkt: Summe über den letzten Index von `a` und den
/// ersten von `b`, `out[i…, k…] = Σⱼ a[i…, j] · b[j, k…]`.
///
/// Für Stufe 1 und 1 das Skalarprodukt, für 2 und 1 Matrix mal Vektor, für
/// 2 und 2 das Matrixprodukt. Form des Ergebnisses: die Form von `a` ohne
/// den letzten Index, gefolgt von der Form von `b` ohne den ersten.
///
/// Fehler, wenn eine Eingabe die Stufe null hat, die beiden summierten
/// Indizes verschieden lang oder null lang sind oder `out` die falsche
/// Länge hat.
pub fn dot<T: Algebra, const RA: usize, const RB: usize>(
    a: TensorView<'_, T, RA>,
    b: TensorView<'_, T, RB>,
    out: &mut [T],
) -> Result<(), LinalgError> {
    let (Some(&inner), Some(&first)) = (a.shape().last(), b.shape().first()) else {
        return Err(LinalgError::DimensionMismatch);
    };
    if inner != first || inner == 0 {
        return Err(LinalgError::DimensionMismatch);
    }
    let (left, right) = (a.data(), b.data());
    let rows = left.len() / inner;
    let cols = right.len() / inner;
    if rows.checked_mul(cols) != Some(out.len()) {
        return Err(LinalgError::DimensionMismatch);
    }
    if cols == 0 {
        return Ok(());
    }
    for (out_row, left_row) in out.chunks_exact_mut(cols).zip(left.chunks_exact(inner)) {
        out_row.fill(T::ZERO);
        for (&x, right_row) in left_row.iter().zip(right.chunks_exact(cols)) {
            for (slot, &y) in out_row.iter_mut().zip(right_row) {
                *slot += x * y;
            }
        }
    }
    Ok(())
}

/// Kontraktion: Summe über die Diagonale der Indizes `first` und `second`
/// desselben Tensors. Für Stufe 2 die Spur.
///
/// Form des Ergebnisses: die Form von `a` ohne die beiden Indizes, in
/// unveränderter Reihenfolge.
///
/// Fehler, wenn die Indizes gleich sind oder außerhalb liegen, verschieden
/// lang sind oder `out` die falsche Länge hat.
// Die Schleifen brauchen den Index für Form und Position zugleich.
#[allow(clippy::needless_range_loop)]
pub fn contract<T: Algebra, const RANK: usize>(
    a: TensorView<'_, T, RANK>,
    first: usize,
    second: usize,
    out: &mut [T],
) -> Result<(), LinalgError> {
    let shape = a.shape();
    let (Some(&dim_first), Some(&dim_second)) = (shape.get(first), shape.get(second)) else {
        return Err(LinalgError::DimensionMismatch);
    };
    if first == second || dim_first != dim_second {
        return Err(LinalgError::DimensionMismatch);
    }
    let mut expected = Some(1_usize);
    for axis in 0..RANK {
        if axis != first && axis != second {
            expected = expected.and_then(|count| count.checked_mul(shape[axis]));
        }
    }
    if expected != Some(out.len()) {
        return Err(LinalgError::DimensionMismatch);
    }

    out.fill(T::ZERO);
    for (flat, &value) in a.data().iter().enumerate() {
        // Position in die einzelnen Indizes zerlegen, letzter Index zuerst.
        let mut rest = flat;
        let mut index = [0_usize; RANK];
        for axis in (0..RANK).rev() {
            index[axis] = rest % shape[axis];
            rest /= shape[axis];
        }
        if index[first] != index[second] {
            continue;
        }
        let mut target = 0_usize;
        for axis in 0..RANK {
            if axis != first && axis != second {
                target = target * shape[axis] + index[axis];
            }
        }
        if let Some(slot) = out.get_mut(target) {
            *slot += value;
        }
    }
    Ok(())
}
