//! Matrixfreier Operator über ein Netz und Zusammensetzen der rechten Seite.

use super::error::FemError;
use crate::linalg::{LinalgError, LinearOperator, Matrix, Vector};
use crate::scalar::Real;

/// Die Koordinaten der Knoten eines Elements; `None` bei ungültiger Nummer.
fn element_nodes<T: Real, const D: usize, const V: usize>(
    nodes: &[[T; D]],
    element: &[usize; V],
) -> Option<[[T; D]; V]> {
    let mut corners = [[T::ZERO; D]; V];
    for (corner, &index) in corners.iter_mut().zip(element) {
        *corner = *nodes.get(index)?;
    }
    Some(corners)
}

/// Die aus allen Elementmatrizen zusammengesetzte Matrix als Operator, ohne
/// sie zu speichern.
///
/// - `D`: Raumdimension, `V`: Knoten je Element,
/// - `M`: Größe der Elementmatrix; `M / V` ist die Zahl der Freiheitsgrade
///   je Knoten (1 für skalare Probleme, 2 für ebene Elastizität).
///
/// Freiheitsgrad `k` des Knotens `n` hat die Nummer `n · (M / V) + k`.
///
/// Mit einer Maske festgehaltener Freiheitsgrade (Dirichlet-Randbedingung)
/// wirkt der Operator dort wie die Einheitsmatrix und koppelt sie nicht an
/// die übrigen. So bleibt er symmetrisch und für konjugierte Gradienten
/// geeignet; die vorgegebenen Werte trägt [`lift`](Self::lift) in die rechte
/// Seite ein.
#[derive(Debug)]
pub struct ElementOperator<'a, T, F, const D: usize, const V: usize, const M: usize> {
    nodes: &'a [[T; D]],
    elements: &'a [[usize; V]],
    fixed: Option<&'a [bool]>,
    element_matrix: F,
    dofs_per_node: usize,
    size: usize,
}

impl<'a, T, F, const D: usize, const V: usize, const M: usize> ElementOperator<'a, T, F, D, V, M>
where
    T: Real,
    F: Fn(&[[T; D]; V]) -> Result<Matrix<T, M, M>, FemError>,
{
    /// Erzeugt den Operator und prüft Netz und Elemente.
    ///
    /// `fixed` markiert festgehaltene Freiheitsgrade (`true`) und muss so
    /// lang sein wie die Zahl der Freiheitsgrade. Fehler bei ungültigen
    /// Knotennummern, unpassender Matrixgröße, falscher Maskenlänge oder
    /// einem entarteten Element.
    pub fn new(
        nodes: &'a [[T; D]],
        elements: &'a [[usize; V]],
        fixed: Option<&'a [bool]>,
        element_matrix: F,
    ) -> Result<Self, FemError> {
        if V == 0 || M == 0 || M % V != 0 {
            return Err(FemError::InvalidMesh);
        }
        let dofs_per_node = M / V;
        let size = nodes.len().checked_mul(dofs_per_node).ok_or(FemError::InvalidMesh)?;
        if fixed.is_some_and(|mask| mask.len() != size) {
            return Err(FemError::LengthMismatch);
        }
        for element in elements {
            let corners = element_nodes(nodes, element).ok_or(FemError::InvalidMesh)?;
            element_matrix(&corners)?;
        }
        Ok(Self { nodes, elements, fixed, element_matrix, dofs_per_node, size })
    }

    /// `true`, wenn Freiheitsgrad `dof` festgehalten ist.
    #[inline]
    fn is_fixed(&self, dof: usize) -> bool {
        self.fixed.is_some_and(|mask| mask.get(dof).copied().unwrap_or(false))
    }

    /// Nummern der Freiheitsgrade eines Elements.
    fn element_dofs(&self, element: &[usize; V]) -> [usize; M] {
        core::array::from_fn(|local| {
            element[local / self.dofs_per_node] * self.dofs_per_node + local % self.dofs_per_node
        })
    }

    /// Schleife über alle Elemente: `y += K·x`, nur in freie Zeilen.
    /// `skip_fixed_columns` behandelt festgehaltene Einträge von `x` als null.
    fn accumulate(&self, x: &[T], y: &mut [T], skip_fixed_columns: bool) -> Result<(), LinalgError> {
        if x.len() != self.size || y.len() != self.size {
            return Err(LinalgError::DimensionMismatch);
        }
        for element in self.elements {
            // Beim Erzeugen geprüft; ein Fehler hier hieße, dass sich die
            // Elementfunktion anders verhält als damals.
            let matrix = element_nodes(self.nodes, element)
                .and_then(|corners| (self.element_matrix)(&corners).ok())
                .ok_or(LinalgError::DimensionMismatch)?;
            let dofs = self.element_dofs(element);
            let local = Vector(dofs.map(|dof| {
                if skip_fixed_columns && self.is_fixed(dof) { T::ZERO } else { x[dof] }
            }));
            let image = matrix * local;
            for (&dof, &value) in dofs.iter().zip(image.iter()) {
                if !self.is_fixed(dof) {
                    y[dof] += value;
                }
            }
        }
        Ok(())
    }

    /// Trägt vorgegebene Werte festgehaltener Freiheitsgrade in die rechte
    /// Seite ein.
    ///
    /// `values` enthält an den festgehaltenen Stellen die Vorgaben (andere
    /// Einträge werden ignoriert). Danach gilt für das System mit diesem
    /// Operator: Die Lösung nimmt an den festgehaltenen Stellen die Vorgaben
    /// an, und in den freien Gleichungen ist ihr Einfluss berücksichtigt.
    /// `work` braucht so viele Elemente wie es Freiheitsgrade gibt.
    pub fn lift(&self, values: &[T], rhs: &mut [T], work: &mut [T]) -> Result<(), FemError> {
        if values.len() != self.size || rhs.len() != self.size || work.len() != self.size {
            return Err(FemError::LengthMismatch);
        }
        // work = nur die Vorgaben, sonst null; rhs −= K·work in freien Zeilen.
        for (dof, (slot, &value)) in work.iter_mut().zip(values).enumerate() {
            *slot = if self.is_fixed(dof) { -value } else { T::ZERO };
        }
        self.accumulate(work, rhs, false).map_err(|_| FemError::LengthMismatch)?;
        for (dof, (slot, &value)) in rhs.iter_mut().zip(values).enumerate() {
            if self.is_fixed(dof) {
                *slot = value;
            }
        }
        Ok(())
    }
}

impl<T, F, const D: usize, const V: usize, const M: usize> LinearOperator<T>
    for ElementOperator<'_, T, F, D, V, M>
where
    T: Real,
    F: Fn(&[[T; D]; V]) -> Result<Matrix<T, M, M>, FemError>,
{
    #[inline]
    fn rows(&self) -> usize {
        self.size
    }

    #[inline]
    fn cols(&self) -> usize {
        self.size
    }

    fn apply(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        if y.len() != self.size {
            return Err(LinalgError::DimensionMismatch);
        }
        y.fill(T::ZERO);
        self.accumulate(x, y, true)?;
        for (dof, (out, &input)) in y.iter_mut().zip(x).enumerate() {
            if self.is_fixed(dof) {
                *out = input;
            }
        }
        Ok(())
    }
}

/// Setzt die rechte Seite aus den Elementvektoren zusammen: `out += Σ fₑ`.
///
/// `element_vector` liefert zu den Knotenkoordinaten eines Elements dessen
/// Vektor mit `M` Einträgen (`M / V` Freiheitsgrade je Knoten). `out` muss so
/// viele Elemente haben wie es Freiheitsgrade gibt und wird nicht vorher
/// gelöscht.
pub fn assemble_load<T, F, const D: usize, const V: usize, const M: usize>(
    nodes: &[[T; D]],
    elements: &[[usize; V]],
    element_vector: F,
    out: &mut [T],
) -> Result<(), FemError>
where
    T: Real,
    F: Fn(&[[T; D]; V]) -> Result<Vector<T, M>, FemError>,
{
    if V == 0 || M == 0 || M % V != 0 {
        return Err(FemError::InvalidMesh);
    }
    let dofs_per_node = M / V;
    if nodes.len().checked_mul(dofs_per_node) != Some(out.len()) {
        return Err(FemError::LengthMismatch);
    }
    for element in elements {
        let corners = element_nodes(nodes, element).ok_or(FemError::InvalidMesh)?;
        let vector = element_vector(&corners)?;
        for (local, &value) in vector.iter().enumerate() {
            let dof = element[local / dofs_per_node] * dofs_per_node + local % dofs_per_node;
            out[dof] += value;
        }
    }
    Ok(())
}
