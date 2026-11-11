//! Elementmatrizen und -vektoren der linearen Standardelemente.
//!
//! „Steifigkeit“ ist die Matrix zu `∫ ∇φᵢ·∇φⱼ` (Laplace-Operator), „Masse“
//! die zu `∫ φᵢ·φⱼ`, „Last“ der Vektor zu `∫ φᵢ` (rechte Seite mit
//! Quelldichte eins; für andere Dichten skalieren).

use super::error::FemError;
use crate::linalg::{Matrix, Vector};
use crate::scalar::{Real, RealElementary};

/// Prüft ein Maß (Länge, Fläche, Volumen) und gibt seinen Betrag zurück.
fn measure<T: Real>(signed: T) -> Result<T, FemError> {
    let size = signed.abs();
    if size > T::ZERO && size.is_finite() {
        Ok(size)
    } else {
        Err(FemError::DegenerateElement)
    }
}

/// Massenmatrix eines Simplex mit `V` Ecken und Maß `size`:
/// `size / ((V)(V + 1)) · (1 + δᵢⱼ)`.
fn simplex_mass<T: Real, const V: usize>(size: T) -> Matrix<T, V, V> {
    let corners = T::from_f64(V as f64);
    let unit = size / (corners * (corners + T::ONE));
    Matrix::from_fn(|i, j| if i == j { unit + unit } else { unit })
}

// ------------------------------------------------------------ Strecke ----

/// Steifigkeitsmatrix des linearen Streckenelements.
pub fn line2_stiffness<T: Real>(nodes: &[[T; 1]; 2]) -> Result<Matrix<T, 2, 2>, FemError> {
    let inverse = measure(nodes[1][0] - nodes[0][0])?.recip();
    Ok(Matrix([[inverse, -inverse], [-inverse, inverse]]))
}

/// Massenmatrix des linearen Streckenelements.
pub fn line2_mass<T: Real>(nodes: &[[T; 1]; 2]) -> Result<Matrix<T, 2, 2>, FemError> {
    Ok(simplex_mass(measure(nodes[1][0] - nodes[0][0])?))
}

/// Lastvektor des linearen Streckenelements.
pub fn line2_load<T: Real>(nodes: &[[T; 1]; 2]) -> Result<Vector<T, 2>, FemError> {
    Ok(Vector::splat(measure(nodes[1][0] - nodes[0][0])? * T::from_f64(0.5)))
}

// ------------------------------------------------------------ Dreieck ----

/// Doppelte vorzeichenbehaftete Fläche und die Koeffizienten `b`, `c` der
/// Gradienten: `∇φᵢ = (bᵢ, cᵢ) / (2·Fläche)`.
fn triangle_geometry<T: Real>(nodes: &[[T; 2]; 3]) -> (T, [T; 3], [T; 3]) {
    let [[x0, y0], [x1, y1], [x2, y2]] = *nodes;
    let twice_area = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    (twice_area, [y1 - y2, y2 - y0, y0 - y1], [x2 - x1, x0 - x2, x1 - x0])
}

/// Steifigkeitsmatrix des linearen Dreiecks.
pub fn triangle3_stiffness<T: Real>(nodes: &[[T; 2]; 3]) -> Result<Matrix<T, 3, 3>, FemError> {
    let (twice_area, b, c) = triangle_geometry(nodes);
    let scale = (measure(twice_area)? * T::from_f64(2.0)).recip();
    Ok(Matrix::from_fn(|i, j| (b[i] * b[j] + c[i] * c[j]) * scale))
}

/// Massenmatrix des linearen Dreiecks.
pub fn triangle3_mass<T: Real>(nodes: &[[T; 2]; 3]) -> Result<Matrix<T, 3, 3>, FemError> {
    let (twice_area, _, _) = triangle_geometry(nodes);
    Ok(simplex_mass(measure(twice_area)? * T::from_f64(0.5)))
}

/// Lastvektor des linearen Dreiecks.
pub fn triangle3_load<T: Real>(nodes: &[[T; 2]; 3]) -> Result<Vector<T, 3>, FemError> {
    let (twice_area, _, _) = triangle_geometry(nodes);
    Ok(Vector::splat(measure(twice_area)? / T::from_f64(6.0)))
}

/// Stoffgesetz der ebenen Elastizität.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlaneModel {
    /// Ebener Spannungszustand: dünne Scheibe, keine Spannung senkrecht zur
    /// Ebene.
    Stress,
    /// Ebener Verzerrungszustand: langer Körper, keine Dehnung senkrecht zur
    /// Ebene.
    Strain,
}

/// Steifigkeitsmatrix des linearen Dreiecks für ebene lineare Elastizität
/// (konstante Dehnung im Element).
///
/// Freiheitsgrade: je Knoten die Verschiebungen `(u, v)`, also
/// `(u₀, v₀, u₁, v₁, u₂, v₂)`. `youngs_modulus` ist der Elastizitätsmodul,
/// `poisson_ratio` die Querkontraktionszahl, `thickness` die Dicke.
pub fn triangle3_elasticity<T: Real>(
    nodes: &[[T; 2]; 3],
    youngs_modulus: T,
    poisson_ratio: T,
    thickness: T,
    model: PlaneModel,
) -> Result<Matrix<T, 6, 6>, FemError> {
    let (twice_area, b, c) = triangle_geometry(nodes);
    let area = measure(twice_area)? * T::from_f64(0.5);
    let (one, half, nu) = (T::ONE, T::from_f64(0.5), poisson_ratio);
    // Stoffmatrix D: Spannung = D · Dehnung.
    let (factor, diagonal, coupling, shear) = match model {
        PlaneModel::Stress => (youngs_modulus / (one - nu * nu), one, nu, (one - nu) * half),
        PlaneModel::Strain => (
            youngs_modulus / ((one + nu) * (one - nu - nu)),
            one - nu,
            nu,
            (one - nu - nu) * half,
        ),
    };
    let material = Matrix([
        [diagonal, coupling, T::ZERO],
        [coupling, diagonal, T::ZERO],
        [T::ZERO, T::ZERO, shear],
    ]) * factor;
    // Dehnungs-Verschiebungs-Matrix B (mal 2·Fläche).
    let strain = Matrix::<T, 3, 6>::from_fn(|row, col| {
        let (node, component) = (col / 2, col % 2);
        match (row, component) {
            (0, 0) | (2, 1) => b[node],
            (1, 1) | (2, 0) => c[node],
            _ => T::ZERO,
        }
    });
    let scale = thickness / (twice_area * twice_area) * area;
    Ok(strain.transpose() * material * strain * scale)
}

// ---------------------------------------------------------- Tetraeder ----

/// Volumen und Gradienten der vier Formfunktionen.
fn tetra_geometry<T: Real>(nodes: &[[T; 3]; 4]) -> Result<(T, [Vector<T, 3>; 4]), FemError> {
    let origin = Vector(nodes[0]);
    // Spalten: Kanten von Knoten 0 zu den Knoten 1, 2, 3.
    let jacobian =
        Matrix::<T, 3, 3>::from_fn(|row, col| nodes[col + 1][row] - origin[row]);
    let volume = measure(jacobian.det())? / T::from_f64(6.0);
    let inverse = jacobian.inverse().map_err(|_| FemError::DegenerateElement)?;
    // Zeile i der Inversen ist der Gradient der Formfunktion i + 1.
    let rows = inverse.0.map(Vector);
    let first = -(rows[0] + rows[1] + rows[2]);
    Ok((volume, [first, rows[0], rows[1], rows[2]]))
}

/// Steifigkeitsmatrix des linearen Tetraeders.
pub fn tetra4_stiffness<T: Real>(nodes: &[[T; 3]; 4]) -> Result<Matrix<T, 4, 4>, FemError> {
    let (volume, gradients) = tetra_geometry(nodes)?;
    Ok(Matrix::from_fn(|i, j| gradients[i].dot(gradients[j]) * volume))
}

/// Massenmatrix des linearen Tetraeders.
pub fn tetra4_mass<T: Real>(nodes: &[[T; 3]; 4]) -> Result<Matrix<T, 4, 4>, FemError> {
    Ok(simplex_mass(tetra_geometry(nodes)?.0))
}

/// Lastvektor des linearen Tetraeders.
pub fn tetra4_load<T: Real>(nodes: &[[T; 3]; 4]) -> Result<Vector<T, 4>, FemError> {
    Ok(Vector::splat(tetra_geometry(nodes)?.0 * T::from_f64(0.25)))
}

// ------------------------------------------------------------ Viereck ----

/// Integriert über das bilineare Viereck mit 2×2 Gauß-Punkten.
/// `integrand(Formfunktionen, Gradienten)` liefert die Matrix am Punkt.
fn quad4_integrate<T: RealElementary>(
    nodes: &[[T; 2]; 4],
    integrand: impl Fn([T; 4], [[T; 2]; 4]) -> Matrix<T, 4, 4>,
) -> Result<Matrix<T, 4, 4>, FemError> {
    // Ecken des Referenzquadrats [−1, 1]², gegen den Uhrzeigersinn.
    let corners = [[-T::ONE, -T::ONE], [T::ONE, -T::ONE], [T::ONE, T::ONE], [-T::ONE, T::ONE]];
    let quarter = T::from_f64(0.25);
    let point = (T::ONE / T::from_f64(3.0)).sqrt();
    let mut sum = Matrix::<T, 4, 4>::ZERO;
    for xi in [-point, point] {
        for eta in [-point, point] {
            let shape = corners.map(|[cx, cy]| (T::ONE + cx * xi) * (T::ONE + cy * eta) * quarter);
            let local = corners.map(|[cx, cy]| {
                [cx * (T::ONE + cy * eta) * quarter, cy * (T::ONE + cx * xi) * quarter]
            });
            // Jacobi-Matrix der Abbildung vom Referenzquadrat.
            let mut jacobian = Matrix::<T, 2, 2>::ZERO;
            for (derivative, node) in local.iter().zip(nodes) {
                for row in 0..2 {
                    for col in 0..2 {
                        jacobian[(row, col)] += derivative[row] * node[col];
                    }
                }
            }
            let determinant = jacobian.det();
            if determinant <= T::ZERO || !determinant.is_finite() {
                return Err(FemError::DegenerateElement);
            }
            let inverse = jacobian.inverse().map_err(|_| FemError::DegenerateElement)?;
            let gradients = local.map(|derivative| (inverse * Vector(derivative)).0);
            sum += integrand(shape, gradients) * determinant;
        }
    }
    Ok(sum)
}

/// Steifigkeitsmatrix des bilinearen Vierecks (Knoten gegen den
/// Uhrzeigersinn).
pub fn quad4_stiffness<T: RealElementary>(nodes: &[[T; 2]; 4]) -> Result<Matrix<T, 4, 4>, FemError> {
    quad4_integrate(nodes, |_, gradients| {
        Matrix::from_fn(|i, j| Vector(gradients[i]).dot(Vector(gradients[j])))
    })
}

/// Massenmatrix des bilinearen Vierecks (Knoten gegen den Uhrzeigersinn).
pub fn quad4_mass<T: RealElementary>(nodes: &[[T; 2]; 4]) -> Result<Matrix<T, 4, 4>, FemError> {
    quad4_integrate(nodes, |shape, _| Matrix::from_fn(|i, j| shape[i] * shape[j]))
}
