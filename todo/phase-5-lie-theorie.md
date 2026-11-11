# Phase 5 – Lie-Theorie (Detailplan)

Status: **umgesetzt und geprüft** (volle Prüfmatrix grün).

## Ziel

Drehungen und starre Bewegungen in der Ebene und im Raum als Lie-Gruppen,
mit ihren Lie-Algebren, der Lie-Klammer, Exponential- und
Logarithmusabbildung und der adjungierten Darstellung. Dazu Kommutator und
Exponential für kleine Matrizen.

Meilenstein: ausreichend für Robotik-Kinematik.

## Nicht-Ziele dieser Phase

- Keine Jacobi-Matrizen der Gruppenoperationen (linke/rechte Jacobi-Matrix).
- Keine Euler-Winkel.
- Keine Interpolation auf SE(3), keine Mittelung von Drehungen.
- Keine weiteren Gruppen (Sim(3), SL(n), ...).
- Kein Matrix-Logarithmus.

## Abhängigkeiten und Features

Keine neue Abhängigkeit, kein neues Feature, kein `cfg` im Bibliothekscode.
Alles läuft in `no_std` ohne `alloc`.

## Entscheidungen

| Thema | Entscheidung |
|-------|--------------|
| Typnamen | `So2`, `So3`, `Se2`, `Se3` (Rust-Schreibweise statt `SO2`) |
| Reihenfolge im Tangentialvektor | erst Translation `rho`, dann Drehung |
| Darstellung | `So2`: Kosinus und Sinus; `So3`: Einheitsquaternion; `Se2`/`Se3`: Drehung plus Translation |
| Kleine Winkel | unter der Maschinengenauigkeit gelten die Grenzwerte der Reihen |

## Umfang

Rund 700 Zeilen Code und Tests in einem Commit.

## Schritt A – Traits

```rust
pub trait LieAlgebra: Copy + Add<Output = Self> + Sub<Output = Self> + Neg<Output = Self> {
    type Scalar: Real;
    const ZERO: Self;
    fn scale(self, factor: Self::Scalar) -> Self;
    fn bracket(self, other: Self) -> Self;          // Lie-Klammer [self, other]
}

pub trait LieGroup: Copy {
    type Algebra: LieAlgebra;
    fn identity() -> Self;
    fn compose(self, other: Self) -> Self;
    fn inverse(self) -> Self;
    fn exp(tangent: Self::Algebra) -> Self;
    fn log(self) -> Self::Algebra;
    fn adjoint(self, tangent: Self::Algebra) -> Self::Algebra;   // Ad_g(ξ)
}
```

## Schritt B – Lie-Algebren

| Typ | Felder | Klammer |
|-----|--------|---------|
| `So2Tangent<T>` | `angle` | null |
| `So3Tangent<T>` | `phi: Vector<T, 3>` | Kreuzprodukt; zusätzlich `hat()` → schiefsymmetrische Matrix |
| `Se2Tangent<T>` | `rho: Vector<T, 2>`, `theta` | gemischt aus Drehung und Translation |
| `Se3Tangent<T>` | `rho`, `phi: Vector<T, 3>` | gemischt aus Kreuzprodukten |

## Schritt C – Lie-Gruppen

| Typ | Zusätzlich zu `LieGroup` |
|-----|--------------------------|
| `So2<T>` | `from_angle`, `angle`, `rotate`, `matrix` |
| `So3<T>` | `from_quaternion -> Option`, `from_axis_angle -> Option`, `quaternion`, `rotate`, `matrix` |
| `Se2<T>`, `Se3<T>` | `new`, öffentliche Felder `rotation`, `translation`, `act` (Punkt transformieren) |

Alle brauchen `T: RealElementary`.

## Schritt D – Matrizen

```rust
impl Matrix<T, N, N> {
    pub fn commutator(self, other: Self) -> Self;   // A·B − B·A
    pub fn exp(self) -> Self;                       // Skalieren und Quadrieren mit Taylor-Reihe
}
```

## Dateien

```
src/lie.rs                       Deklarationen
src/lie/traits.rs  tangent.rs  so2.rs  so3.rs  se2.rs  se3.rs  tests.rs
src/linalg/matrix/expm.rs        commutator, exp
src/linalg/matrix.rs, src/lib.rs angepasst
```

## Tests

Nur über mathematische Identitäten mit frei gewählten Eingaben; die KI
rechnet keine Erwartungswerte aus.

- `exp(log(g)) ≈ g` und `log(exp(ξ)) ≈ ξ` für alle vier Gruppen.
- `g · g⁻¹ ≈ Identität`, Assoziativität der Verknüpfung.
- Adjungierte: `g · exp(ξ) · g⁻¹ ≈ exp(Ad_g ξ)`.
- Klammer: antisymmetrisch, Jacobi-Identität.
- Drehungen erhalten Längen; `matrix() · v ≈ rotate(v)`.
- `Se2`/`Se3`: `(a · b).act(p) ≈ a.act(b.act(p))`.
- Matrizen: `exp(0) = I`, `exp(A) · exp(−A) ≈ I`,
  `exp(hat(φ)) ≈ So3::exp(φ).matrix()`, Kommutator antisymmetrisch.

Jeder Test läuft für `f32` und `f64`.

## Risiko

Die Formeln für `exp`/`log` von SE(2)/SE(3) und die Klammern stammen aus dem
Gedächtnis und sind ungeprüft. Die Identitäts-Tests decken Fehler darin auf.
