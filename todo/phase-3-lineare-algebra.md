# Phase 3 – Lineare Algebra fester Größe (Detailplan)

Status: **umgesetzt und geprüft** (volle Prüfmatrix grün).

## Abweichungen vom Plan bei der Umsetzung

- `Real` setzt jetzt auch `Conjugate<Real = Self>` voraus. So genügt in
  generischem Code die Schranke `T: Real`, um Normen und Zerlegungen zu
  nutzen.
- Schritte A und B sowie D bis F sind jeweils in einem Commit
  zusammengefasst.
- Die internen Makros liegen in `src/macros.rs` statt in `src/number/`,
  weil auch `linalg` sie braucht.
- `Matrix` hat zusätzlich `rows`, `cols`, `as_rows` und `zip_map`;
  `Vector` zusätzlich `IntoIterator` für Referenzen.
- `Cholesky` und `Lu` bieten `inverse()`; `Lu::solve_matrix` löst mehrere
  rechte Seiten.
- Hilfsdatei `src/linalg/decomp/util.rs`.

## Ziel

Vektoren und Matrizen fester Größe auf dem Stack, dazu die direkten
Zerlegungen LU, Cholesky und QR mit Lösen linearer Gleichungssysteme,
Determinante und Inverse. Alles generisch über den Elementtyp.

## Nicht-Ziele dieser Phase

- Keine dynamischen Größen, keine Sichten auf fremde Puffer (Phase 6).
- Keine Eigenwerte, keine SVD (Phase 11).
- Keine dünnbesetzten Matrizen, keine iterativen Löser (Phase 6).
- Keine Rotationsmatrizen/Euler-Winkel für Quaternionen (Phase 5).
- Keine Größenarithmetik im Typ (Blockmatrizen, Anhängen von Zeilen): auf
  stable nicht ausdrückbar.
- Keine Optimierung auf Geschwindigkeit (SIMD, Blockung).

## Abhängigkeiten und Features

- Keine neue Abhängigkeit, kein neues Feature, kein `cfg` im
  Bibliothekscode. Alles läuft in `no_std` ohne `alloc`.

## Entscheidungen

| Thema | Entscheidung |
|-------|--------------|
| Vektor | eigener Typ `Vector<T, N>`, kein Alias auf Matrix |
| Speicherordnung | zeilenweise: `Matrix<T, R, C>` über `[[T; C]; R]` |
| Indexzugriff | `get`/`get_mut` → `Option` **und** `v[i]`, `m[(i, j)]` wie Rust-Arrays (Panic bei falschem Index, dokumentiert) |
| Elementtypen der Zerlegungen | `Field + Conjugate`: reell, komplex (hermitesch), dual |
| Singularität | exakt null als Pivot → Fehler; keine Toleranzschwelle |

## Schritt A – Nachträge an den Skalar-Traits

- `Conjugate` bekommt `fn real_part(self) -> Self::Real` (Realteil; für
  Cholesky nötig).
- `Dual<T>`/`DualN<T, N>` implementieren `Conjugate` für `T: Real` mit
  `Real = Self` (Konjugation ohne Wirkung). Damit laufen Ableitungen durch
  alle Zerlegungen.

## Schritt B – `Vector<T, N>`

```rust
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector<T, const N: usize>(pub [T; N]);
```

| Für | Angebot |
|-----|---------|
| alle `T` | `new`, `From<[T; N]>` und zurück, `as_array`, `as_slice`, `iter`, `get`, `get_mut`, `Index`/`IndexMut<usize>`, `map`, `zip_map` |
| `T: Algebra` | `ZERO`, `splat`, `unit(i) -> Option`, `+ - ` Negation, `* T`, Zuweisungsformen, `dot` (bilinear), `component_mul`, `sum`, `outer` → Matrix |
| `T: Field` | `/ T` |
| `f32`, `f64` | `T * Vector` |
| `T: Conjugate` | `inner` (hermitesch: erster Faktor konjugiert), `norm_sqr` |
| `T::Real: RealElementary` | `norm`, `normalize() -> Option`, `distance` |
| `N = 3` | `cross` |

## Schritt C – `Matrix<T, R, C>`

```rust
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Matrix<T, const R: usize, const C: usize>(pub [[T; C]; R]);
```

| Für | Angebot |
|-----|---------|
| alle `T` | `new`, `from_fn`, `From<[[T; C]; R]>`, `get`, `get_mut`, `Index`/`IndexMut<(usize, usize)>`, `row`, `col` → `Option<Vector>`, `transpose`, `map` |
| `T: Algebra` | `ZERO`, `+ -` Negation, `* T`, `Matrix * Matrix` (Größen im Typ geprüft), `Matrix * Vector`, Zuweisungsformen |
| `T: Field` | `/ T` |
| `T: Conjugate` | `adjoint` (konjugiert transponiert), `norm_sqr` (Frobenius), `norm` |
| quadratisch, `T: Algebra` | `identity`, `from_diagonal`, `diagonal`, `trace` |
| quadratisch, `T: Ring` | `Algebra` und `Ring` für die Matrix selbst (damit z. B. `powi`-artige Algorithmen und Matrizen als Skalare funktionieren) |

## Schritt D – LU-Zerlegung

Kern-Algorithmen sind gegen `IndexMut<(usize, usize)>` geschrieben, nicht
gegen `Matrix`. Sichten und dynamische Matrizen (Phase 6) nutzen später
denselben Code.

```rust
pub struct Lu<T, const N: usize> { /* L und U kompakt, Permutation, Vorzeichen */ }

impl<T: Field + Conjugate, const N: usize> Matrix<T, N, N> {
    pub fn lu(self) -> Result<Lu<T, N>, LinalgError>;
    pub fn det(self) -> T;                                  // singulär → 0
    pub fn inverse(self) -> Result<Self, LinalgError>;
    pub fn solve(self, b: Vector<T, N>) -> Result<Vector<T, N>, LinalgError>;
}

impl<T: Field + Conjugate, const N: usize> Lu<T, N> {
    pub fn solve(&self, b: Vector<T, N>) -> Vector<T, N>;
    pub fn solve_matrix<const K: usize>(&self, b: Matrix<T, N, K>) -> Matrix<T, N, K>;
    pub fn det(&self) -> T;
    pub fn inverse(&self) -> Matrix<T, N, N>;
    pub fn l(&self) -> Matrix<T, N, N>;
    pub fn u(&self) -> Matrix<T, N, N>;
    pub fn permutation(&self) -> [usize; N];
}
```

- Spaltenpivotisierung nach größtem Betragsquadrat (`norm_sqr`).
- `det` nutzt für N ≤ 3 die geschlossene Formel, sonst LU.

## Schritt E – Cholesky-Zerlegung

Für hermitesch positiv definite Matrizen, `A = L·Lᴴ`.

```rust
pub struct Cholesky<T, const N: usize> { /* L */ }

impl<T: Field + Conjugate, const N: usize> Matrix<T, N, N>
where T::Real: RealElementary
{
    pub fn cholesky(self) -> Result<Cholesky<T, N>, LinalgError>;
}
// Cholesky: l(), solve(b), det(), inverse()
```

Gelesen wird nur das untere Dreieck. Nicht positiv definit → Fehler.

## Schritt F – QR-Zerlegung

Householder-Spiegelungen, auch komplex, für `R × C`-Matrizen.

```rust
pub struct Qr<T, const R: usize, const C: usize> { /* Q (R×R), R (R×C) */ }

impl<T: Field + Conjugate, const R: usize, const C: usize> Matrix<T, R, C>
where T::Real: RealElementary
{
    pub fn qr(self) -> Qr<T, R, C>;
}
// Qr: q(), r(), solve(b) -> Result<Vector<T, C>, _>  (kleinste Quadrate, R ≥ C)
```

## Fehlerbehandlung

```rust
#[non_exhaustive]
pub enum LinalgError { Singular, NotPositiveDefinite, DimensionMismatch }
```

Mit `Display` und `core::error::Error`. Einzige Panics: die Index-Operatoren
bei falschem Index (siehe Entscheidungen).

## Dateien

```
src/linalg.rs                       Deklarationen/Re-Exporte
src/linalg/error.rs
src/linalg/vector.rs                Deklarationen
src/linalg/vector/repr.rs  ops.rs  products.rs  tests.rs
src/linalg/matrix.rs                Deklarationen
src/linalg/matrix/repr.rs  ops.rs  square.rs  tests.rs
src/linalg/decomp.rs                Deklarationen
src/linalg/decomp/lu.rs  cholesky.rs  qr.rs  tests.rs
src/scalar/conjugate.rs, number/*   (Schritt A, angepasst)
src/lib.rs                          `pub mod linalg;` + Re-Exporte
```

## Tests

Erwartungswerte sind von Hand exakt herleitbar; Löser werden zusätzlich über
Rückeinsetzen (Residuum) geprüft. Jeder Test läuft für `f32` und `f64`.

- **Vector**: Summe, `dot([1,2,3],[4,5,6]) = 32`,
  `cross([1,2,3],[4,5,6]) = [−3,6,−3]`, `norm([3,4]) = 5`, Normieren der
  Null → `None`, falscher Index → `None`, komplex: `inner` ≠ `dot`.
- **Matrix**: `[[1,2],[3,4]]·[[5,6],[7,8]] = [[19,22],[43,50]]`,
  rechteckige Produkte, Transponieren, Matrix·Vektor, Einheitsmatrix, Spur.
- **LU**: `det` für 2×2, 3×3 und eine 4×4-Matrix mit nötigem Zeilentausch
  (Vorzeichen), `P·A = L·U`, Lösen mit Residuum, `A·A⁻¹ ≈ I`, singuläre
  Matrix → Fehler, komplexe Matrix, **Ableitung durch `solve` mit
  Dualzahlen** gegen die Formel von Hand.
- **Cholesky**: `[[4,2],[2,3]]` → `L = [[2,0],[1,√2]]`, `L·Lᴴ ≈ A` (auch
  komplex hermitesch), nicht positiv definit → Fehler.
- **QR**: `Q·R ≈ A`, `Qᴴ·Q ≈ I`, `R` obere Dreiecksform, kleinste Quadrate
  mit exakt lösbarem und mit überbestimmtem System (gegen die
  Normalgleichungen von Hand).

## Reihenfolge

Schritte A bis F, je ein Commit, nach jedem Schritt ein Prüflauf der KI
(Clippy und Tests mit Standard-Features und ohne Features). Am Ende volle
Prüfmatrix durch den Menschen.
