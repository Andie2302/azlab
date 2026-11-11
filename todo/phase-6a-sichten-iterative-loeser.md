# Phase 6, Teil 1 – Sichten, dünnbesetzte Matrizen, iterative Löser (Detailplan)

Status: **umgesetzt, noch ungeprüft**. Die KI hat nichts ausgeführt und
nichts nachgerechnet; geprüft wird nur durch `scripts/pruefen.sh` beim
Menschen.

## Warum zwei Teile

Phase 6 ist groß. Ohne Prüflauf dazwischen steigt mit der Blockgröße die
Zahl unentdeckter Fehler. Teil 1 reicht bereits für große lineare Systeme
(Navier-Stokes, Gewebe). Teil 2 folgt nach einem Prüflauf: GMRES,
Tensor-Sichten und -Operationen sowie die `alloc`-Typen `DMatrix`/`DTensor`.

## Ziel

Lineare Algebra, deren Größe erst zur Laufzeit feststeht – weiterhin ohne
Heap: Alle Daten liegen in Puffern des Aufrufers.

## Nicht-Ziele dieses Teils

- Kein GMRES, keine Tensoren, keine `alloc`-Typen (Teil 2).
- Keine Vorkonditionierer.
- Keine komplexen iterativen Löser (nur reelle Zahlen und Dualzahlen).
- Kein Aufbau dünnbesetzter Matrizen aus Tripeln; CSR-Daten kommen fertig
  vom Aufrufer.
- Keine Cholesky- und QR-Lösung auf Sichten (nur LU).

## Abhängigkeiten und Features

Keine neue Abhängigkeit, kein neues Feature, kein `cfg` im Bibliothekscode.

## Umfang

Rund 1000 Zeilen Code und Tests in einem Commit.

## Schritt A – Matrix-Sichten (`linalg`)

```rust
pub struct MatrixView<'a, T> { /* &[T], Zeilen, Spalten; zeilenweise */ }
pub struct MatrixViewMut<'a, T> { /* &mut [T], Zeilen, Spalten */ }
```

| Für | Angebot |
|-----|---------|
| beide | `new(data, rows, cols) -> Result` (Länge muss passen), `rows`, `cols`, `get`, `m[(i, j)]` |
| `MatrixViewMut` | `get_mut`, schreibender Index, `as_view`, `lu_in_place(perm)`, `cholesky_in_place()` |
| `MatrixView` | `lu_solve(perm, b, x)`: löst mit einer zuvor in-place zerlegten Matrix |

Die Zerlegungen nutzen dieselben Kernalgorithmen wie die Matrizen fester
Größe.

## Schritt B – `LinearOperator` (`linalg`)

```rust
pub trait LinearOperator<T> {
    fn rows(&self) -> usize;
    fn cols(&self) -> usize;
    /// y = A·x
    fn apply(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError>;
}
```

Implementiert für `Matrix<T, R, C>`, `MatrixView`, `MatrixViewMut` und
`CsrMatrix`. Matrixfreie Verfahren implementieren den Trait selbst.

## Schritt C – Dünnbesetzte Matrix im CSR-Format (`linalg`)

```rust
pub struct CsrMatrix<'a, T> { /* row_ptr, col_idx, values als Slices */ }
// new(rows, cols, row_ptr, col_idx, values) -> Result  (prüft die Struktur)
// get(i, j) -> Option<T>, nnz()
```

## Schritt D – Iterative Löser (`solve`)

```rust
pub fn conjugate_gradient<T, A>(a: &A, b: &[T], x: &mut [T], work: &mut [T], tol: &Tolerance<T>)
    -> Result<Solution<(), T>, SolveError<(), T>>;   // symmetrisch positiv definit; work: 3·n
pub fn bicgstab<T, A>(a: &A, b: &[T], x: &mut [T], work: &mut [T], tol: &Tolerance<T>)
    -> Result<Solution<(), T>, SolveError<(), T>>;   // allgemein; work: 6·n
```

- `x` enthält beim Aufruf den Startwert und danach die Lösung.
- Konvergiert, wenn `‖b − A·x‖ <= abs_f + rel_x · ‖b‖`.
- Zusammenbruch des Verfahrens → `Singular`, falsche Größen →
  `InvalidInput`, sonst `NotConverged` mit letztem Restfehler.

## Dateien

```
src/linalg/view.rs  operator.rs  sparse.rs  dynamic_tests.rs
src/linalg/decomp.rs             Kernalgorithmen crate-weit sichtbar
src/solve/iterative.rs           Deklarationen
src/solve/iterative/util.rs  cg.rs  bicgstab.rs  tests.rs
src/linalg.rs, src/solve.rs, src/lib.rs   angepasst
```

## Tests

Nur Strukturprüfungen, Fehlerfälle und Vergleiche zweier Rechenwege; die KI
rechnet keine Erwartungswerte aus.

- Sichten: Zugriff trifft das richtige Element des Puffers; falsche Länge →
  Fehler; LU auf einer Sicht liefert dieselbe Lösung wie LU auf der Matrix
  fester Größe.
- `LinearOperator`: alle Implementierungen liefern dasselbe wie
  `Matrix * Vector`; falsche Längen → Fehler.
- CSR: Strukturfehler werden erkannt; `get` stimmt mit der dichten Matrix
  überein.
- CG und BiCGSTAB: Restfehler klein, Lösung nahe der LU-Lösung; falsche
  Größen → `InvalidInput`.

## Risiko

Ungeprüfter Code. BiCGSTAB stammt aus dem Gedächtnis.
