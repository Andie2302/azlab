# Phase 6, Teil 2 – GMRES, Tensoren, `alloc`-Typen (Detailplan)

Status: **umgesetzt, noch ungeprüft**. Die KI hat nichts ausgeführt; geprüft
wird nur durch `scripts/pruefen.sh` beim Menschen.

## Ziel

- GMRES als dritter iterativer Löser.
- Tensoren beliebiger Stufe als Sichten auf Puffer des Aufrufers, mit
  Umformen, äußerem Produkt, Kontraktion und verallgemeinertem Produkt.
- Hinter dem Feature `alloc`: besitzende Typen `DMatrix` und `DTensor`.

## Nicht-Ziele dieses Teils

- Keine Vorkonditionierer, kein komplexes GMRES.
- Keine Achsenvertauschung (Transponieren) von Tensoren, kein Broadcasting,
  keine Einstein-Summenschreibweise.
- `DMatrix`: nur Lösen über LU; keine eigenen Cholesky-/QR-Methoden.
- Keine dynamischen Vektortypen (`Vec<T>` genügt).

## Abhängigkeiten und Features

- Keine neue Abhängigkeit.
- Erstmals wird das Feature `alloc` benutzt (`extern crate alloc` in
  `lib.rs`, nur mit Feature). Ohne `alloc` ändert sich nichts.

## Entscheidungen

| Thema | Entscheidung |
|-------|--------------|
| Stufe eines Tensors | im Typ: `TensorView<T, RANK>` mit Form `[usize; RANK]` |
| Speicherordnung | zeilenweise verallgemeinert: der letzte Index läuft am schnellsten |
| Ergebnisse stufenändernder Operationen | in einen Puffer des Aufrufers; die Zielstufe gibt der Aufrufer an (stable Rust kann sie nicht ausrechnen) |
| GMRES | mit Neustart nach `restart` Schritten; Arbeitsspeicher vom Aufrufer, Größe über `gmres_work_len` |

## Umfang

Rund 1300 Zeilen Code und Tests in einem Commit.

## Schritt A – GMRES (`solve`)

```rust
pub fn gmres_work_len(size: usize, restart: usize) -> Option<usize>;
pub fn gmres<T, A>(a: &A, b: &[T], x: &mut [T], restart: usize, work: &mut [T], tol: &Tolerance<T>)
    -> LinearResult<T>;
```

Arnoldi-Verfahren mit modifiziertem Gram-Schmidt, Givens-Drehungen für das
kleine Ausgleichsproblem. Quelle für den Aufbau: Wikipedia, „Generalized
minimal residual method“, Abschnitt „Example code“.

## Schritt B – Tensor-Sichten (`tensor`)

```rust
pub struct TensorView<'a, T, const RANK: usize> { /* &[T], Form */ }
pub struct TensorViewMut<'a, T, const RANK: usize> { /* &mut [T], Form */ }
```

| Für | Angebot |
|-----|---------|
| beide | `new(data, shape) -> Result`, `shape`, `len`, `is_empty`, `get(index)`, `t[[i, j, k]]` |
| `TensorView` | `data`, `reshape(shape) -> Result` (auch auf andere Stufe) |
| `TensorViewMut` | `get_mut`, schreibender Index, `as_view`, `map_in_place`, `zip_in_place` |

## Schritt C – Tensor-Operationen (`tensor`)

```rust
pub fn outer(a, b, out: &mut [T]) -> Result<(), LinalgError>;            // äußeres Produkt
pub fn dot(a, b, out: &mut [T]) -> Result<(), LinalgError>;              // letzter Index von a mit erstem von b
pub fn contract(a, first, second, out: &mut [T]) -> Result<(), LinalgError>; // zwei Indizes desselben Tensors
```

`dot` verallgemeinert Skalarprodukt, Matrix·Vektor und Matrix·Matrix.

## Schritt D – Besitzende Typen (Feature `alloc`)

```rust
pub struct DMatrix<T> { /* Vec<T>, Zeilen, Spalten */ }
// zeros, identity, from_vec, from_fn, rows, cols, get, Index, view, view_mut,
// transpose, mul, mul_vec, solve; LinearOperator

pub struct DTensor<T, const RANK: usize> { /* Vec<T>, Form */ }
// zeros, from_vec, shape, view, view_mut, into_vec, outer, dot, contract
```

## Dateien

```
src/solve/iterative/gmres.rs
src/tensor.rs                    Deklarationen
src/tensor/view.rs  ops.rs  owned.rs (alloc)  tests.rs
src/linalg/dmatrix.rs            (alloc)
src/linalg/view.rs, linalg.rs, solve.rs, solve/iterative.rs, lib.rs   angepasst
```

## Tests

Struktur, Fehlerfälle und Vergleiche zweier Rechenwege; keine von der KI
ausgerechneten Erwartungswerte.

- GMRES: Lösung nahe der LU-Lösung, auch mit Neustart kleiner als die
  Systemgröße; Rückeinsetzen; Fehlerfälle.
- Tensor-Sichten: Lage der Elemente im Puffer, Umformen, falsche Formen.
- `outer` gegen `Vector::outer`, `dot` gegen `Matrix * Matrix` und
  `Matrix * Vector`, `contract` gegen die Spur; Spur des äußeren Produkts
  zweier Vektoren gleich ihrem Skalarprodukt.
- `DMatrix`: `mul`, `mul_vec`, `solve` gegen die Matrix fester Größe.
- `DTensor`: dieselben Ergebnisse wie die Operationen auf Sichten.
