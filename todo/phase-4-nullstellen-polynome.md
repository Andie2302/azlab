# Phase 4 – Nullstellen, Polynome, erste Löser (Detailplan)

Status: **umgesetzt und geprüft** (volle Prüfmatrix grün).

## Abweichungen vom Plan bei der Umsetzung

- Alles steckt in einem Commit statt in zehn: Ohne Prüfläufe ließ sich nicht
  sicherstellen, dass jeder Zwischenstand für sich übersetzt.
- `Real` bekommt zusätzlich die Konstante `NEG_INFINITY` (für
  `Interval::ENTIRE`).
- Mehrfache Nullstellen in den geschlossenen Formeln: Sie erscheinen nur
  dann einmal, wenn die Rechnung sie exakt gleich liefert. Durch Rundung
  kann eine doppelte Nullstelle bei Grad 3 und 4 auch als zwei sehr nahe
  Werte oder gar nicht erscheinen. Die Garantie liefert Schritt H.
- `Interval` erweitert jedes berechnete Ergebnis um eine Maschinenzahl nach
  außen, auch wenn die Rechnung exakt war. `powi` kommt aus `Field`
  (wiederholtes Quadrieren): gültiger, aber nicht engster Einschluss; für
  das Quadrat gibt es das engere `sqr`.
- `complex_roots` verlangt, dass der Ausgabepuffer genau so lang ist wie der
  Grad (zu kurz: `CapacityExceeded`, zu lang: `InvalidInput`).
- Autodiff-Hilfen haben zusätzlich `variables` (macht aus einem Punkt die
  unabhängigen Variablen).
- Newton für Systeme: Lässt sich der Restfehler nicht mehr verkleinern,
  obwohl der Schritt schon winzig ist, gilt das als konvergiert
  (Rechengenauigkeit erreicht), nicht als Fehler.
- Typ-Aliase `ScalarResult<T>` und `SystemResult<T, N>` für die Ergebnisse.
- Die Garantie-Tests für Grad 7 nutzen `x·(x²−1)·(x²−4)·(x²−9)` statt der
  Nullstellen 1 bis 7: Bei jenem Polynom wäre der Einschluss so grob, dass
  sehr viele unaufgelöste Klammern entstünden.

## Bekannte Grenze

Die Einschließung arbeitet mit Intervall-Horner und Mittelwertform. Bei
Polynomen mit sehr großen Koeffizienten und dicht liegenden Nullstellen
(Wilkinson-artig) bleibt sie korrekt, meldet aber viele unaufgelöste
Klammern. Für die Polynome der Kollisionserkennung (Grad 3) ist das kein
Problem.

## Ziel

Alles, was man zum Lösen von Gleichungen in einer und in mehreren
Unbekannten braucht, solange die Größe zur Übersetzungszeit feststeht:

- Polynome auswerten und ihre Nullstellen finden – auch **garantiert
  lückenlos** in einem Intervall, für jeden Grad.
- Nullstellen beliebiger Funktionen einer Veränderlichen.
- Nichtlineare Gleichungssysteme mit Newton, Jacobi-Matrix aus Autodiff.
- Intervallarithmetik als Grundlage für beweisbare Einschlüsse.

Meilenstein: ausreichend für ein IAPWS-95-Crate und – mit Phase 3 – für
kontinuierliche Kollisionserkennung in `azphysics`.

## Nicht-Ziele dieser Phase

- Keine Elementarfunktionen auf Intervallen (`sin`, `exp`, ...): die
  Backends garantieren dafür keine Rundungsrichtung.
- Keine Polynom-Multiplikation, -Division, kein ggT (Ergebnisgröße ist auf
  stable nicht im Typ ausdrückbar; bei Bedarf später über Puffer).
- Keine Polynome in mehreren Veränderlichen.
- Keine Ausgleichsrechnung für nichtlineare Probleme (Gauß-Newton,
  Levenberg-Marquardt) und keine Optimierung (Phase 10).
- Keine Systeme variabler Größe (Phase 6).
- Keine komplexen Newton-Systeme.

## Abhängigkeiten, Features, Toolchain

- Keine neue Abhängigkeit, kein neues Feature. Alles läuft in `no_std`
  ohne `alloc`; Ergebnisse variabler Anzahl gehen in Puffer des Aufrufers.
- `rust-version` steigt von 1.85 auf **1.86** (`next_up`/`next_down` für die
  Rundung nach außen).

## Entscheidungen

| Thema | Entscheidung |
|-------|--------------|
| Nicht-Konvergenz | Fehler, der die letzte Näherung, den Restfehler und die Iterationszahl enthält |
| Reelle Polynom-Nullstellen im Intervall | garantiert lückenlos für jeden Grad |
| `Interval<T>` | implementiert `Algebra`, `Ring`, `Commutative`, `Field`; Gesetze gelten als Einschluss |
| Koeffizientenreihenfolge | aufsteigend: `c[0] + c[1]·x + c[2]·x² + …` |
| Mehrfache Nullstellen | geschlossene Formeln melden jede verschiedene reelle Nullstelle einmal, aufsteigend sortiert |

## Schritt A – Nachträge

- `Real` bekommt `next_up(self)` und `next_down(self)` (nächste darstellbare
  Zahl). Für Dualzahlen wirkt das nur auf den Wert.
- `rust-version = "1.86"` im Workspace.

## Schritt B – Gemeinsame Typen der Löser (`solve`)

```rust
/// Abbruchkriterien. Standard: Toleranzen aus der Maschinengenauigkeit.
#[non_exhaustive]
pub struct Tolerance<R> { pub abs_x: R, pub rel_x: R, pub abs_f: R, pub max_iterations: u32 }
// Tolerance::default(), .with_abs_x(..), .with_rel_x(..), .with_abs_f(..), .with_max_iterations(..)

/// Ergebnis eines Lösers: Lösung, Restfehler |f(x)| bzw. ‖F(x)‖, Iterationen.
pub struct Solution<X, R> { pub x: X, pub residual: R, pub iterations: u32 }

#[non_exhaustive]
pub enum SolveError<X, R> {
    /// Nicht konvergiert; enthält den letzten Stand.
    NotConverged(Solution<X, R>),
    /// Die Funktionswerte an den Intervallenden haben kein verschiedenes Vorzeichen.
    NoBracket,
    /// Ein Zwischenwert war unendlich oder NaN; enthält den letzten Stand.
    NonFinite(Solution<X, R>),
    /// Die Jacobi-Matrix war singulär; enthält den letzten Stand.
    Singular(Solution<X, R>),
    /// Ungültige Eingabe (z. B. leeres Intervall, Leitkoeffizient null).
    InvalidInput,
    /// Der Ausgabepuffer des Aufrufers ist zu klein.
    CapacityExceeded,
}
```

Mit `Display` und `core::error::Error`.

## Schritt C – Nullstellen einer Veränderlichen (`solve`)

Alle generisch über `T: Real`; keiner braucht ein Rechen-Backend.

```rust
pub fn bisection<T, F>(f: F, a: T, b: T, tol: &Tolerance<T>) -> Result<Solution<T, T>, SolveError<T, T>>;
pub fn brent<T, F>(f: F, a: T, b: T, tol: &Tolerance<T>) -> Result<…>;      // Standardverfahren mit Klammer
pub fn secant<T, F>(f: F, x0: T, x1: T, tol: &Tolerance<T>) -> Result<…>;
pub fn newton<T, F>(f: F, x0: T, tol: &Tolerance<T>) -> Result<…>;          // F: Fn(T) -> (T, T) = (Wert, Ableitung)
pub fn newton_auto<T, F>(f: F, x0: T, tol: &Tolerance<T>) -> Result<…>;     // F: Fn(Dual<T>) -> Dual<T>
pub fn newton_bracketed<T, F>(f: F, a: T, b: T, tol: &Tolerance<T>) -> Result<…>; // Newton, fällt auf Bisektion zurück
```

- `brent` und `newton_bracketed` verlassen die Klammer nie und konvergieren
  immer, wenn die Klammer gültig ist.
- `newton_auto` ist der bequeme Weg: Funktion einmal generisch schreiben,
  Ableitung kommt aus den Dualzahlen.

## Schritt D – Intervallarithmetik (`interval`)

```rust
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interval<T> { lo: T, hi: T }   // Invariante: lo <= hi
```

| Bereich | Angebot |
|---------|---------|
| Erzeugen | `new(lo, hi) -> Option`, `point(x)`, `ENTIRE`, `hull_of(a, b)` |
| Abfragen | `lo`, `hi`, `mid`, `width`, `contains(x)`, `contains_zero`, `is_subset_of`, `intersect -> Option`, `hull` |
| Vergleiche | `certainly_lt/gt`, `possibly_lt/gt` (keine totale Ordnung) |
| Rechnen | `+ − · /` mit Rundung nach außen, `sqr`, `abs`, `powi`; Division durch ein Intervall mit Null ergibt `ENTIRE` |
| Traits | `Algebra`, `Ring`, `Commutative`, `Field` |

- Garantie: Das Ergebnis jeder Operation enthält das exakte Ergebnis für
  alle Punkte der Eingangsintervalle.
- NaN in der Eingabe führt zu `ENTIRE`, nie zu einem ungültigen Intervall.
- Kein `Real`, kein `Conjugate`: Intervalle haben keine totale Ordnung.
  Zerlegungen mit Pivotsuche laufen deshalb nicht auf Intervallen.

## Schritt E – Polynome (`poly`)

```rust
/// Polynom mit N Koeffizienten (Grad höchstens N − 1), aufsteigend.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Polynomial<T, const N: usize>(pub [T; N]);
```

| Für | Angebot |
|-----|---------|
| `Polynomial<T, N>`, `T: Algebra` | `eval(x)` (Horner), `eval_with_derivative(x)`, `derivative()` (gleiches `N`, höchster Koeffizient null), `degree() -> Option<usize>`, `+ −` Negation, `* T`, `map` (Koeffizienten umwandeln, z. B. nach `Dual`, `Complex`, `Interval`) |
| Slices | `horner(&[T], x)`, `horner_with_derivative(&[T], x)`, `clenshaw(&[T], x)` (Tschebyschow-Reihe) |

Auswertung an einer Stelle anderen Typs geht über `map`: z. B.
`p.map(Interval::point).eval(intervall)` liefert einen garantierten
Einschluss des Wertebereichs.

## Schritt F – Geschlossene Formeln bis Grad 4 (`poly`)

Brauchen `T: RealElementary`.

```rust
/// Bis zu N verschiedene reelle Nullstellen, aufsteigend, ohne Heap.
pub struct RealRoots<T, const N: usize> { /* [T; N], Anzahl */ }   // as_slice(), len(), is_empty(), iter()

pub fn solve_linear<T>(c: [T; 2]) -> RealRoots<T, 1>;
pub fn solve_quadratic<T>(c: [T; 3]) -> RealRoots<T, 2>;   // auslöschungsfreie Form
pub fn solve_cubic<T>(c: [T; 4]) -> RealRoots<T, 3>;       // trigonometrisch/Cardano, danach Newton-Politur
pub fn solve_quartic<T>(c: [T; 5]) -> RealRoots<T, 4>;     // Ferrari, danach Newton-Politur
pub fn solve_quadratic_complex<T>(c: [Complex<T>; 3]) -> [Complex<T>; 2];
```

- Ist der Leitkoeffizient null, wird der niedrigere Grad gelöst.
- Sind alle Koeffizienten null, ist das Ergebnis leer (jede Zahl wäre
  Nullstelle; das ist dokumentiert).

## Schritt G – Alle komplexen Nullstellen (`poly`)

Aberth-Ehrlich-Verfahren, für jeden Grad, Ausgabe in einen Puffer.

```rust
pub fn complex_roots<T: RealElementary>(
    coefficients: &[Complex<T>],
    roots: &mut [Complex<T>],          // Länge = Grad
    tol: &Tolerance<T>,
) -> Result<u32, SolveError<(), T>>;   // Iterationen
```

Schnell und für jeden Grad, aber **ohne Garantie** (Näherungsverfahren).

## Schritt H – Garantierte reelle Nullstellen im Intervall (`poly`)

Für jeden Grad, rundungssicher durch Intervallarithmetik, ohne Heap.

```rust
pub enum BracketKind {
    /// Genau eine Nullstelle (Polynom dort streng monoton, Vorzeichenwechsel).
    Unique,
    /// Mindestens eine mögliche Nullstelle, nicht weiter trennbar
    /// (z. B. mehrfache oder sehr dicht liegende Nullstellen).
    Unresolved,
}
pub struct RootBracket<T> { pub interval: Interval<T>, pub kind: BracketKind }

pub fn isolate_real_roots<T: Real>(
    coefficients: &[T],
    search: Interval<T>,
    min_width: T,
    out: &mut [RootBracket<T>],
) -> Result<usize, SolveError<(), T>>;            // Anzahl gefundener Klammern, aufsteigend

pub fn refine_root<T: Real>(coefficients: &[T], bracket: &RootBracket<T>, tol: &Tolerance<T>)
    -> Result<Solution<T, T>, SolveError<T, T>>;  // Brent innerhalb der Klammer
```

Verfahren: Der Wertebereich des Polynoms und seiner Ableitung wird auf
einem Teilintervall per Intervall-Horner eingeschlossen.

1. Enthält der Wertebereich die Null nicht → dort sicher keine Nullstelle.
2. Enthält der Wertebereich der Ableitung die Null nicht → streng monoton,
   also höchstens eine Nullstelle; Vorzeichen an den Enden entscheidet.
3. Sonst halbieren. Unter `min_width` wird das Stück als `Unresolved`
   gemeldet.

**Garantie:** Jede reelle Nullstelle im Suchintervall liegt in einer der
gemeldeten Klammern. Es kann zu viel gemeldet werden (`Unresolved`), nie zu
wenig. Die Halbierung nutzt einen festen Stapel auf dem Stack; ist der
Ausgabepuffer zu klein, gibt es `CapacityExceeded` statt stillem Verlust.

## Schritt I – Autodiff-Hilfen (`autodiff`)

```rust
pub fn derivative<T, F>(f: F, x: T) -> (T, T);                         // (f, f')      F: Fn(Dual<T>) -> Dual<T>
pub fn second_derivative<T, F>(f: F, x: T) -> (T, T, T);               // (f, f', f'') F: Fn(Dual<Dual<T>>) -> …
pub fn gradient<T, F, const N: usize>(f: F, x: Vector<T, N>) -> (T, Vector<T, N>);
pub fn jacobian<T, F, const N: usize, const M: usize>(f: F, x: Vector<T, N>)
    -> (Vector<T, M>, Matrix<T, M, N>);
pub fn hessian<T, F, const N: usize>(f: F, x: Vector<T, N>) -> (T, Vector<T, N>, Matrix<T, N, N>);
```

Die Funktion `f` wird jeweils mit Dualzahlen aufgerufen; man schreibt sie
einmal generisch über die Skalar-Traits.

## Schritt J – Newton für Gleichungssysteme (`solve`)

```rust
pub fn newton_system<T, F, const N: usize>(f: F, x0: Vector<T, N>, tol: &Tolerance<T>)
    -> Result<Solution<Vector<T, N>, T>, SolveError<Vector<T, N>, T>>
where F: Fn(Vector<DualN<T, N>, N>) -> Vector<DualN<T, N>, N>;

pub fn newton_system_with_jacobian<T, F, const N: usize>(f: F, x0: Vector<T, N>, tol: &Tolerance<T>)
    -> Result<…>
where F: Fn(Vector<T, N>) -> (Vector<T, N>, Matrix<T, N, N>);
```

- Jacobi-Matrix aus `DualN`, Schritt über die LU-Zerlegung aus Phase 3.
- Gedämpft: Rückwärts-Liniensuche auf `‖F‖²`, bis der Restfehler sinkt.
- Braucht `T: RealElementary` nur für die Norm; singuläre Jacobi-Matrix
  und nicht endliche Werte sind Fehler mit letztem Stand.

## Dateien

```
Cargo.toml                           rust-version 1.86
src/scalar/real.rs, primitives.rs    next_up/next_down
src/number/dual/traits.rs            next_up/next_down

src/interval.rs                      Deklarationen
src/interval/repr.rs  ops.rs  traits.rs  tests.rs

src/solve.rs                         Deklarationen
src/solve/error.rs  tolerance.rs  solution.rs
src/solve/bisection.rs  brent.rs  secant.rs  newton.rs  system.rs  tests.rs

src/poly.rs                          Deklarationen
src/poly/repr.rs  ops.rs  eval.rs  real_roots.rs  closed.rs
src/poly/aberth.rs  isolate.rs  tests.rs

src/autodiff.rs                      Deklarationen
src/autodiff/scalar.rs  vector.rs  tests.rs

src/lib.rs                           neue Module + Re-Exporte
```

Abhängigkeiten zwischen den neuen Modulen (nur nach unten):
`interval` → Skalar-Traits; `solve` → `number`, `linalg`; `autodiff` →
`number`, `linalg`; `poly` → `interval`, `solve`, `number`.

## Tests

Erwartungswerte von Hand herleitbar; Löser zusätzlich über Rückeinsetzen.

- **Skalare Löser**: `x² − 2` auf `[1, 2]` → Residuum unter Toleranz und
  `x² ≈ 2`; `x³ − x − 2`; Klammer ohne Vorzeichenwechsel → `NoBracket`;
  `max_iterations = 1` → `NotConverged` mit letztem Stand in der Klammer;
  `newton_auto` gleich `newton` mit Ableitung von Hand; Newton auf
  `x² + 1` → `NotConverged`; alle Verfahren liefern dieselbe Nullstelle.
- **Interval**: Einschluss-Eigenschaft für `+ − · /` an Beispielen, deren
  exaktes Ergebnis bekannt ist (`[1,2]·[−3,4] = [−6,8]`); `0.1 + 0.2`
  enthält den `f64`-Wert von `0.3`; Division durch `[−1, 1]` → `ENTIRE`;
  `sqr([−2, 1]) = [0, 4]` (enger als `x·x`); `new(2, 1)` → `None`.
- **Polynomial**: `1 + 2x + 3x²` bei `x = 2` → `17`, Ableitung `14`;
  `derivative`, `degree`; Auswertung mit `Dual`, `Complex` und `Interval`
  über `map`; Clenshaw gegen `T₂(x) = 2x² − 1`.
- **Geschlossene Formeln**: `(x−1)(x−2)`, `(x−1)(x−2)(x−3)`,
  `(x−1)(x−2)(x−3)(x−4)` → genau diese Nullstellen; doppelte Nullstelle
  einmal; kein reelles Ergebnis für `x² + 1`; Auslöschungsfall
  `x² − 10⁸x + 1`; Leitkoeffizient null; komplex: `x² + 1` → `±i`.
- **Aberth**: Nullstellen von `x⁴ − 1` (`±1, ±i`) und eines Polynoms
  5. Grades mit bekannten Nullstellen, jeweils über Rückeinsetzen.
- **Isolation**: `(x−1)(x−2)(x−3)` auf `[0, 4]` → drei `Unique`-Klammern,
  jede enthält genau ihre Nullstelle; Suchintervall ohne Nullstelle → 0;
  doppelte Nullstelle `(x−1)²` → `Unresolved`-Klammer, die `1` enthält;
  Polynom 7. Grades mit bekannten Nullstellen: keine fehlt; zu kleiner
  Puffer → `CapacityExceeded`; `refine_root` liefert die Nullstelle.
- **Autodiff-Hilfen**: `x³` → `(8, 12, 12)` bei `x = 2`; Gradient von
  `x²·y`; Jacobi-Matrix von `(x·y, x + y)`; Hesse-Matrix von `x²·y`.
- **Newton-Systeme**: `x² + y² = 4`, `x·y = 1` (Residuum); lineares System
  in einem Schritt; singuläre Jacobi-Matrix → `Singular` mit letztem
  Stand; beide Varianten liefern dasselbe.

Jeder Test läuft für `f32` und `f64`.

## Reihenfolge

Schritte A bis J, je ein Commit. Die KI führt nichts aus; nach der
Umsetzung ein Prüflauf durch den Menschen, Korrekturen anhand von
`pruefung/ergebnis.txt`.

## Risiken

- Ohne eigene Prüfläufe sind bei diesem Umfang mehrere Korrekturrunden zu
  erwarten.
- Die geschlossene Formel für Grad 4 ist numerisch heikel; die
  Newton-Politur fängt das ab, garantiert aber nichts. Wo es auf
  Vollständigkeit ankommt, ist Schritt H das richtige Werkzeug.
