# Phase 7, Teil 1 – Gewöhnliche Differentialgleichungen (Detailplan)

Status: **umgesetzt, noch ungeprüft**. Freigegeben vom Menschen („weiter mit
dem nächsten Schritt“). Die KI hat nichts ausgeführt; geprüft wird nur durch
`scripts/pruefen.sh` beim Menschen.

Teil 2 (Quadratur und Interpolation) folgt nach einem Prüflauf.

## Ziel

Anfangswertprobleme `y' = f(t, y)` mit fester Zahl von Gleichungen lösen:
explizit mit fester und mit gesteuerter Schrittweite, implizit für steife
Systeme. Alles ohne Heap.

## Nicht-Ziele dieses Teils

- Keine dichte Ausgabe (Interpolation zwischen den Schritten); Zwischenwerte
  gibt es über einen Beobachter nach jedem angenommenen Schritt.
- Keine Ereigniserkennung (Nullstellen während der Integration).
- Keine Mehrschrittverfahren (BDF, Adams), keine Rosenbrock-Verfahren.
- Keine Systeme variabler Größe.
- Keine symplektischen Verfahren.

## Abhängigkeiten und Features

Keine neue Abhängigkeit, kein neues Feature.

## Umfang

Rund 800 Zeilen Code und Tests in einem Commit.

## Gemeinsame Typen

```rust
pub struct OdeState<T, const N: usize> { pub t: T, pub y: Vector<T, N>, pub steps: u32, pub rejected: u32 }

#[non_exhaustive]
pub enum OdeError<T, const N: usize> {
    InvalidInput,
    NonFinite(OdeState<T, N>),      // Lösung wurde unendlich oder NaN
    TooManySteps(OdeState<T, N>),
    StepTooSmall(OdeState<T, N>),
    SolverFailed(OdeState<T, N>),   // Newton im impliziten Schritt gescheitert
}
pub type OdeResult<T, const N: usize> = Result<OdeState<T, N>, OdeError<T, N>>;
```

Fehler tragen den letzten erreichten Zustand.

## Schritt A – Feste Schrittweite (`T: Real`, kein Backend nötig)

```rust
pub fn euler_step(f, t, y, h) -> Vector<T, N>;
pub fn rk4_step(f, t, y, h) -> Vector<T, N>;
pub fn euler(f, t0, y0, t1, steps) -> OdeResult<T, N>;
pub fn rk4(f, t0, y0, t1, steps) -> OdeResult<T, N>;
```

`f: FnMut(T, Vector<T, N>) -> Vector<T, N>`.

## Schritt B – Gesteuerte Schrittweite

```rust
#[non_exhaustive]
pub struct OdeOptions<T> { pub rel_tol: T, pub abs_tol: T, pub initial_step: Option<T>, pub min_step: T, pub max_steps: u32 }

pub fn dormand_prince(f, t0, y0, t1, options) -> OdeResult<T, N>;
pub fn dormand_prince_observed(f, t0, y0, t1, options, observer) -> OdeResult<T, N>;
```

Dormand-Prince 5(4): Verfahren 5. Ordnung mit eingebettetem Fehlerschätzer
4. Ordnung. Koeffizienten nach Dormand & Prince (1980), abgeglichen mit
Wikipedia, „Dormand–Prince method“. Rückwärts integrieren (`t1 < t0`) ist
erlaubt.

## Schritt C – Implizit, für steife Systeme

```rust
pub fn implicit_euler(f, t0, y0, t1, steps, tol) -> OdeResult<T, N>;   // 1. Ordnung, sehr stabil
pub fn trapezoid(f, t0, y0, t1, steps, tol) -> OdeResult<T, N>;        // 2. Ordnung (Crank-Nicolson)
```

`f` wird mit Dualzahlen aufgerufen
(`FnMut(T, Vector<DualN<T, N>, N>) -> Vector<DualN<T, N>, N>`); die
Jacobi-Matrix für den Newton-Schritt entsteht automatisch. `tol` steuert
das Newton-Verfahren.

## Dateien

```
src/ode.rs                     Deklarationen
src/ode/error.rs  fixed.rs  adaptive.rs  implicit.rs  tests.rs
src/lib.rs                     angepasst
```

## Tests

Vergleiche mit der bekannten Lösung, die zur Laufzeit über die
Elementarfunktionen berechnet wird, und Eigenschaften der Verfahren; keine
von der KI ausgerechneten Erwartungswerte.

- `y' = −y`: Euler wird mit mehr Schritten genauer; RK4 ist genauer als
  Euler; Dormand-Prince trifft `e^(−t)` innerhalb der Toleranz.
- Schwingung `y'' = −y`: RK4 und Dormand-Prince gegen `cos`/`sin`.
- Hin und zurück integrieren ergibt den Startwert.
- Beobachter: wird je angenommenem Schritt einmal gerufen, Zeiten wachsen,
  letzter Wert ist das Ziel.
- Steifes Problem `y' = −1000·y` mit grober Schrittweite: implizites
  Euler-Verfahren bleibt beschränkt, explizites nicht.
- Trapezregel ist genauer als implizites Euler-Verfahren.
- Fehlerfälle: ungültige Eingaben, zu wenige erlaubte Schritte (Fehler trägt
  den letzten Zustand).
