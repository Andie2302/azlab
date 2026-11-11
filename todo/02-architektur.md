# Architektur von azmath

## 1. Schichten

Jede Schicht darf nur Schichten unter sich benutzen. Das hält azmath
erweiterbar und erlaubt, Schichten später in eigene Crates auszulagern.

```
6  Löser & Numerik    Nullstellen, lineare/nichtlineare Systeme, ODE,
                      Quadratur, Interpolation, Optimierung
5  Autodiff           vorwärts (Dualzahlen), rückwärts (Tape)
4  Lie-Theorie        Lie-Algebren (Klammer), Lie-Gruppen (exp/log)
3  Lineare Algebra    Vektor, Matrix, Zerlegungen, Tensoren
2  Zahlentypen        Complex, Dual, Cayley-Dickson, Quaternion
1  Skalar-Traits      Ring, Field, Real, Elementary, ...
0  Backend (privat)   libm | std | weitere
```

## 2. Skalar-Traits (Schicht 1)

Kleine Stufen, damit jeder Zahlentyp nur zusagt, was er kann, und jeder
Algorithmus nur fordert, was er braucht.

| Trait | Zusage | Typen |
|-------|--------|-------|
| `Ring` | `+ - *`, `ZERO`, `ONE`; Multiplikation nicht zwingend kommutativ | Ganzzahlen, Gleitkomma, Complex, Quaternion |
| `Field: Ring` | Division, `recip`, `powi` | Gleitkomma, Complex, Dual, Quaternion (Schiefkörper) |
| `Elementary: Field` | `exp`, `ln`, `sqrt`, `sin`, `cos`, `tan`, Umkehr- und Hyperbelfunktionen, `powf` | Gleitkomma, Complex, Dual |
| `Real: Field + PartialOrd` | Ordnung, `abs`, `min`/`max`, Konstanten (`PI`, `E`, `EPSILON`), `from_f64` | Gleitkomma, Dual über Real |
| `RealElementary: Real + Elementary` | `atan2`, `hypot`, `cbrt`, `floor`/`ceil`/`round`/`trunc`, `mul_add` | Gleitkomma, Dual über Real |
| `Conjugate` | Konjugation, Normquadrat | Complex, Quaternion, Cayley-Dickson; für reelle Typen die Identität |

Wichtig: `Elementary` setzt **keine** Ordnung voraus, damit komplexe Zahlen
es implementieren können. Nicht-assoziative Algebren (Oktonionen) erfüllen
`Ring` nicht und bekommen nur die Operatoren plus `Conjugate`.

`Real`, `Ring`, `Field` brauchen kein Backend. Nur die Implementierungen von
`Elementary`/`RealElementary` für `f32`/`f64` kommen aus dem Backend. Ohne
Backend-Feature fehlen genau diese Implementierungen, alles andere baut
(Regel 1).

## 3. Backend (Schicht 0)

- Privates Modul `backend`. Jede Backend-Datei enthält nichts als die
  Implementierung von `Elementary`/`RealElementary` für `f32` und `f64`.
- Auswahl per Feature, Rangfolge fest: `std` > `libm` > (weitere). Es wird
  immer genau eine Backend-Datei übersetzt, Features bleiben additiv.
- Ein Tausch ändert keine öffentliche Signatur. Aufbauende Crates reichen die
  Features nur durch.
- Neues Backend = neue Datei + neues Feature + Eintrag in der Rangfolge.
  Kandidaten: schnelle Näherungen für Mikrocontroller, später eine eigene
  Implementierung (nötig für `f16`/`f128`, sobald stabil).

## 4. Zahlentypen (Schicht 2)

- `Complex<T>`: `Field`, `Elementary`, `Conjugate`.
- `Dual<T>`: Wert + Ableitung; implementiert alle Skalar-Traits, die `T` hat.
  Dadurch liefert jeder generische Algorithmus automatisch Ableitungen.
  `DualN<T, const N: usize>` für Gradienten ohne Heap; `HyperDual<T>` bzw.
  Schachtelung `Dual<Dual<T>>` für zweite Ableitungen.
- `CayleyDickson<A>`: generische Verdopplung (reell → komplex → Quaternion →
  Oktonion → ...) über `Conjugate`.
- `Quaternion<T>`: eigener Typ mit Rotations-API (Robotik), umwandelbar
  von/nach `CayleyDickson<Complex<T>>`.

## 5. Lineare Algebra und Tensoren (Schicht 3)

- Feste Größen: `Vector<T, const N>`, `Matrix<T, const R, const C>` auf dem
  Stack. Auf stable sind Größenrechnungen im Typ (`R + C`) nicht möglich;
  solche Operationen werden nicht oder mit explizit angegebener Zielgröße
  angeboten.
- Große/variable Größen ohne Heap: Sichten auf Puffer des Aufrufers
  (`MatrixView`, `MatrixViewMut`, Tensor-Sicht mit `[usize; RANK]` als Form).
- Mit Feature `alloc`: besitzende dynamische Typen (`DMatrix`, `DTensor`).
- **Algorithmen sind gegen Zugriffs-Traits geschrieben** (lesender und
  schreibender Matrixzugriff), nicht gegen konkrete Typen. LU, QR, Cholesky
  usw. existieren genau einmal und laufen auf festen, geliehenen und
  dynamischen Matrizen.
- Dünnbesetzte Matrizen und matrixfreie Verfahren über einen Trait
  `LinearOperator` (Anwenden auf einen Vektor).

## 6. Lie-Theorie (Schicht 4)

- Trait `LieAlgebra`: Vektorraum + `bracket`, `hat`/`vee`.
- Trait `LieGroup`: `compose`, `inverse`, `exp`, `log`, `adjoint`, Wirkung.
- Konkret: `SO2`, `SO3` (über Quaternion), `SE2`, `SE3`; allgemeiner
  Matrix-Kommutator für quadratische Matrizen.

## 7. Autodiff (Schicht 5)

- Vorwärts: über `Dual`/`DualN`, ohne Heap. Hilfsfunktionen `derivative`,
  `gradient`, `jacobian`.
- Rückwärts: Tape. Variante mit fester Kapazität (ohne Heap) und wachsende
  Variante hinter `alloc`. Der Tape-Skalar implementiert dieselben
  Skalar-Traits.

## 8. Löser und Numerik (Schicht 6)

Alle generisch über die Skalar-Traits, Arbeitsspeicher kommt vom Aufrufer
oder liegt in const-generischen Größen.

- Nullstellen 1D: Bisektion, Sekante, Newton, Brent.
- Nichtlineare Systeme: Newton (Jacobi-Matrix über Autodiff), gedämpft.
- Lineare Systeme direkt: LU, Cholesky, QR; später Eigenwerte, SVD.
- Lineare Systeme iterativ: CG, BiCGSTAB, GMRES über `LinearOperator`.
- ODE: Euler, RK4, eingebettete RK (Schrittweitensteuerung), implizite
  Verfahren für steife Systeme.
- Quadratur: Gauß-Legendre, adaptiv.
- Interpolation: linear, kubische Splines, Polynome (Horner, Tschebyschow).
- Optimierung: Gradientenabstieg, (L-)BFGS, Levenberg-Marquardt, Nelder-Mead.

## 9. Fehlerbehandlung

- Keine Panics in Bibliotheksfunktionen. Indexzugriffe prüfend (`get`) oder
  durch const generics zur Übersetzungszeit abgesichert.
- Fehler-Enums pro Schicht, `#[non_exhaustive]`, implementieren
  `core::error::Error`.
- Iterative Verfahren geben ein Ergebnis mit Iterationszahl und Restfehler
  zurück; Nicht-Konvergenz ist ein Fehlerwert.

## 10. Dateistruktur (ohne `mod.rs`)

```
crates/azmath/src/
  lib.rs              nur Crate-Doku, #![no_std], mod-Deklarationen, pub use
  backend.rs          Auswahl per cfg
  backend/libm.rs
  backend/std.rs
  scalar.rs           Deklaration der Untermodule + Re-Exporte
  scalar/ring.rs  field.rs  real.rs  elementary.rs  conjugate.rs  primitives.rs
  number.rs
  number/complex.rs  dual.rs  cayley_dickson.rs  quaternion.rs
  linalg.rs
  linalg/vector.rs  matrix.rs  view.rs  access.rs  operator.rs  decomp.rs
  linalg/decomp/lu.rs  cholesky.rs  qr.rs
  tensor.rs  tensor/...
  lie.rs     lie/...
  autodiff.rs  autodiff/...
  solve.rs   solve/roots.rs  nonlinear.rs  iterative.rs
  ode.rs  quad.rs  interp.rs  optim.rs  poly.rs
  error.rs
```

Modul-Dateien wie `scalar.rs` enthalten ebenfalls nur Deklarationen und
Re-Exporte; Code steht in den Blatt-Dateien.

## 11. Regeln für Erweiterbarkeit

1. Algorithmen nennen nie `f32`/`f64`, nur Traits. Konstanten über
   `T::from_f64` oder Trait-Konstanten.
2. Neuer Zahlentyp = Traits implementieren, sonst nichts anfassen.
3. Neue Speicherform = Zugriffs-Traits implementieren.
4. Features sind additiv; kein Feature entfernt API oder ändert Verhalten
   außer Genauigkeit/Geschwindigkeit des Backends.
5. Alles hinter `alloc`/`std` ist Zusatz. Jede Fähigkeit hat zuerst eine
   Variante ohne Heap.
6. Öffentliche Enums und Konfigurations-Structs sind `#[non_exhaustive]`.
7. Traits bekommen neue Methoden nur mit Standardimplementierung.
8. Abhängigkeiten zeigen nur nach unten (Abschnitt 1).
9. Fachwissen (Physik, Stoffdaten) gehört nie in azmath.
