# Ziele und Nicht-Ziele

## Workspace azlab

- Ein Cargo-Workspace mit eigenständigen Crates unter `crates/`.
- `azmath` ist das Fundament. Fachliche Crates bauen darauf auf und hängen
  nie umgekehrt voneinander ab, ohne dass es hier geplant ist.
- Mögliche spätere Crates (nur zur Einordnung, noch nicht geplant):
  Stoffdaten (IAPWS-95), Strömung (Navier-Stokes), Stoffsimulation,
  Robotik, KI.

## Ziele von azmath

1. **Reine Mathematik**: rechnet, stellt Rechenmethoden und Löser bereit.
   Keine Fachdomäne, keine Einheiten, keine Stoffkonstanten.
2. **`no_std` ohne `alloc` als Normalfall**: der Kern arbeitet mit festen
   Größen (const generics) und mit Puffern, die der Aufrufer bereitstellt.
3. **Generisch über den Zahlentyp**: jeder Algorithmus läuft mit `f32`, `f64`,
   Dualzahlen, wo sinnvoll mit komplexen Zahlen, später mit `f16`/`f128` –
   ohne Änderung am Algorithmus.
4. **Austauschbares Rechen-Backend**: wer `sqrt`, `exp`, `sin` berechnet, ist
   von außen unsichtbar und per Feature wählbar.
5. **Erweiterbar ohne Brüche**: neue Zahlentypen, Speicherformen, Löser und
   Backends kommen hinzu, ohne bestehende API zu ändern
   (Regeln in [02-architektur.md](02-architektur.md)).
6. **Keine Panics**: Fehler als `Result`, `unsafe` ist verboten
   (Workspace-Lints).

## Was azmath als Grundlage leisten muss

| Anwendung | Benötigt aus azmath |
|-----------|--------------------|
| IAPWS-95 | Elementarfunktionen, Polynom-/Reihenauswertung, Autodiff (Ableitungen der Helmholtz-Energie), Nullstellen- und Newton-Löser für Systeme |
| Navier-Stokes | Vektoren/Matrizen, Gitter-/Tensor-Sichten auf Puffer, dünnbesetzte und iterative lineare Löser, Zeitintegration |
| Stoffsimulation | Vektoren, Zeitintegration (auch implizit), lineare Löser |
| Robotik | Quaternionen, Lie-Gruppen SO(3)/SE(3), Lie-Algebren mit Klammer, exp/log, Jacobi-Matrizen |
| KI | Tensoren, Matrizen, Autodiff (vorwärts und rückwärts), Optimierer |
| Physik allgemein | komplexe Zahlen, Cayley-Dickson-Algebren, ODE-Löser, Quadratur, Interpolation |

## Nicht-Ziele von azmath

- Keine Physik, keine Einheiten, keine Stoff- oder Naturkonstanten.
- Keine Ein-/Ausgabe, kein Dateiformat, keine Visualisierung.
- Kein Nightly-Rust, keine instabilen Sprachfeatures.
- Kein `unsafe`, damit vorerst keine handgeschriebenen SIMD-Intrinsics.
- Keine GPU-Anbindung, kein Multithreading im Kern.
- Kein Computeralgebra-System (keine symbolische Rechnung).
- Keine beliebig genaue Arithmetik (BigInt/BigFloat) – höchstens später als
  eigenes Crate, das die azmath-Traits implementiert.
- Kein Ersatz für BLAS/LAPACK in Sachen Höchstleistung; Korrektheit und
  Portabilität gehen vor.
- Keine Pflicht-Abhängigkeit: ohne Features hat azmath null Abhängigkeiten.
