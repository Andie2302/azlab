# Roadmap azmath

Reihenfolge folgt den Schichten aus [02-architektur.md](02-architektur.md).
Jede Phase wird vor der Umsetzung hier verfeinert (API, Dateien, Tests) und
freigegeben. Die Phasen 0–2 sind bewusst klein: sie legen die Verträge fest,
auf denen alles Weitere steht.

## Phase 0 – Workspace und Gerüst

- [ ] `Cargo.toml` (Workspace): `[workspace.dependencies]`, `version`, `rust-version`
- [ ] `crates/azmath/Cargo.toml`: Features `libm`, `std`, `alloc`; `[lints] workspace = true`
- [ ] `lib.rs`: Crate-Doku, `#![no_std]`, Modul-Deklarationen (siehe Frage F2)
- [ ] Regeln aus `lib.rs` zusätzlich in `todo/README.md` (erledigt) – Verbleib im Quelltext klären
- [ ] Prüfmatrix als Skript oder CI (siehe 05-qualitaet.md, Frage F6)

## Phase 1 – Skalar-Traits und Backend

- [ ] `Ring`, `Field` (`recip`, `powi`)
- [ ] `Real` (Ordnung, `abs`, Konstanten, `from_f64`)
- [ ] `Elementary`, `RealElementary`
- [ ] `Conjugate`
- [ ] Implementierungen für `f32`, `f64`; `Ring` für vorzeichenbehaftete Ganzzahlen
- [ ] Backend `libm`, Backend `std`, Auswahl mit Rangfolge
- [ ] Tests: generische Testfunktionen, die für jeden Typ und jedes Backend laufen

Ergebnis: Vertrag, gegen den alle späteren Algorithmen geschrieben werden.

## Phase 2 – Zahlentypen

- [ ] `Complex<T>` inkl. `Elementary`
- [ ] `Dual<T>` inkl. aller Skalar-Traits (erste Ableitung)
- [ ] `DualN<T, N>` (Gradient), zweite Ableitungen (`HyperDual` oder Schachtelung)
- [ ] `CayleyDickson<A>`, Oktonionen als Probe
- [ ] `Quaternion<T>` mit Rotationen, Slerp, Umwandlungen

## Phase 3 – Lineare Algebra fester Größe

- [ ] `Vector<T, N>`, `Matrix<T, R, C>`: Arithmetik, Skalar-/Kreuzprodukt, Normen, Transponieren
- [ ] Zugriffs-Traits für Matrizen (lesend/schreibend)
- [ ] Determinante, Inverse kleiner Matrizen
- [ ] LU mit Pivotisierung, Cholesky, QR; Lösen linearer Systeme
- [ ] Fehler-Enum der Schicht

## Phase 4 – Nullstellen, Polynome, erste Löser

- [ ] Polynomauswertung (Horner), Tschebyschow
- [ ] Bisektion, Sekante, Newton, Brent
- [ ] Newton für Systeme (Jacobi-Matrix über `DualN`)
- [ ] `derivative`, `gradient`, `jacobian` als Autodiff-Hilfen

Meilenstein: **ausreichend für ein IAPWS-95-Crate.**

## Phase 5 – Lie-Theorie

- [ ] Traits `LieAlgebra`, `LieGroup`
- [ ] `SO2`, `SO3`, `SE2`, `SE3` mit exp/log, Adjungierte, Klammer
- [ ] Matrix-Kommutator, Matrix-Exponential kleiner Matrizen

Meilenstein: **ausreichend für Robotik-Kinematik.**

## Phase 6 – Sichten, Tensoren, große Systeme

- [ ] Matrix- und Tensor-Sichten auf Puffer des Aufrufers
- [ ] `LinearOperator`, dünnbesetzte Formate (CSR) als Sicht
- [ ] CG, BiCGSTAB, GMRES
- [ ] Tensor-Grundoperationen (Kontraktion, äußeres Produkt, Umformen)
- [ ] Feature `alloc`: `DMatrix`, `DTensor`

## Phase 7 – ODE, Quadratur, Interpolation

- [ ] Euler, RK4, eingebettete RK mit Schrittweitensteuerung
- [ ] Implizite Verfahren für steife Systeme
- [ ] Gauß-Legendre, adaptive Quadratur
- [ ] Lineare Interpolation, kubische Splines

Meilenstein: **ausreichend für Navier-Stokes- und Stoffsimulations-Crates.**

## Phase 8 – Rückwärts-Autodiff und Optimierung

- [ ] Tape mit fester Kapazität, wachsendes Tape (`alloc`)
- [ ] Gradientenabstieg, (L-)BFGS, Levenberg-Marquardt, Nelder-Mead

Meilenstein: **ausreichend für KI-Grundlagen.**

## Phase 9 – Ausbau

- [ ] Eigenwerte, SVD
- [ ] FFT
- [ ] Weitere Backends
- [ ] `f16`/`f128`, sobald stabil
- [ ] Spezielle Funktionen (Gamma, Bessel, erf) nach Bedarf
