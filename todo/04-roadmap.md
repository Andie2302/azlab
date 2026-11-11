# Roadmap azmath

Reihenfolge folgt den Schichten aus [02-architektur.md](02-architektur.md).
Jede Phase wird vor der Umsetzung hier verfeinert (API, Dateien, Tests) und
freigegeben. Die Phasen 0–2 sind bewusst klein: sie legen die Verträge fest,
auf denen alles Weitere steht.

## Phase 0 – Workspace und Gerüst

- [x] `Cargo.toml` (Workspace): `[workspace.dependencies]`, `version`, `rust-version`
- [x] `crates/azmath/Cargo.toml`: Features `libm`, `std`, `alloc`; `[lints] workspace = true`
- [x] `lib.rs`: Crate-Doku, `#![no_std]`, Modul-Deklarationen (Frage F2 entschieden)
- [x] Regeln in `CLAUDE.md` und `todo/README.md`
- [x] Prüfmatrix als Skript `scripts/pruefen.sh` (ungetestet, erster Lauf durch den Menschen)

## Phase 1 – Skalar-Traits und Backend

Stand: als Entwurf im Code vorhanden, **vom Menschen noch zu prüfen**
(Prüfmatrix in 05-qualitaet.md).

- [x] `Ring`, `Field` (`recip`, `powi`)
- [x] `Real` (Ordnung, `abs`, Konstanten, `from_f64`)
- [x] `Elementary`, `RealElementary`
- [x] `Conjugate` (mit Phase 2 umgesetzt)
- [x] Implementierungen für `f32`, `f64`; `Ring` für vorzeichenbehaftete Ganzzahlen
- [x] Backend `libm`, Backend `std`, Auswahl mit Rangfolge
- [x] Tests: generische Testfunktionen, die für jeden Typ und jedes Backend laufen

Ergebnis: Vertrag, gegen den alle späteren Algorithmen geschrieben werden.

## Phase 2 – Zahlentypen

Detailplan: [phase-2-zahlentypen.md](phase-2-zahlentypen.md) – umgesetzt.
Von der KI mit Standard-Features und ohne Features geprüft; die volle
Prüfmatrix (`scripts/pruefen.sh`) steht noch aus.

- [x] Trait-Umbau: `Algebra`, `Commutative`, `Conjugate`; `from_f64` nach `Field`

- [x] `Complex<T>` inkl. `Elementary`
- [x] `Dual<T>` inkl. aller Skalar-Traits (erste Ableitung)
- [x] `DualN<T, N>` (Gradient); zweite Ableitungen über Schachtelung `Dual<Dual<T>>`
- [x] `CayleyDickson<A>`, Oktonionen als Probe
- [x] `Quaternion<T>` mit Rotationen, Slerp, Umwandlungen

## Phase 3 – Lineare Algebra fester Größe

- [ ] `Vector<T, N>`, `Matrix<T, R, C>`: Arithmetik, Skalar-/Kreuzprodukt, Normen, Transponieren
- [ ] Zugriffs-Traits für Matrizen (lesend/schreibend)
- [ ] Determinante, Inverse kleiner Matrizen
- [ ] LU mit Pivotisierung, Cholesky, QR; Lösen linearer Systeme
- [ ] Fehler-Enum der Schicht

## Phase 4 – Nullstellen, Polynome, erste Löser

- [ ] Polynomauswertung (Horner), Tschebyschow
- [ ] Bisektion, Sekante, Newton, Brent
- [ ] Gleichungen höherer Ordnung: Formeln Grad 2–4, alle Polynom-Nullstellen (komplex), reelle Nullstellen im Intervall
- [ ] Intervallarithmetik `Interval<T>` mit Rundung nach außen (für beweisbar konservative Tests, z. B. CCD)
- [ ] Newton für Systeme (Jacobi-Matrix über `DualN`), Dämpfung und Liniensuche
- [ ] `derivative`, `gradient`, `jacobian` als Autodiff-Hilfen

Meilenstein: **ausreichend für ein IAPWS-95-Crate** und, zusammen mit
Phase 3, für PBD/XPBD-Gewebe mit kontinuierlicher Kollisionserkennung in
`azphysics`.

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

## Phase 8 – PDE auf Gittern

- [ ] Regelmäßige Gitter als Sicht, Randbedingungen
- [ ] Finite-Differenzen- und Finite-Volumen-Operatoren als `LinearOperator`
- [ ] Linienmethode (Kopplung an die ODE-Löser)
- [ ] Mehrgitter
- [ ] Referenzprobleme: Poisson, Wärmeleitung, Advektion, Welle

Meilenstein: **ausreichend für ein Navier-Stokes-Crate.**

## Phase 9 – Finite Elemente

- [ ] Netz-Sicht (Knoten, Elemente) auf Puffer des Aufrufers
- [ ] Referenzelemente, Ansatzfunktionen, Quadraturregeln
- [ ] Assemblierung in dünnbesetzte Matrizen, Randbedingungen
- [ ] Referenzprobleme: Poisson, lineare Elastizität

Meilenstein: **ausreichend für Gewebe- und Festkörpermodelle auf
Finite-Elemente-Basis.**

## Phase 10 – Rückwärts-Autodiff und Optimierung

- [ ] Tape mit fester Kapazität, wachsendes Tape (`alloc`)
- [ ] Gradientenabstieg, (L-)BFGS, Levenberg-Marquardt, Nelder-Mead

Meilenstein: **ausreichend für KI-Grundlagen.**

## Phase 11 – Ausbau

- [ ] Eigenwerte, SVD
- [ ] FFT
- [ ] Weitere Backends
- [ ] `f16`/`f128`, sobald stabil
- [ ] Spezielle Funktionen (Gamma, Bessel, erf) nach Bedarf
