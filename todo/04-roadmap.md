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
Volle Prüfmatrix grün.

- [x] Trait-Umbau: `Algebra`, `Commutative`, `Conjugate`; `from_f64` nach `Field`

- [x] `Complex<T>` inkl. `Elementary`
- [x] `Dual<T>` inkl. aller Skalar-Traits (erste Ableitung)
- [x] `DualN<T, N>` (Gradient); zweite Ableitungen über Schachtelung `Dual<Dual<T>>`
- [x] `CayleyDickson<A>`, Oktonionen als Probe
- [x] `Quaternion<T>` mit Rotationen, Slerp, Umwandlungen

## Phase 3 – Lineare Algebra fester Größe

Detailplan: [phase-3-lineare-algebra.md](phase-3-lineare-algebra.md) – umgesetzt.
Volle Prüfmatrix grün.

- [x] `Vector<T, N>`, `Matrix<T, R, C>`: Arithmetik, Skalar-/Kreuzprodukt, Normen, Transponieren
- [x] Zugriff für Algorithmen: über `IndexMut<(usize, usize)>` statt eigener Traits
- [x] Determinante, Inverse kleiner Matrizen
- [x] LU mit Pivotisierung, Cholesky, QR; Lösen linearer Systeme
- [x] Fehler-Enum der Schicht

## Phase 4 – Nullstellen, Polynome, erste Löser

Detailplan: [phase-4-nullstellen-polynome.md](phase-4-nullstellen-polynome.md) – umgesetzt.
Volle Prüfmatrix grün.

- [x] Polynomauswertung (Horner), Tschebyschow
- [x] Bisektion, Sekante, Newton, Brent
- [x] Gleichungen höherer Ordnung: Formeln Grad 2–4, alle Polynom-Nullstellen (komplex), reelle Nullstellen im Intervall
- [x] Intervallarithmetik `Interval<T>` mit Rundung nach außen (für beweisbar konservative Tests, z. B. CCD)
- [x] Newton für Systeme (Jacobi-Matrix über `DualN`), Dämpfung und Liniensuche
- [x] `derivative`, `gradient`, `jacobian` als Autodiff-Hilfen

Meilenstein: **ausreichend für ein IAPWS-95-Crate** und, zusammen mit
Phase 3, für PBD/XPBD-Gewebe mit kontinuierlicher Kollisionserkennung in
`azphysics`.

## Phase 5 – Lie-Theorie

Detailplan: [phase-5-lie-theorie.md](phase-5-lie-theorie.md) – umgesetzt. Volle Prüfmatrix grün.

- [x] Traits `LieAlgebra`, `LieGroup`
- [x] `SO2`, `SO3`, `SE2`, `SE3` mit exp/log, Adjungierte, Klammer
- [x] Matrix-Kommutator, Matrix-Exponential kleiner Matrizen

Meilenstein: **ausreichend für Robotik-Kinematik.**

## Phase 6 – Sichten, Tensoren, große Systeme

Teil 1: [phase-6a-sichten-iterative-loeser.md](phase-6a-sichten-iterative-loeser.md) – umgesetzt, volle Prüfmatrix grün.
Teil 2: [phase-6b-gmres-tensoren-alloc.md](phase-6b-gmres-tensoren-alloc.md) – umgesetzt, noch ungeprüft.

- [x] Matrix-Sichten auf Puffer des Aufrufers (Teil 1)
- [x] Tensor-Sichten auf Puffer des Aufrufers (Teil 2)
- [x] `LinearOperator`, dünnbesetzte Formate (CSR) als Sicht (Teil 1)
- [x] CG, BiCGSTAB (Teil 1)
- [x] GMRES (Teil 2)
- [x] Tensor-Grundoperationen (Kontraktion, äußeres Produkt, Umformen) (Teil 2)
- [x] Feature `alloc`: `DMatrix`, `DTensor` (Teil 2)

## Phase 7 – ODE, Quadratur, Interpolation

Umgesetzt und geprüft (volle Prüfmatrix grün). Teil 1:
[phase-7a-ode.md](phase-7a-ode.md).

- [x] Euler, RK4, Dormand-Prince 5(4) mit Schrittweitensteuerung
- [x] Implizites Euler-Verfahren, Trapezregel
- [x] Trapez-, Simpsonregel, Gauß-Legendre (Stützstellen berechnet), adaptive Simpson-Quadratur
- [x] Lineare Interpolation, natürliche kubische Splines, baryzentrische Polynominterpolation, Tschebyschow-Reihen

## Phase 8 – PDE auf Gittern

Umgesetzt und geprüft (volle Prüfmatrix grün). Modul `pde`.

- [x] Regelmäßige Gitter (`UniformGrid`), Randbedingungen (`Boundary`: Dirichlet, Neumann, periodisch)
- [x] Finite-Differenzen-Operatoren als `LinearOperator`: `Laplacian`, `Derivative` (zentral, vorwärts, rückwärts; die einseitigen Formeln sind die Aufwind- bzw. Finite-Volumen-Flüsse 1. Ordnung)
- [x] Linienmethode: explizit über `ode::rk4_slice`/`euler_slice`, implizit über `pde::theta_step`
- [x] Mehrgitter für die Poisson-Gleichung in 1D und 2D
- [x] Referenzprobleme als Tests: Poisson, Wärmeleitung, Transport, Welle
- [x] `AffineOperator` (`a·I + b·A`) in `linalg`

Nicht enthalten: Finite-Volumen-Verfahren höherer Ordnung mit Limitern,
Mehrgitter in 3D oder für andere Operatoren, nichtlineare PDE. Das gehört in
die Fach-Crates bzw. kommt bei Bedarf.

Meilenstein: **ausreichend für ein Navier-Stokes-Crate.**

## Phase 9 – Finite Elemente

Umgesetzt und geprüft (volle Prüfmatrix grün). Modul `fem`.

- [x] Netz als Knoten- und Elementliste des Aufrufers
- [x] Elementmatrizen: Strecke, Dreieck, Viereck (bilinear, 2×2 Gauß), Tetraeder; Steifigkeit, Masse, Last
- [x] Ebene lineare Elastizität auf dem Dreieck (Spannungs- und Verzerrungszustand)
- [x] `ElementOperator`: zusammengesetzte Matrix als `LinearOperator` ohne Speicherung, ein oder mehrere Freiheitsgrade je Knoten
- [x] Dirichlet-Randbedingungen über Maske und `lift`, rechte Seite über `assemble_load`
- [x] Feature `alloc`: `assemble_csr`
- [x] Referenzproblem: Poisson auf dem Einheitsquadrat mit Dreiecksnetz

Nicht enthalten: Elemente höherer Ordnung, Elastizität in 3D, Netzerzeugung.

Meilenstein: **ausreichend für Gewebe- und Festkörpermodelle auf
Finite-Elemente-Basis.**

## Phase 10 – Rückwärts-Autodiff und Optimierung

Umgesetzt und geprüft (volle Prüfmatrix grün). Module `autodiff`, `optim`.

- [x] `Tape` mit fester Kapazität (ohne Heap), `GrowTape` (Feature `alloc`)
- [x] `Var`: implementiert alle Skalar-Traits, läuft durch jeden generischen Algorithmus
- [x] `reverse_gradient` als bequemer Einstieg
- [x] Gradientenabstieg, BFGS, L-BFGS (Verlauf fester Länge, ohne Heap)
- [x] Nelder-Mead
- [x] Levenberg-Marquardt, auch mit automatischer Jacobi-Matrix

Meilenstein: **ausreichend für KI-Grundlagen.**

## Phase 11 – Ausbau

- [ ] Eigenwerte, SVD
- [ ] FFT
- [ ] Weitere Backends
- [ ] `f16`/`f128`, sobald stabil
- [ ] Spezielle Funktionen (Gamma, Bessel, erf) nach Bedarf
