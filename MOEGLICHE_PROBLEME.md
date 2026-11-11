# Mögliche Probleme

| Priorität | Fundstelle | Problem und Auswirkung | Stand |
|---|---|---|---|
| Hoch | [LU-Zerlegung](https://github.com/Andie2302/azlab/blob/master/crates/azmath/src/linalg/decomp/lu.rs) | Die Pivotbewertung verwendet `norm_sqr()`. Bei `f64` wird beispielsweise das Quadrat von `1e-200` zu null. Dadurch kann bereits eine invertierbare 1×1-Matrix als singulär abgelehnt werden. | **Behoben.** Pivotsuche und Singularitätsprüfung nutzen `Conjugate::abs_sum` (Summe der Beträge, ohne Quadrieren); ebenso die Rangprüfung in `Qr::solve`. Regressionstest `extremes` in `linalg/decomp/tests.rs`. |
| Hoch | [Vektornorm](https://github.com/Andie2302/azlab/blob/master/crates/azmath/src/linalg/vector/products.rs) | `sqrt(sum(x²))` kann über- oder unterlaufen, obwohl die tatsächliche Norm darstellbar ist. Für `[1e200, 0]` entsteht unendlich; `normalize()` liefert dann `None`. Eine skalierte Normberechnung wäre robuster. | **Behoben.** `Vector::norm`, `Matrix::norm` und `Quaternion::norm` teilen vor dem Quadrieren durch den größten Betrag. Regressionstest `extremes` in `linalg/vector/tests.rs`. |
| Mittel | [Bisektion](https://github.com/Andie2302/azlab/blob/master/crates/azmath/src/solve/bisection.rs) | `lo + (hi - lo) * 0.5` kann bei endlichen Grenzen wie `−1e308` und `1e308` überlaufen. Hier braucht es eine gegen Überlauf geschützte Mittelpunktberechnung. | **Behoben.** Mittelpunkt als `lo/2 + hi/2`; dieselbe Korrektur in `brent` und `newton_bracketed`. Regressionstest `extremes` in `solve/tests.rs`. |
| Mittel | [Newton](https://github.com/Andie2302/azlab/blob/master/crates/azmath/src/solve/newton.rs) | Nichtendliche Ableitungen werden nicht konsequent abgefangen. Liefert die Funktion etwa `(1.0, INFINITY)`, wird der Schritt null und das Verfahren meldet wegen der Schrittweite Erfolg — trotz Restfehler 1. | **Behoben.** `newton` und `secant` melden `NonFinite` mit dem letzten Stand; `newton_bracketed` halbiert in diesem Fall. Regressionstest `extremes` in `solve/tests.rs`. |
| Niedrig | [README](https://github.com/Andie2302/azlab/blob/master/README.md) | Mehrere Bereiche stehen noch unter „geplant“, obwohl bereits Implementierungen existieren: etwa Lie-Gruppen, iterative Löser, ODE, Quadratur und Interpolation. Auch die Aussage zum lediglich vorbereiteten `alloc`-Feature ist überholt. | **Behoben.** Die README beschreibt den tatsächlichen Umfang; ihre Beispiele laufen als Tests mit. |

## Bekannte verbleibende Grenzen

- Die Householder-Schritte der QR-Zerlegung und die Cholesky-Zerlegung
  bilden weiterhin Quadratsummen. Für Matrizen mit Einträgen nahe den
  Grenzen des Zahlbereichs können sie über- oder unterlaufen.
- `CayleyDickson::norm` ist nicht skaliert.



# Diese Punkte würde ich als Nächstes korrigieren:

**1. Rückwärts-Autodiff verwechselt Variablen verschiedener Bänder**

In [`var.rs`](https://github.com/Andie2302/azlab/blob/40e018cbdbe27e97b886571a095099fea72b9680/crates/azmath/src/autodiff/var.rs) werden nur die Knotenindizes weitergereicht. Ob zwei Variablen tatsächlich zum selben Band gehören, wird nicht geprüft.

Konkretes Beispiel:

```rust
let a = Tape::<f64, 16>::new();
let b = Tape::<f64, 16>::new();

let x = Var::new(&a, 2.0); // Index 0 auf a
let y = Var::new(&b, 3.0); // Index 0 auf b
```

`x.gradient([y])` liefert nach dem aktuellen Code `Ok([1.0])`. Dokumentiert ist für eine fremde Variable jedoch ein Fehler. Auch `x + y` vermischt die beiden Knoten mit Index 0.

Zusätzlich können alte Variablen nach `clear()` wieder als gültig erscheinen, sobald ihre Indizes neu belegt wurden. **Bandidentität und eine Generation für jeden Löschzyklus** würden beide Probleme abdecken.

**2. Die Normberechnung hat noch zwei Lücken**

In [`scaled_norm`](https://github.com/Andie2302/azlab/blob/40e018cbdbe27e97b886571a095099fea72b9680/crates/azmath/src/linalg/vector/products.rs) wird mit `scale.recip()` multipliziert. Bei einem sehr kleinen subnormalen Wert wie `1e-320_f64` läuft dieser Kehrwert über. Die Norm eines Vektors mit diesem einen Eintrag wäre aber weiterhin darstellbar.

Außerdem kann `abs_sum()` selbst überlaufen: Für `Complex::new(1e308, 1e308)` ergibt die Betragssumme unendlich, obwohl die euklidische Norm von ungefähr `1.414e308` darstellbar ist.

Hier würde ich durch die Skalierung **direkt dividieren** und als Skalierungsmaß den größten Komponentenbetrag verwenden. Die Zusage in `Conjugate::abs_sum`, dass kein Überlauf auftreten könne, stimmt derzeit nicht.

**3. Die neue SVD kann eine falsche Zerlegung als erfolgreich zurückgeben**

In [`svd.rs`](https://github.com/Andie2302/azlab/blob/40e018cbdbe27e97b886571a095099fea72b9680/crates/azmath/src/linalg/decomp/svd.rs) entstehen weiterhin unskalierte Quadratsummen und das Produkt `alpha * beta`.

Für

```rust
Matrix::new([
    [1e100_f64, 1e100],
    [0.0,       1e100],
])
```

sind `alpha` und `beta` endlich, ihr Produkt läuft jedoch über. Die Orthogonalitätsprüfung akzeptiert dadurch die nichtorthogonalen Spalten und das Verfahren meldet Konvergenz.

Das ist besonders relevant, weil anschließend Rang, Pseudoinverse und Ausgleichslösung auf dieser Zerlegung beruhen. Ich würde die Rotationsberechnung skalieren und Tests für **Rekonstruktion und Orthogonalität über verschiedene Größenordnungen** ergänzen.

**4. Newton kann nach einem Schritt weiterhin einen unendlichen Restfehler akzeptieren**

Die neue Prüfung greift vor dem Schritt. Nach der erneuten Funktionsauswertung wird aber nur `fx.is_nan()` geprüft; danach kann eine kleine Schrittweite bereits zu `Ok` führen.

Beispielsweise ergibt ein Start bei `x = 1`, mit zunächst `(fx, slope) = (1, 1e20)` und anschließend `(INFINITY, 1)`, einen winzigen Schritt und danach einen erfolgreichen Rückgabewert mit unendlichem Residuum.

Die Prüfung auf endliche Werte sollte deshalb auch **nach der Auswertung und vor der Erfolgsmeldung** stehen.

**Die Trait-Frage von vorher bleibt ebenfalls offen:** `det()` ist weiterhin für `Field + Conjugate` verfügbar, damit auch für Quaternionen. Für die üblichen Determinantenformeln fehlt hier die Voraussetzung `Commutative`.

Insbesondere die Bandzugehörigkeit sowie die SVD priorisieren. Beide können derzeit plausible Ergebnisse liefern, die mathematisch falsch sind.

## Stand zu diesen fünf Punkten

Alle fünf sind behoben, jeweils mit Regressionstest; die volle Prüfmatrix ist
grün.

| Punkt | Behebung | Regressionstest |
|---|---|---|
| 1. Bandzugehörigkeit | `Var` merkt sich die Generation seines Bands; `clear` zählt sie hoch. Verknüpfungen prüfen, dass beide Operanden zum selben Band in seinem aktuellen Stand gehören. Sonst ist das Ergebnis „vergiftet“: Der Wert stimmt, `gradient` meldet `TapeError::UnknownVariable`. `gradient` prüft zusätzlich alle Beteiligten. | `autodiff/reverse_tests.rs`, Abschnitt „Zwei Bänder“ |
| 2. Normberechnung | `scaled_norm` dividiert direkt durch den Maßstab (neu: `Conjugate::div_real`). Als Maßstab dient der größte Komponentenbetrag (`Conjugate::max_abs`, ersetzt `abs_sum`); er läuft nicht über. | `limits` in `linalg/vector/tests.rs` |
| 3. SVD | Die Matrix wird vor dem Verfahren durch ihren größten Betrag geteilt, die Singulärwerte danach zurückskaliert. | `magnitudes` in `linalg/decomp/spectral_tests.rs` (Rekonstruktion und Orthogonalität von 1e-200 bis 1e200) |
| 4. Newton | `newton` und `secant` prüfen den Funktionswert auch nach dem Schritt auf Endlichkeit. `bisection`, `brent` und `newton_bracketed` melden nie Erfolg mit unendlichem Restfehler. | `extremes` in `solve/tests.rs` |
| 5. `det()` | `Matrix::det`, `Lu::det` und `Cholesky::det` verlangen jetzt `Commutative`; für Quaternion-Matrizen gibt es sie nicht mehr. | Übersetzungszeit |

Die Zusage in der Doku von `max_abs` („läuft weder über noch unter“) gilt
jetzt, weil nicht mehr summiert wird.

