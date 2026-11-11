# Mögliche Probleme

| Priorität | Fundstelle | Problem und Auswirkung | Stand |
|---|---|---|---|
| Hoch | [LU-Zerlegung](https://github.com/Andie2302/azlab/blob/master/crates/azmath/src/linalg/decomp/lu.rs) | Die Pivotbewertung verwendet `norm_sqr()`. Bei `f64` wird beispielsweise das Quadrat von `1e-200` zu null. Dadurch kann bereits eine invertierbare 1×1-Matrix als singulär abgelehnt werden. | **Behoben.** Pivotsuche und Singularitätsprüfung nutzen `Conjugate::abs_sum` (Summe der Beträge, ohne Quadrieren); ebenso die Rangprüfung in `Qr::solve`. Regressionstest `extremes` in `linalg/decomp/tests.rs`. |
| Hoch | [Vektornorm](https://github.com/Andie2302/azlab/blob/master/crates/azmath/src/linalg/vector/products.rs) | `sqrt(sum(x²))` kann über- oder unterlaufen, obwohl die tatsächliche Norm darstellbar ist. Für `[1e200, 0]` entsteht unendlich; `normalize()` liefert dann `None`. Eine skalierte Normberechnung wäre robuster. | **Behoben.** `Vector::norm`, `Matrix::norm` und `Quaternion::norm` teilen vor dem Quadrieren durch den größten Betrag. Regressionstest `extremes` in `linalg/vector/tests.rs`. |
| Mittel | [Bisektion](https://github.com/Andie2302/azlab/blob/master/crates/azmath/src/solve/bisection.rs) | `lo + (hi - lo) * 0.5` kann bei endlichen Grenzen wie `−1e308` und `1e308` überlaufen. Hier braucht es eine gegen Überlauf geschützte Mittelpunktberechnung. | **Behoben.** Mittelpunkt als `lo/2 + hi/2`; dieselbe Korrektur in `brent` und `newton_bracketed`. Regressionstest `extremes` in `solve/tests.rs`. |
| Mittel | [Newton](https://github.com/Andie2302/azlab/blob/master/crates/azmath/src/solve/newton.rs) | Nichtendliche Ableitungen werden nicht konsequent abgefangen. Liefert die Funktion etwa `(1.0, INFINITY)`, wird der Schritt null und das Verfahren meldet wegen der Schrittweite Erfolg — trotz Restfehler 1. | **Behoben.** `newton` und `secant` melden `NonFinite` mit dem letzten Stand; `newton_bracketed` halbiert in diesem Fall. Regressionstest `extremes` in `solve/tests.rs`. |
| Niedrig | [README](https://github.com/Andie2302/azlab/blob/master/README.md) | Mehrere Bereiche stehen noch unter „geplant“, obwohl bereits Implementierungen existieren: etwa Lie-Gruppen, iterative Löser, ODE, Quadratur und Interpolation. Auch die Aussage zum lediglich vorbereiteten `alloc`-Feature ist überholt. | Offen; wird mit dem Abschluss von azmath überarbeitet. |

## Bekannte verbleibende Grenzen

- Die Householder-Schritte der QR-Zerlegung und die Cholesky-Zerlegung
  bilden weiterhin Quadratsummen. Für Matrizen mit Einträgen nahe den
  Grenzen des Zahlbereichs können sie über- oder unterlaufen.
- `CayleyDickson::norm` ist nicht skaliert.
