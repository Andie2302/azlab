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


# Ein spontaner Reviewer hat gemeldet:

**Ja, die zuletzt gemeldeten Probleme sind im neuen Stand gezielt korrigiert.** Bandidentität und Löschgenerationen werden geprüft, Normen dividieren direkt durch den Maßstab, Newton prüft nach dem Schritt auf Endlichkeit, und `det()` verlangt jetzt `Commutative`.

Ich habe die Änderungen bis [`75b06a8`](https://github.com/Andie2302/azlab/commit/75b06a892a55245feca742bfa35a38240014ea53) gelesen. Das eingecheckte Prüfprotokoll zeigt mit allen Features **118 bestandene Tests und fünf Dokumentationstests**, dazu erfolgreiche Clippy- und Embedded-Builds. Eigene Testläufe habe ich nicht ausgeführt; die folgenden Befunde stammen aus der Codeprüfung.

**Vier konkrete Lücken sehe ich noch:**

**1. LU: Der vorab berechnete Kehrwert kann weiterhin überlaufen**

In [`Lu::solve`](https://github.com/Andie2302/azlab/blob/75b06a892a55245feca742bfa35a38240014ea53/crates/azmath/src/linalg/decomp/lu.rs) steht:

```rust
x[i] = row[i].recip() * rest;
```

Für dieses System ist die Lösung genau 1:

```rust
let a = Matrix::new([[1e-320_f64]]);
let b = Vector::new([1e-320_f64]);
```

Der Kehrwert des Diagonalelements wird jedoch unendlich, wodurch auch die berechnete Lösung unendlich wird. Dasselbe Muster steht in `MatrixView::lu_solve`.

Für kommutative Skalare hilft direkte Division. Für Quaternionen muss dabei die Reihenfolge erhalten bleiben: Gesucht ist `diagonal⁻¹ · rest`, weshalb ein allgemeiner Austausch gegen `rest / diagonal` nicht korrekt wäre.

**2. Komplexe Division: `z / z` kann NaN ergeben**

In der [komplexen Division](https://github.com/Andie2302/azlab/blob/75b06a892a55245feca742bfa35a38240014ea53/crates/azmath/src/number/complex/ops.rs) kann trotz Smith-Verfahren der Nenner überlaufen:

```rust
let z = Complex::new(1e308_f64, 1e308);
let result = z / z; // Erwartet: 1 + 0i
```

Hier werden sowohl der Nenner als auch der Zähler des Realteils unendlich. Das ergibt NaN. Die Dokumentationszusage „vermeidet Überlauf, solange das Ergebnis darstellbar ist“ ist daher zu stark. Es braucht zusätzliche Skalierung.

**3. SVD: Einheitliches Skalieren schützt nicht vor stark unterschiedlichen Spaltengrößen**

Deine Korrektur behebt das vorherige Gegenbeispiel. Bei dieser Matrix bleiben allerdings Quadratsummen kleiner Spalten problematisch:

```rust
Matrix::new([
    [1.0_f64, 0.0,    0.0],
    [0.0,     1e-200, 1e-200],
    [0.0,     0.0,    1e-200],
])
```

Der globale Maßstab ist 1. Für die beiden kleinen Spalten werden `alpha`, `beta` und `gamma` durch Unterlauf null. Die [SVD](https://github.com/Andie2302/azlab/blob/75b06a892a55245feca742bfa35a38240014ea53/crates/azmath/src/linalg/decomp/svd.rs) überspringt deren Rotation, obwohl sie nicht orthogonal sind. Entsprechend stimmen die kleinen Singulärwerte und die zugehörigen linken Singulärvektoren nicht.

Hier helfen Skalierung innerhalb jedes Spaltenpaars und eine geschützte Korrelationsberechnung. Die bisherigen Tests verändern den **gemeinsamen Maßstab** einer Matrix; ergänzen würde ich Tests mit **unterschiedlichen Größenordnungen innerhalb derselben Matrix**.

**4. Spezielle Funktionen: Große endliche Argumente können NaN liefern**

[`erf` und `erfc`](https://github.com/Andie2302/azlab/blob/75b06a892a55245feca742bfa35a38240014ea53/crates/azmath/src/special/error_function.rs) bilden zunächst `x * x`. Bei `1e200_f64` entsteht unendlich. Die anschließende Gamma-Berechnung enthält dann unter anderem `∞ − ∞` und liefert NaN.

Für diese Eingabe sollten die gerundeten Ergebnisse sein:

```rust
erf(1e200_f64)  // 1.0
erfc(1e200_f64) // 0.0
```

Die Grenzwerte für unendliche Argumente sollten in `gamma_p`/`gamma_q` ausdrücklich behandelt werden.

**Zwei weitere Verbesserungen lohnen sich:**

- **Mindestversion prüfen:** Du deklarierst Rust 1.86, das Prüfprotokoll verwendet 1.97. Ergänze einen Build mit 1.86, damit die zugesagte Mindestversion tatsächlich abgesichert ist.
- **Unvollständige Gammafunktionen:** Reihe und Kettenbruch enden spätestens nach 500 Schritten, geben aber auch ohne erreichte Konvergenz kommentarlos einen Wert zurück. Hier wäre ein erkennbarer Fehler oder eine dokumentierte Einschränkung sinnvoll.

Die bisherigen Reparaturen würde ich beibehalten. **Als Nächstes würde ich LU und komplexe Division korrigieren und anschließend die Tests um gemischte Größenordnungen erweitern.** Damit prüfst du eine andere Fehlerklasse als mit dem bloßen Vergrößern oder Verkleinern einer ganzen Matrix.

## Stand zu den Befunden des Reviewers

Alle sechs Punkte sind bearbeitet, jeweils mit Regressionstest; die volle
Prüfmatrix ist grün.

| Punkt | Behebung | Regressionstest |
|---|---|---|
| 1. LU: Kehrwert des Diagonalelements | Neue Hilfsfunktion `left_divide` berechnet `d⁻¹ · rest`, indem sie beide Seiten zuerst durch den Betrag von `d` teilt; die Reihenfolge bleibt erhalten, es stimmt also auch für Quaternionen. Verwendet in `Lu::solve`, `MatrixView::lu_solve`, `Qr::solve` und `Cholesky::solve`. | `tiny_pivots` in `linalg/decomp/tests.rs` (reell, komplex, Quaternion, Sicht auf Puffer) |
| 2. Komplexe Division | Zähler und Nenner werden je auf die Größenordnung eins gebracht, dann nach Smith geteilt; die Maßstäbe kommen am Ende dazu. `z / z` ist auch am Rand des Zahlbereichs eins. Die Doku nennt die verbleibende Grenze: Nur ein Ergebnis, das der größten darstellbaren Zahl bis auf den Faktor √2 nahe kommt, kann noch vorzeitig überlaufen. | `division_at_the_limits` in `number/complex/tests.rs` |
| 3. SVD: verschiedene Spaltengrößen | Je Spaltenpaar werden die Längen skaliert berechnet und der Kosinus des Winkels aus den normierten Spalten; die Drehung folgt aus dem Verhältnis der Längen. Es werden keine Quadratsummen der Spalten mehr gebildet. | `mixed_magnitudes` in `linalg/decomp/spectral_tests.rs` (große Spalte vorn, in der Mitte, hinten) |
| 4. `erf`, `erfc` für große Argumente | `gamma_p` und `gamma_q` behandeln unendliche Argumente ausdrücklich (`P(a, ∞) = 1`, `Q(a, ∞) = 0`); damit liefern `erf` und `erfc` die Grenzwerte, sobald `x²` nicht mehr darstellbar ist. | `special/tests.rs` |
| 5. Mindestversion | `scripts/pruefen.sh` übersetzt azmath und azphysics zusätzlich mit der im Manifest genannten Mindestversion (derzeit 1.86), sofern sie installiert ist: `rustup toolchain install 1.86 --profile minimal`. Der Stand besteht das. | Prüfprotokoll |
| 6. Unvollständige Gammafunktionen | Reihe und Kettenbruch dürfen jetzt 10 000 Glieder brauchen (reicht bis `a` um 10⁶). Ohne Konvergenz ist das Ergebnis NaN statt eines ungenauen Werts; ebenso, wenn der Zahltyp `a + 1` nicht mehr von `a` unterscheidet. Beides steht in der Doku. | `special/tests.rs` |

Nebenbei gefunden und behoben (azphysics): Die kontinuierliche
Kollisionserkennung meldete für ein schnell fliegendes ebenes Tuch
Berührungen mit sich selbst. Der zugehörige Test bestand, ohne dass sich das
Tuch bewegte; er prüft jetzt, dass es am Hindernis ankommt.

## Teilchenflüssigkeit: Wasser kam an den Wänden nicht zur Ruhe

**Behoben.** In `fluid_dam_break` bekam das Wasser an beiden Rändern des
Beckens dauerhaft „Stöße“. Schon ein ruhendes, eben gefülltes Becken zeigte
es: Die Spalte direkt an der Wand sank über die ganze Höhe mit etwa 0,3 m/s
ab, oben floss Wasser nach – eine Umwälzung an jeder Wand, die kein
Parameter beseitigte.

Ursache: Die Wandteilchen spiegelten nur den Druck des jeweiligen
Flüssigkeitsteilchens. In der Senkrechten heben sich ihre Beiträge dann auf;
die Wand trägt kein Druckgefälle mit der Tiefe. Dem Teilchen an der Wand
fehlt damit die Stütze gegen die Schwerkraft von der ganzen Wandseite.

Behebung in `fluid/particles.rs`: Wandteilchen bekommen eine eigene Dichte
und einen eigenen Druckfaktor, der mit der Tiefe wächst wie in der
Flüssigkeit; die Wand ist dafür drei statt zwei Lagen dick.

| | vorher | nachher |
|---|---|---|
| Ruhendes Becken, mittlere Geschwindigkeit an der Wand nach 20 s | 0,20 m/s | 0,000 m/s |
| Dammbruch, mittlere Geschwindigkeit an den Wänden nach 50 s | 0,08 bis 0,13 m/s | 0,002 m/s |
| Dammbruch, schnellstes Teilchen nach 50 s | 0,5 m/s | 0,06 m/s |

Regressionstest: `fluid::tests::pool` prüft die Ruhe getrennt an beiden
Wänden und in der Mitte.

## Gewebe: frei fliegendes Tuch trieb in x und z ab

**Behoben.** Ein Tuch ohne Schwerkraft und Dämpfung, das mit
(0,15; 0; 0,05) m/s fliegt, änderte stetig seine Geschwindigkeit in x und z;
y blieb exakt. Mit `f64` war der Fehler am kleinsten.

Die Ursache liegt nicht in der Kollisionserkennung, der Abstoßung oder den
Verbindungen, sondern im Zeitschritt selbst – der Fehler trat auch ganz ohne
Korrekturen auf (`iterations = 0`, CCD und Selbstkollision aus):

1. **Die Lage verliert in jedem Teilschritt denselben Rest.** Die Lage ist
   auf das Raster des Zahltyps gerundet; bei 9 m ist ein Rasterschritt in
   `f32` knapp 10⁻⁶ m, der Weg eines Teilschritts aber nur 3·10⁻⁴ m. Der
   abgeschnittene Rest ist jedes Mal gleich, also ein Dauerfehler.
2. **Die Geschwindigkeit wurde aus dem gerundeten Weg zurückgerechnet** und
   rastete damit auf Vielfache von Rasterschritt durch Teilschritt ein: bei
   9 m sind das 4,6·10⁻⁴ m/s – genau die gemessene Größe.
3. **Nachbarn werden verschieden gerundet.** Die Verbindungen sahen eine
   Scheinverzerrung von einem Rasterschritt und korrigierten sie mit Stößen
   in Rastergröße; das hat die Drift mit Stoffbedingungen noch verstärkt.

Das erklärt alle drei Beobachtungen: y ist exakt, weil dort Weg und
Geschwindigkeit null sind und nichts zu runden ist. `f64` hat denselben
Fehler, nur zehn Größenordnungen kleiner. Und der Fehler wächst mit der
Zeit, weil das Tuch sich vom Ursprung entfernt und das Raster dort gröber
wird.

Behebung in `cloth/step.rs`: Die Lage wird als fehlerfreie Summe
weitergerückt und der Rest je Teilchen gemerkt; die Verbindungen rechnen mit
Lage plus Rest; die Geschwindigkeit wird nur noch um die Korrekturen des
Teilschritts geändert, nicht mehr aus dem Weg zurückgerechnet.

Freier Flug in `f32` nach 60 s, alle Korrekturen eingeschaltet:

| | vorher | nachher |
|---|---|---|
| Mittlerer Geschwindigkeitsfehler | bis 4·10⁻³ m/s | exakt 0 |
| Streuung der Geschwindigkeit zwischen Teilchen | 1,2·10⁻³ m/s | exakt 0 |
| Ortsfehler des Schwerpunkts | bis 10 cm | 10⁻⁷ m |

Regressionstest: `cloth::tests::free_flight_does_not_drift`. Der Test
`zero_iterations_preserves_free_particle_path` bildete die alte Rechnung
bitgenau nach; er prüft jetzt, dass freie Teilchen ihre Geschwindigkeit
bitgenau behalten.

In den Beispielen stand an drei Stellen `std::array` statt `core::array`
(vom Workspace-Lint verboten); das ist mit korrigiert.

## Fehlersuche mit gezielten Versuchen

Fünf Versuche am Gewebe und einer an der Strömung, alle als Rechnung mit
Zahlen statt im Fenster. Zwei Fehler gefunden und behoben, eine Lücke
geschlossen, drei Eigenschaften des Verfahrens vermessen.

### Behoben: ein ruhig hängendes Tuch schaukelte sich auf

Ein Tuch aus 24 × 24 Teilchen, senkrecht in seiner Ruhelage aufgehängt, mit
den Standardwerten: nach einer halben Sekunde mehrere Meter je Sekunde, dann
in der eigenen Ebene zusammengefaltet (Dehnung über 60 %) oder bis zum
Überlauf davongelaufen. Mit 12 × 12 trat es nicht auf, deshalb fiel es
bisher nicht auf. Der Stand vor der Drift-Korrektur zeigte dasselbe.

Ursache: Die Verbindungen wurden stets in derselben Reihenfolge korrigiert.
Das ist als Ganzes unsymmetrisch und kann Fehler im Kreis drehen, statt sie
nur zu verkleinern. Jetzt laufen die Durchläufe paarweise vor und zurück;
eine ungerade Zahl wird aufgerundet.

| Hängendes Tuch, 10 s | vorher | nachher |
|---|---|---|
| 24 × 24, Standardwerte: Dehnung am Ende | 70 % | 1,6 % |
| 24 × 24: Höchsttempo | über 5 m/s | 0,14 m/s |
| 40 × 40: Resttempo | 6 m/s | 0 |
| 48 × 48, 20 Teilschritte, 1 Durchlauf | Überlauf | Dehnung 0,7 %, Ruhe |

### Behoben: nichts haftete

Die Reibung kürzte das Gleiten um einen festen Anteil je Teilschritt. Ein
mittig auf die Kugel gefallenes Tuch kroch deshalb langsam herunter (nach
15 s lag der Schwerpunkt 30 bis 50 cm neben der Mitte). Jetzt gilt Coulomb:
haften unter der Grenzsteigung, sonst rutschen mit `g·(sin α − µ·cos α)`.
Auf der schiefen Ebene stimmt der Weg auf 0,2 % mit der Formel überein; das
Tuch bleibt auf der Kugel liegen (Schwerpunkt 2 bis 5 mm neben der Mitte).

`ClothParams::friction` ist damit ein Reibwert, kein Anteil mehr.

### Lücke geschlossen: Reibung an Dreiecken und zwischen Tuchlagen

Dort gab es gar keine. Ein Tuch rutschte von jeder geneigten Platte aus
Dreiecken reibungsfrei ab. Jetzt gilt dieselbe Coulomb-Reibung.

### Vermessen, kein Fehler

- **Drehender freier Flug:** Impuls exakt erhalten, die Sicherung greift
  nicht ein. Der Drehimpuls sinkt bei 3 rad/s um 1,8 % je Sekunde,
  proportional zur Schrittweite – die bekannte Dämpfung des Verfahrens.
- **Liegen auf einer Platte aus Dreiecken, 20 s:** ruhig, genau im Abstand
  der Tuchdicke, kein Eingriff der Sicherung, keine Durchdringung.
- **Abhängigkeit von der Schrittweite:** Die Dehnung unter Last halbiert
  sich etwa mit doppelt so vielen Teilschritten.
- **Kanalströmung (Gitter):** Durchfluss exakt, Parabelprofil auf 0,8 %
  genau (Test `fluid::tests::poiseuille`).

Einzelheiten und Zahlen: `todo/azphysics.md`, Abschnitte „Durchlaufrichtung“,
„Reibung“ und „Was das Verfahren nicht kann“.
