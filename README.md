# azlab

**Ein Rust-Workspace für mathematische Grundlagen und darauf aufbauende Fachbibliotheken.**

Das Fundament ist **azmath**: eine generische Mathematik-Bibliothek für Zahlentypen,
lineare Algebra, automatische Differentiation und Gleichungslöser. Ihr Kern
arbeitet mit `no_std` ohne Heap-Allokationen und eignet sich damit auch für
Embedded-Anwendungen.

Fachliche Modelle – beispielsweise Stoffdaten, Physik oder Robotik – sollen
später in eigenen Crates darauf aufbauen. azmath selbst enthält keine Einheiten,
Stoffkonstanten oder Fachmodelle.

> **Entwicklungsstand:** Das Projekt befindet sich im Aufbau. Derzeit enthält
> der Workspace das Crate `azmath`. Die API kann sich noch ändern; eine
> Veröffentlichung auf crates.io ist aktuell nicht vorgesehen
> (`publish = false`).

## Grundsätze

- **`no_std` zuerst:** Der Kern benötigt weder Betriebssystem noch Heap.
- **Generische Algorithmen:** Skalar-Traits beschreiben die benötigten
  Rechenfähigkeiten. Neben `f32` und `f64` können passende Algorithmen auch mit
  komplexen Zahlen oder Dualzahlen arbeiten.
- **Feste Größen und eigene Puffer:** Vektoren und Matrizen verwenden const
  generics; etwa für Polynom-Nullstellen stellt der Aufrufer Ausgabepuffer bereit.
- **Austauschbare Rechen-Backends:** Elementarfunktionen werden über `libm`
  oder `std` bereitgestellt, ohne die öffentliche API des Algorithmus zu ändern.
- **Kein `unsafe`:** Im Workspace durch einen Rust-Lint verboten.
- **Explizite Fehler:** Fallible Löser liefern `Result`; iterative Verfahren
  geben Näherung, Restfehler und Iterationszahl zurück. Bei Nicht-Konvergenz
  bleibt der letzte Stand im Fehlerwert verfügbar.

Die Indexoperatoren von Vektoren und Matrizen folgen der Array-Konvention und
können bei ungültigen Indizes eine Panic auslösen. Für geprüften Zugriff stehen
`get` und `get_mut` bereit.

## Bereits implementiert

| Bereich | Umfang |
|---------|--------|
| Skalar-Traits | `Algebra`, `Ring`, `Field`, `Commutative`, `Conjugate`, `Real`, `Elementary`, `RealElementary` |
| Zahlentypen | `Complex<T>`, `Dual<T>`, `DualN<T, N>`, `CayleyDickson<A>`, `Octonion`, `Quaternion<T>` |
| Rotationen | Quaternion-Rotationen, Interpolation mit Slerp und Umwandlungen |
| Lineare Algebra | `Vector<T, N>`, `Matrix<T, R, C>`, Produkte, Normen, Transponieren und Adjungieren |
| Matrixzerlegungen | LU mit Pivotisierung, Cholesky, QR, lineare Systeme, Determinante und Inverse |
| Autodiff | Funktionswert und Ableitungen, Gradient, Jacobi- und Hesse-Matrix über Dualzahlen |
| Nullstellen | Bisektion, Sekante, Newton, Newton mit Klammer, Brent |
| Nichtlineare Systeme | Newton mit vorgegebener oder automatisch berechneter Jacobi-Matrix, Dämpfung und Liniensuche |
| Polynome | Horner- und Tschebyschow-Auswertung, Formeln für Grad 2–4, komplexe Nullstellen mit Aberth-Ehrlich |
| Intervalle | `Interval<T>` mit Rundung nach außen und reelle Polynom-Nullstellenisolation |

Welche Zahlentypen ein Verfahren unterstützt, bestimmen seine Trait-Grenzen.
Nicht jede Kombination aus Algorithmus und Zahlentyp ist sinnvoll oder verfügbar.

## Einbinden

Vorausgesetzt wird stable Rust. Das Workspace-Manifest deklariert Rust
`1.86` als Mindestversion und verwendet die Edition `2024`.

Da das Projekt noch nicht auf crates.io veröffentlicht wird, kann azmath als
Git-Abhängigkeit eingebunden werden:

```toml
[dependencies]
azmath = { git = "https://github.com/Andie2302/azlab.git", branch = "master" }
```

Für reproduzierbare Abhängigkeiten kann statt `branch` ein konkreter Commit
über `rev` festgelegt werden.

Bei einem lokalen Checkout ist alternativ eine Pfad-Abhängigkeit möglich.
Den Pfad an die eigene Verzeichnisstruktur anpassen:

```toml
[dependencies]
azmath = { path = "../azlab/crates/azmath" }
```

### Features und Backends

| Konfiguration | Wirkung |
|---------------|---------|
| Standard: `libm` | `no_std` ohne Heap; Elementarfunktionen über die optionale Abhängigkeit `libm` |
| `default-features = false` | Keine Abhängigkeiten; keine Implementierungen von `Elementary` und `RealElementary` für `f32`/`f64` |
| `std` | Elementarfunktionen über `std`; aktiviert zusätzlich `alloc` |
| `alloc` | Vorbereiteter Schalter für spätere Heap-basierte Typen; dynamische Matrizen und Tensoren sind noch geplant |

Sind `std` und `libm` gleichzeitig aktiv, verwendet azmath das
`std`-Backend. Die Features bleiben additiv.

Ohne Rechen-Backend:

```toml
[dependencies]
azmath = { git = "https://github.com/Andie2302/azlab.git", default-features = false }
```

Mit `std`-Backend:

```toml
[dependencies]
azmath = { git = "https://github.com/Andie2302/azlab.git", default-features = false, features = ["std"] }
```

## Beispiele

Die folgenden Beispiele zeigen die Verwendung in einer normalen Rust-Anwendung.
`println!` gehört zur Anwendung; der mathematische Kern benötigt dafür kein
`std`. Die Beispiele wurden für diese README nicht ausgeführt.

### Ein lineares Gleichungssystem lösen

Die Matrix und die rechte Seite werden mit fester Größe angelegt.
Fehler wie eine singuläre Matrix werden über `Result` weitergegeben.

```rust
use azmath::{LinalgError, Matrix, Vector};

fn main() -> Result<(), LinalgError> {
    let a = Matrix::new([[2.0_f64, 1.0], [1.0, 3.0]]);
    let b = Vector::new([3.0_f64, 5.0]);

    let x = a.solve(b)?;
    println!("Lösung: {x:?}");
    Ok(())
}
```

Für mehrere rechte Seiten kann die Matrix einmal mit `a.lu()?` zerlegt und
die Zerlegung anschließend über `Lu::solve` wiederverwendet werden.

### Eine generische Funktion automatisch ableiten

Die Funktion wird über einen Skalar-Trait formuliert. Dadurch lässt sie sich
sowohl mit einem Gleitkommatyp als auch mit Dualzahlen aufrufen.

```rust
use azmath::autodiff::derivative;
use azmath::scalar::Elementary;

fn model<T: Elementary>(x: T) -> T {
    x.sin() * x.exp()
}

fn main() {
    let (value, slope) = derivative(model, 1.0_f64);
    println!("Funktionswert: {value}, Ableitung: {slope}");
}
```

Dieses Beispiel benötigt ein Backend für die Elementarfunktionen, etwa das
standardmäßig aktivierte `libm`. Für mehrere Veränderliche stehen
`gradient`, `jacobian` und `hessian` zur Verfügung.

### Eine Nullstelle mit Brent suchen

Brent nutzt ein Intervall mit Vorzeichenwechsel. Abbruchkriterien werden über
`Tolerance` vorgegeben; das Ergebnis enthält auch Restfehler und Iterationszahl.

```rust
use azmath::solve::brent;
use azmath::{SolveError, Tolerance};

fn main() -> Result<(), SolveError<f64, f64>> {
    let tolerance = Tolerance::<f64>::default();
    let solution = brent(
        |x: f64| x * x - 2.0,
        1.0,
        2.0,
        &tolerance,
    )?;

    println!(
        "Näherung: {}, Restfehler: {}, Iterationen: {}",
        solution.x, solution.residual, solution.iterations,
    );
    Ok(())
}
```

Brent und die übrigen skalaren Nullstellenverfahren benötigen selbst kein
Rechen-Backend. Eine übergebene Funktion kann je nach verwendeten Operationen
dennoch eines benötigen.

## Aufbau des Repositorys

| Pfad | Inhalt |
|------|--------|
| [crates/azmath](crates/azmath) | Mathematische Kernbibliothek |
| [todo](todo/README.md) | Ziele, Architektur, Entscheidungen und Detailpläne |
| [scripts/pruefen.sh](scripts/pruefen.sh) | Lokale Prüfmatrix |
| [pruefung/ergebnis.txt](pruefung/ergebnis.txt) | Eingechecktes Prüfprotokoll |
| [CLAUDE.md](CLAUDE.md) | Verbindliche Projektregeln und Entwicklungsablauf |

## Dokumentation und lokale Prüfung

Die API-Dokumentation lässt sich lokal erzeugen:

```bash
cargo doc -p azmath --no-deps --all-features --open
```

Für die vollständige Prüfmatrix führt der Entwickler im Repository aus:

```bash
scripts/pruefen.sh
```

Das Skript prüft Tests und Clippy für die Kombinationen ohne Features, mit
`libm`, mit `std`, mit `alloc` und mit allen Features. Zusätzlich erzeugt es
die API-Dokumentation und baut, sofern das Target installiert ist, für
`thumbv7em-none-eabihf` mit und ohne Standard-Features.

Das Embedded-Target kann zuvor installiert werden:

```bash
rustup target add thumbv7em-none-eabihf
```

Das Ergebnis wird nach `pruefung/ergebnis.txt` geschrieben.
Prüfläufe erfolgen gemäß den Projektregeln durch den Menschen; die Befehle
hier dienen als Anleitung.

## Geplanter Ausbau

Die folgenden Bereiche sind geplant und noch nicht als fertige Funktionen
verfügbar:

- Lie-Gruppen und Lie-Algebren für geometrische und robotische Anwendungen.
- Matrix- und Tensor-Sichten auf Puffer, dünnbesetzte und iterative Löser.
- Dynamische Matrizen und Tensoren hinter dem Feature `alloc`.
- ODE-Verfahren, Quadratur und Interpolation.
- PDE-Verfahren auf Gittern und Finite Elemente.
- Rückwärts-Autodiff und Optimierung.

Mögliche spätere Fach-Crates sollen darauf aufbauen. Dazu gehören beispielsweise
Physik-Simulationen und Stoffdatenmodelle; sie sind derzeit kein Bestandteil von
azmath.

Details stehen in der [Roadmap](todo/04-roadmap.md), die Abgrenzung in den
[Zielen und Nicht-Zielen](todo/01-ziele.md).

## Entwicklung

Erweiterungen werden vor der Implementierung in `todo/` geplant und freigegeben.
Neue Abhängigkeiten und Features werden dort ebenfalls dokumentiert.
Dokumentation und Commit-Nachrichten sind deutsch, Bezeichner englisch.

Die vollständigen Regeln stehen in [CLAUDE.md](CLAUDE.md).

## Lizenz

Im Workspace-Manifest ist `MIT OR Apache-2.0` als Lizenz angegeben.
