# azlab – Planung

Dieser Ordner ist die verbindliche Planung des Workspace. **Code entsteht erst,
wenn der zugehörige Abschnitt hier geplant und freigegeben ist.**

## Arbeitsweise

Die früheren Projektregeln (`CLAUDE.md`) sind vom Menschen aufgehoben: Die KI
arbeitet selbstständig, darf ausführen, testen und bestehenden Code ändern. Als Entwurfsgrundsätze bleiben:

- `no_std` ohne `alloc` ist der Normalfall; alles andere steht hinter
  Feature-Flags.
- Kein `mod.rs`; `lib.rs` und Modul-Dateien mit Unterordner enthalten nur
  Doku, Deklarationen und Re-Exporte.
- Nur stable Rust. Doku und Commits deutsch, Bezeichner englisch.
- Jeder Stand auf `master` besteht `scripts/pruefen.sh`.

## Dokumente

| Datei | Inhalt | Status |
|-------|--------|--------|
| [01-ziele.md](01-ziele.md) | Ziele und Nicht-Ziele | gültig |
| [02-architektur.md](02-architektur.md) | Schichten, Traits, Backend, Erweiterbarkeit | gültig; Dateistruktur siehe Quelltext |
| [03-features-abhaengigkeiten.md](03-features-abhaengigkeiten.md) | Feature-Flags, Abhängigkeiten, Toolchain | gültig |
| [04-roadmap.md](04-roadmap.md) | Phasen mit Aufgabenlisten | alle Phasen umgesetzt |
| [05-qualitaet.md](05-qualitaet.md) | Tests, Prüfmatrix, Doku | gültig |
| [06-offene-fragen.md](06-offene-fragen.md) | Entscheidungen, die noch fehlen | keine offen |
| `phase-*.md` | Detailpläne der Phasen 2 bis 7 | umgesetzt und geprüft |
| [azphysics.md](azphysics.md) | Aufbau, Zusagen und Grenzen von azphysics und azviz | erste Fassung umgesetzt |

Die Phasen 8 bis 11 und der zweite Teil von Phase 7 wurden ohne eigenen
Detailplan umgesetzt; ihr Umfang steht in der Roadmap und in der Doku der
Module.

## Bereits getroffene Entscheidungen

| Thema | Entscheidung |
|-------|--------------|
| Toolchain | nur stable; `f16`/`f128` erst nach Stabilisierung |
| Elementarfunktionen | austauschbares Backend (`libm`, `std`, später weitere); weder azmath-Nutzer noch aufbauende Crates sehen, was rechnet |
| „He-Klammern/He-Algebra“ | gemeint sind Lie-Klammern und Lie-Algebren |
| Sprache | Plan, Doku-Kommentare und Commits deutsch; Bezeichner englisch |
| Standard-Features | `default = ["libm"]` |
| Inhalt von `lib.rs` | Doku, Deklarationen, `pub use`, crate-weite Attribute/Feature-Schalter; kein sonstiger Code |
| „Stoffsimulation“ | Gewebe/Kleidung inklusive Kollision |
| Gleichungslöser | Nullstellen, lineare Systeme, Gleichungen höherer Ordnung, ODE und PDE gehören zu azmath |
| PDE-Umfang | Gitterverfahren **und** Finite Elemente |
| Crate-Schnitt | ein Crate `azmath`, Schichten als Module |
| Quaternionen | eigener Typ `Quaternion<T>` plus generische Cayley-Dickson-Konstruktion |
| Ort der Regeln | nur `CLAUDE.md` im Wurzelverzeichnis |
| Prüfmatrix | Shell-Skript `scripts/pruefen.sh`; der Mensch führt es lokal aus und pusht `pruefung/ergebnis.txt` |
| Vorzeichenlose Ganzzahlen | vorerst nicht in den Traits |
| Kollisionsgeometrie | nicht in azmath, sondern im Crate `azphysics` |
| Anzeige der Beispiele | eigenes Crate `azviz`: Fenster über `macroquad`, Dateien über `tiny-skia` und `gif`; nur Entwicklungs-Abhängigkeit von azphysics |
| Aufruf der Beispiele | Kurzbefehl `cargo demo <name>` in `.cargo/config.toml` |
| Heap in Simulationen/FEM | erlaubt, wo die Simulation ihn braucht, aber immer hinter Feature-Flag |
| Gewebe-Verfahren | PBD/XPBD, kontinuierliche Kollisionserkennung (Vertex-Face, Edge-Edge), Abstoßungskräfte |
| Trait-Umbau | Basis-Trait `Algebra` unter `Ring`; `from_f64` in `Field` |
| Dual-Vergleich | nur der Wert zählt, Ableitungsanteil wird ignoriert |
| Komplexe Zahlen | eigener Typ `Complex<T>` plus Umwandlung von/nach Cayley-Dickson |
| Division nicht-kommutativer Typen | `a / b` ist Rechtsdivision `a · b⁻¹` |
| Komplexe Zweige | Hauptzweige wie C99 |
| Rust-Mindestversion | darf bis zur jeweils neuesten stabilen Version gehen (Stand: 1.99.1); bei Bedarf anheben |
| Arbeitszweig | `master` |
| Quaternion-Rotationen | Kern in Phase 2, Matrizen/Euler-Winkel in Phase 3/5 |
