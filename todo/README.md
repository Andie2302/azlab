# azlab – Planung

Dieser Ordner ist die verbindliche Planung des Workspace. **Code entsteht erst,
wenn der zugehörige Abschnitt hier geplant und freigegeben ist.**

## Regeln

Die Projektregeln stehen ausschließlich in [`CLAUDE.md`](../CLAUDE.md) im
Wurzelverzeichnis.

## Dokumente

| Datei | Inhalt | Status |
|-------|--------|--------|
| [01-ziele.md](01-ziele.md) | Ziele und Nicht-Ziele | Entwurf |
| [02-architektur.md](02-architektur.md) | Schichten, Traits, Backend, Erweiterbarkeit | Entwurf |
| [03-features-abhaengigkeiten.md](03-features-abhaengigkeiten.md) | Feature-Flags, Abhängigkeiten, Toolchain | Entwurf |
| [04-roadmap.md](04-roadmap.md) | Phasen mit Aufgabenlisten | Entwurf |
| [05-qualitaet.md](05-qualitaet.md) | Tests, Prüfmatrix, Doku | Entwurf |
| [06-offene-fragen.md](06-offene-fragen.md) | Entscheidungen, die noch fehlen | offen |
| [phase-2-zahlentypen.md](phase-2-zahlentypen.md) | Detailplan Phase 2 | umgesetzt und geprüft |
| [phase-3-lineare-algebra.md](phase-3-lineare-algebra.md) | Detailplan Phase 3 | umgesetzt und geprüft |
| [phase-4-nullstellen-polynome.md](phase-4-nullstellen-polynome.md) | Detailplan Phase 4 | umgesetzt und geprüft |
| [phase-5-lie-theorie.md](phase-5-lie-theorie.md) | Detailplan Phase 5 | umgesetzt, noch ungeprüft |
| [phase-6a-sichten-iterative-loeser.md](phase-6a-sichten-iterative-loeser.md) | Detailplan Phase 6, Teil 1 | umgesetzt, noch ungeprüft |

Status-Stufen: `Entwurf` → `freigegeben` → `umgesetzt`.

## Ablauf pro Phase

1. Phase in `04-roadmap.md` verfeinern: öffentliche API (Typen, Traits,
   Signaturen), Dateien, Tests, Referenzquellen.
2. Offene Fragen der Phase in `06-offene-fragen.md` klären.
3. Freigabe durch den Menschen.
4. Umsetzung, danach Prüfmatrix durch den Menschen.
5. Roadmap abhaken, Status aktualisieren.

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
| Kollisionsgeometrie | nicht in azmath, sondern im späteren Crate `azphysics` |
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
