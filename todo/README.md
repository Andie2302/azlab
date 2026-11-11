# azlab – Planung

Dieser Ordner ist die verbindliche Planung des Workspace. **Code entsteht erst,
wenn der zugehörige Abschnitt hier geplant und freigegeben ist.**

## Regeln (aus `crates/azmath/src/lib.rs`, gelten für den ganzen Workspace)

1. Alles, was nicht hinter Feature-Flags steht, MUSS in `no_std` ohne Features laufen.
2. Neuer Rust-Stil: nirgends `mod.rs`; in `lib.rs` steht kein Code.
3. Erst planen (Abhängigkeiten, Feature-Flags, Ziele, Nicht-Ziele), dann Code. Gilt besonders für KI.
4. KI führt von sich aus keinen Code aus – nicht zum Testen, nicht zum Rechnen, nicht zum Prüfen.

Folgerungen für die Zusammenarbeit:

- Prüfbefehle (`cargo test`, `cargo clippy`, ...) führt der Mensch aus; die KI
  nennt nur, was zu prüfen ist (siehe [05-qualitaet.md](05-qualitaet.md)).
- Zahlenwerte in Tests (Referenzwerte) stammen aus zitierten Quellen, nicht
  aus KI-Rechnungen.

## Dokumente

| Datei | Inhalt | Status |
|-------|--------|--------|
| [01-ziele.md](01-ziele.md) | Ziele und Nicht-Ziele | Entwurf |
| [02-architektur.md](02-architektur.md) | Schichten, Traits, Backend, Erweiterbarkeit | Entwurf |
| [03-features-abhaengigkeiten.md](03-features-abhaengigkeiten.md) | Feature-Flags, Abhängigkeiten, Toolchain | Entwurf |
| [04-roadmap.md](04-roadmap.md) | Phasen mit Aufgabenlisten | Entwurf |
| [05-qualitaet.md](05-qualitaet.md) | Tests, Prüfmatrix, Doku | Entwurf |
| [06-offene-fragen.md](06-offene-fragen.md) | Entscheidungen, die noch fehlen | offen |

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
