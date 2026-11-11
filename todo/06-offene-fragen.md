# Offene Fragen

Jede Frage hat eine Empfehlung. Entscheidung hier eintragen, dann wandert
sie in die Tabelle in `README.md`.

## F1 – Standard-Features von azmath

Soll `default = ["libm"]` gelten oder `default = []`?

- `[]`: entspricht Regel 1 am wörtlichsten; wer azmath einbindet, muss ein
  Backend wählen, sonst haben `f32`/`f64` keine Elementarfunktionen.
- `["libm"]`: funktioniert sofort, bleibt `no_std` ohne `alloc`, bringt aber
  standardmäßig eine Abhängigkeit mit.

Empfehlung: `default = ["libm"]`; Crates im Workspace binden ohnehin mit
`default-features = false` ein.

Entscheidung: _offen_

## F2 – Was genau darf in `lib.rs` stehen?

„Kein Code in `lib.rs`“ – Annahme: erlaubt sind Crate-Doku, `#![no_std]`,
`mod`-Deklarationen und `pub use`; verboten sind Funktionen, Typen, Traits,
Makros. Gilt dasselbe für Modul-Dateien wie `scalar.rs`, neben denen ein
Ordner `scalar/` liegt?

Empfehlung: ja, beide nur Deklarationen und Re-Exporte.

Entscheidung: _offen_

## F3 – „Stoffsimulation“

Gemeint ist Textil/Gewebe (Cloth) oder Stoff im Sinn von Materie/Stoffdaten?
Ändert die Prioritäten in Phase 6/7 (Cloth braucht vor allem implizite
Zeitintegration und dünnbesetzte Löser).

Entscheidung: _offen_

## F4 – Ein Crate oder mehrere?

azmath als ein Crate mit Schichten als Modulen, oder von Anfang an getrennt
(`azmath-core`, `azmath-linalg`, ...)?

Empfehlung: ein Crate. Die Schichtenregel (Abhängigkeiten nur nach unten)
hält eine spätere Aufteilung offen, ohne jetzt Verwaltungsaufwand zu erzeugen.

Entscheidung: _offen_

## F5 – Quaternion: eigener Typ oder Alias auf Cayley-Dickson?

- Alias `CayleyDickson<Complex<T>>`: mathematisch elegant, aber unhandliche
  Felder und Rotations-API.
- Eigener Typ `Quaternion<T>` mit `w, x, y, z` plus Umwandlung.

Empfehlung: eigener Typ, Cayley-Dickson daneben als generische Konstruktion.

Entscheidung: _offen_

## F6 – Prüfmatrix: Skript, CI oder beides?

Da die KI nichts ausführt, braucht es einen bequemen Weg für den Menschen:
ein Shell-Skript im Repo, eine GitHub-Actions-CI oder beides.

Empfehlung: beides; CI führt dasselbe Skript aus.

Entscheidung: _offen_

## F7 – Ganzzahlen in den Traits

Sollen `i8`…`i128` `Ring` implementieren (Ganzzahl-Matrizen, exakte
Rechnung), und brauchen vorzeichenlose Typen eine eigene Stufe ohne Negation
(`Semiring`)?

Empfehlung: vorzeichenbehaftete als `Ring` in Phase 1, `Semiring` erst bei
Bedarf.

Entscheidung: _offen_

## F8 – Verbleib der Regeln in `lib.rs`

Die vier Regeln stehen als Kommentar in `crates/azmath/src/lib.rs`. Sollen
sie dort bleiben, oder nur noch in `todo/README.md` (und ggf. einer
`CLAUDE.md`, damit KI-Sitzungen sie sicher zuerst lesen)?

Empfehlung: `CLAUDE.md` im Repo-Wurzelverzeichnis mit Verweis auf
`todo/README.md`.

Entscheidung: _offen_
