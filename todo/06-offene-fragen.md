# Offene Fragen

Entschiedene Fragen stehen in der Tabelle in [README.md](README.md).

## F12 – Division bei nicht-kommutativen Typen

Für Quaternionen ist `a / b` nicht eindeutig. Vorschlag: der Operator `/`
bedeutet immer Rechtsdivision `a * b⁻¹`; für die andere Reihenfolge schreibt
man `b.recip() * a` aus.

Entscheidung: _offen_ (im Detailplan Phase 2 so angenommen)

## F13 – Zweige der komplexen Funktionen

Vorschlag: Hauptzweige wie in C99/IEEE (Schnitt von `ln` und `sqrt` entlang
der negativen reellen Achse, `arg` in (−π, π]).

Entscheidung: _offen_ (im Detailplan Phase 2 so angenommen)

## F14 – Intervallarithmetik und Rust-Version

Beweisbar konservative Kollisionstests brauchen Rundung nach außen. Auf
stable geht das über `next_up`/`next_down` (ab Rust 1.86). Dafür müsste
`rust-version` von 1.85 auf 1.86 steigen. In Ordnung?

Entscheidung: _offen_ (betrifft erst Phase 4)
