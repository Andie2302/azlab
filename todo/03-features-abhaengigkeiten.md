# Feature-Flags, Abhängigkeiten, Toolchain

## Feature-Flags von azmath

| Feature | Aktiviert | Abhängigkeit | Bemerkung |
|---------|-----------|--------------|-----------|
| *(keins)* | Skalar-Traits, Zahlentypen, lineare Algebra fester Größe, Löser – alles generisch | keine | `f32`/`f64` haben dann kein `Elementary` |
| `libm` | Backend für Elementarfunktionen in `no_std` | `libm` (optional) | |
| `std` | Backend über die Standardbibliothek; schaltet `alloc` mit ein | – | hat Vorrang vor `libm` |
| `alloc` | besitzende dynamische Typen, wachsendes Autodiff-Tape | – | |

Geplant, noch nicht festgelegt: weitere Backend-Features; `serde` (optional)
für Zahlentypen; `nightly` gibt es nicht.

Regeln:

- Features sind additiv und dürfen in beliebiger Kombination bauen.
- Rangfolge der Backends bei mehreren aktiven: `std` > `libm` > weitere.
- Aufbauende Crates binden azmath mit `default-features = false` ein und
  reichen `std`/`libm`/`alloc` als eigene Features durch.
- Standard-Features: `default = ["libm"]` (entschieden).

## Abhängigkeiten

| Crate | Wofür | Art |
|-------|-------|-----|
| `libm` 0.2 | Elementarfunktionen ohne `std` | optional, nur mit Feature `libm` |

Sonst keine. Jede neue Abhängigkeit wird zuerst hier eingetragen, ist
optional, `no_std`-fähig und wird zentral unter `[workspace.dependencies]`
versioniert. Keine Abhängigkeit auf `num-traits`, `nalgebra` o. Ä.: die
Trait-Hierarchie gehört azmath selbst, damit sie frei erweiterbar bleibt.

## Toolchain

- Nur stable Rust, Edition 2024. `rust-version` im Workspace steht derzeit
  auf 1.85 und darf bei Bedarf bis zur neuesten stabilen Version angehoben
  werden (Stand: 1.99.1).
- `f16`/`f128`: sobald stabil, über dieselben Makros wie `f32`/`f64`
  anbinden. Dafür ist ein Backend nötig, das diese Typen kann.
- Folgen von „nur stable“: keine Größenarithmetik in const generics, keine
  Spezialisierung. Die Architektur ist darauf ausgelegt.

## Workspace-Einstellungen (bestehend, bleiben)

- `unsafe_code = "forbid"`, `missing_docs = "warn"`.
- Clippy: `std_instead_of_core`, `std_instead_of_alloc`,
  `alloc_instead_of_core`, `unwrap_used`, `panic` verboten.
- Zu ergänzen bei Umsetzung von Phase 0: `[lints] workspace = true` in
  `crates/azmath/Cargo.toml` (sonst greifen die Lints nicht),
  `[workspace.dependencies]`, gemeinsame `version`/`rust-version`.
