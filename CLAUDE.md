# azlab – Regeln für KI-Sitzungen

Diese Regeln gelten für den gesamten Workspace und haben Vorrang vor
Gewohnheiten der KI. Zuerst lesen, dann [todo/README.md](todo/README.md).

1. **`no_std` zuerst**: Alles, was nicht hinter Feature-Flags steht, MUSS in
   `no_std` ohne Features laufen.
2. **Neuer Rust-Stil**: nirgends `mod.rs`. In `lib.rs` steht kein Code –
   erlaubt sind Doku, Deklarationen (`mod`, `extern crate`), `pub use` und
   crate-weite Attribute bzw. Feature-Schalter (z. B. `alloc` freischalten).
   Dasselbe gilt für Modul-Dateien, neben denen ein gleichnamiger Ordner liegt.
3. **Erst planen, dann Code**: Bevor Code entsteht, muss das Vorhaben im
   Ordner `todo/` geplant sein – inklusive Abhängigkeiten, Feature-Flags,
   Zielen und Nicht-Zielen.
4. **KI führt von sich aus keinen Code aus**: nicht zum Testen, nicht zum
   Rechnen, nicht zum schnellen Prüfen. Kein `cargo build/test/clippy/run`,
   keine Skripte. Prüfbefehle führt der Mensch aus.
5. **Bei Fragen sofort nachfragen**, auch wenn es mehrere sind.

Weitere Festlegungen:

- Nur stable Rust.
- Plan, Doku-Kommentare und Commit-Nachrichten auf Deutsch; Bezeichner englisch.
- Neue Abhängigkeiten und Features zuerst in
  `todo/03-features-abhaengigkeiten.md` eintragen.
