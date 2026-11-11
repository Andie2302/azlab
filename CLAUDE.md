# azlab – Projektregeln

Einziger Ort für die Regeln des Workspace. Sie haben Vorrang vor Gewohnheiten
der KI. Zuerst diese Datei lesen, dann [todo/README.md](todo/README.md).

## Regeln

1. **`no_std` zuerst**: Alles, was nicht hinter Feature-Flags steht, MUSS in
   `no_std` ohne Features laufen. Was Heap oder `std` braucht, darf es
   benutzen, steht dann aber hinter einem Feature-Flag.

2. **Neuer Rust-Stil**: nirgends `mod.rs`. In `lib.rs` steht kein Code –
   erlaubt sind Doku, Deklarationen (`mod`, `extern crate`), `pub use` und
   crate-weite Attribute bzw. Feature-Schalter (z. B. `alloc` freischalten).
   Dasselbe gilt für Modul-Dateien, neben denen ein gleichnamiger Ordner liegt.

3. **Erst planen, dann Code**: Bevor Code entsteht, muss das Vorhaben im
   Ordner `todo/` geplant und freigegeben sein – inklusive Abhängigkeiten,
   Feature-Flags, Zielen und Nicht-Zielen. Gilt ganz besonders für KI.

4. **KI führt selbstständig keinen Code aus**: nicht zum Testen, nicht zum
   Rechnen, nicht zum schnellen Prüfen. Das betrifft `cargo build`, `test`,
   `clippy`, `run`, Skripte, Testschleifen, Hintergrundprozesse und
   Hilfs-Agenten, die Änderungen gegenprüfen.
   - Hintergrund: Testschleifen, parallele Prüf-Agenten und stundenlang
     laufende Skripte haben schon ganze Wochen-Kontingente verbraucht.
   - Erlaubt ist, **vorher zu fragen** – dann mit Angabe, was ausgeführt
     werden soll, wozu und wie lange es voraussichtlich dauert.
   - Normalfall: Der Mensch führt `scripts/pruefen.sh` lokal aus und pusht
     `pruefung/ergebnis.txt`. Die KI pflegt das Skript und liest das Ergebnis.
   - Dateien lesen, schreiben und Git-Befehle (commit, push) sind kein
     „Code ausführen“.

5. **Bei Fragen sofort nachfragen**, auch wenn es mehrere sind.

## Weitere Festlegungen

- Nur stable Rust.
- Plan, Doku-Kommentare und Commit-Nachrichten auf Deutsch; Bezeichner englisch.
- Neue Abhängigkeiten und Features zuerst in
  `todo/03-features-abhaengigkeiten.md` eintragen.
- Zahlenwerte in Tests stammen aus zitierten Quellen oder sind von Hand exakt
  herleitbar, nicht von der KI ausgerechnet.
