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

4. **KI prüft, testet und führt absolut NICHTS aus, ohne vorher zu fragen.**
   Ohne ausdrückliche Erlaubnis des Menschen für genau diesen Fall ist
   verboten:
   - Code ausführen: `cargo build`, `test`, `clippy`, `run`, `doc`, Skripte,
     Testschleifen, Hintergrundprozesse.
   - Hilfs-Agenten oder andere Sitzungen starten, die etwas prüfen.
   - Ersatz-Prüfen: Tests oder Abläufe „von Hand“ oder im Kopf nachrechnen,
     fertigen Code noch einmal prüfend durchgehen.
   - Jede andere Form von Prüfen, Testen, Nachrechnen oder Verifizieren.

   Dazu gilt:
   - Fragen heißt: vorher sagen, **was** geprüft werden soll, **wozu** und
     **wie lange** es voraussichtlich dauert – und die Antwort abwarten.
   - Eine Erlaubnis gilt nur für das, wofür sie gegeben wurde, nicht für
     spätere Schritte oder Phasen.
   - Normalfall: Die KI schreibt und pusht. Der Mensch führt
     `scripts/pruefen.sh` lokal aus und pusht `pruefung/ergebnis.txt`. Die KI
     liest das Ergebnis und behebt die gemeldeten Fehler.
   - Dateien lesen, schreiben und Git-Befehle (commit, pull, push) sind
     erlaubt.
   - Hintergrund: Testschleifen, parallele Prüf-Agenten und stundenlang
     laufende Skripte haben schon ganze Wochen-Kontingente verbraucht.

5. **Bei Fragen sofort nachfragen**, auch wenn es mehrere sind.

## Weitere Festlegungen

- Nur stable Rust; die Mindestversion darf bis zur neuesten stabilen Version gehen.
- Gearbeitet wird direkt auf `master`.
- Plan, Doku-Kommentare und Commit-Nachrichten auf Deutsch; Bezeichner englisch.
- Neue Abhängigkeiten und Features zuerst in
  `todo/03-features-abhaengigkeiten.md` eintragen.
- Zahlenwerte in Tests stammen aus zitierten Quellen oder sind von Hand exakt
  herleitbar, nicht von der KI ausgerechnet.
