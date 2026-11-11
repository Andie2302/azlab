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

4. **ABSOLUTES VERBOT JEDER PRÜFUNG, BERECHNUNG ODER SIMULATION DURCH DIE KI**
   Die KI darf unter keinen Umständen aktiv werden, um Code, Logik oder Zahlen zu verifizieren. 
   Ohne ausdrückliche, vorherige Einzelerlaubnis des Menschen ist **strikt verboten**:
   - **Ausführen von Tools/Skripten**: `cargo build`, `test`, `clippy`, `run`, `doc`, Testschleifen, Hintergrundprozesse.
   - **Agenten & Hilfssitzungen**: Starten anderer Instanzen oder Agenten zur Prüfung.
   - **Mentales oder schriftliches Nachrechnen (STRICT NO-CALCULATION-RULE)**: 
     - Es ist der KI **absolut untersagt**, mathematische Ausdrücke, Testergebnisse, Matrix-Operationen, Formeln oder numerische Werte „per Hand“, im Kopf oder im Ausgabetext auszurechnen, zu simulieren oder zu verifizieren.
     - Wenn Testdaten, Erwartungswerte oder Berechnungen benötigt werden, darf die KI **keinerlei eigene Werte generieren**. Sie muss Platzhalter (z. B. `/* HIER_WERT_AUS_QUEELLE_EINTRAGEN */`) verwenden oder den Menschen explizit darum bitten, die exakten Zahlen bereitzustellen.
   - **Prüfendes Durchgehen**: Fertigen Code gedanklich „durchzuspielen“, um zu schauen, ob er fehlerfrei ist.
   - Die Erlaubnis gilt ausschließlich für den exakten, einmaligen Zweck und darf niemals übertragen werden.

   **Dazu gilt:**
   - Fragen heißt: vorher sagen, **was** geprüft werden soll, **wozu** und **wie lange** es voraussichtlich dauert – und die Antwort abwarten.
   - Normalfall: Die KI schreibt und pusht. Der Mensch führt `scripts/pruefen.sh` lokal aus und pusht `pruefung/ergebnis.txt`. Die KI liest das Ergebnis und behebt die gemeldeten Fehler.
   - Dateien lesen, schreiben und Git-Befehle (`commit`, `pull`, `push`) sind erlaubt.
   - Hintergrund: Testschleifen, parallele Prüf-Agenten und stundenlang laufende Skripte haben schon ganze Wochen-Kontingente verbraucht.

5. **Bei Fragen sofort nachfragen**, auch wenn es mehrere sind.

## Weitere Festlegungen

- Nur stable Rust; die Mindestversion darf bis zur neuesten stabilen Version gehen.
- Gearbeitet wird direkt auf `master`.
- Plan, Doku-Kommentare und Commit-Nachrichten auf Deutsch; Bezeichner englisch.
- Neue Abhängigkeiten und Features zuerst in `todo/03-features-abhaengigkeiten.md` eintragen.
- Zahlenwerte in Tests stammen aus zitierten Quellen oder sind von Hand exakt herleitbar, **niemals** von der KI ausgerechnet oder erfunden. Falls Werte fehlen, werden Platzhalter genutzt.