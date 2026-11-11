# Offene Fragen

Entschiedene Fragen stehen in der Tabelle in [README.md](README.md).

## F6 – Prüfmatrix: Skript, CI oder beides?

Da die KI nichts ausführt, braucht es einen bequemen Weg für den Menschen:
ein Shell-Skript im Repo, eine GitHub-Actions-CI oder beides.

Empfehlung: beides; CI führt dasselbe Skript aus.

Entscheidung: _offen_

## F7 – Ganzzahlen in den Traits

Vorzeichenbehaftete Ganzzahlen implementieren im Entwurf bereits `Ring`.
Brauchen vorzeichenlose Typen eine eigene Stufe ohne Negation (`Semiring`)?

Empfehlung: erst bei Bedarf.

Entscheidung: _offen_

## F9 – Geometrie für Kollisionen

Gewebe-Kollision braucht Geometrie-Bausteine: Abstand Punkt–Dreieck und
Kante–Kante, Schnitt-Tests, kontinuierliche Kollisionserkennung (führt auf
kubische Gleichungen), räumliche Suchstrukturen. Gehört das in azmath oder
in ein eigenes Crate (z. B. `azgeom`) auf azmath?

Empfehlung: eigenes Crate. azmath liefert dafür Vektoren, Polynom-Nullstellen
und Intervall-Hilfen; Suchstrukturen sind keine reine Mathematik.

Entscheidung: _offen_

## F10 – Welches Verfahren ist mit dem „Durchbruch“ gemeint?

Vermutung: *Incremental Potential Contact* (IPC, 2020) für
durchdringungsfreie Kontakte bei Gewebe und Festkörpern; bei Flüssigkeiten
kämen u. a. Material-Point-Methoden in Frage. Die Antwort bestimmt, welche
Löser zuerst gebraucht werden (IPC: Newton mit Liniensuche, dünnbesetzte
symmetrische Löser, zweite Ableitungen über Autodiff).

Entscheidung: _offen_

## F11 – Finite Elemente ohne Heap

Unstrukturierte Netze haben zur Laufzeit bestimmte Größen. Vorschlag: Netz
und Systemmatrix als Sichten auf Puffer des Aufrufers (läuft ohne `alloc`),
bequeme besitzende Typen nur hinter `alloc`. Netzerzeugung (Meshing) ist kein
Ziel von azmath.

Entscheidung: _offen_
