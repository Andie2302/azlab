# azphysics und azviz

Stand: erste Fassung umgesetzt und geprüft.

## Aufbau

| Modul | Datei(en) | Inhalt | Braucht |
|-------|-----------|--------|---------|
| `geometry` | `geometry.rs` | `Aabb`, nächste Punkte Punkt–Strecke, Punkt–Dreieck, Strecke–Strecke, Abstände, Schnitt Strecke–Dreieck | – |
| `ccd` | `ccd.rs` | `vertex_triangle`, `edge_edge`: frühester Berührzeitpunkt bewegter Teile | – |
| `broadphase` | `broadphase.rs` | `SpatialHash`: Raumgitter mit Zählsortierung | `alloc` |
| `cloth` | `cloth/model.rs`, `step.rs`, `collide.rs` | `Cloth`, `ClothParams`, `Collider`, `StepStats` | `alloc` |
| `fluid` | `fluid/grid.rs`, `advect.rs`, `pressure.rs` | `GridFluid`: Navier-Stokes in der Ebene | `alloc` |
| `fluid` | `fluid/particles.rs` | `ParticleFluid`: Flüssigkeit mit freier Oberfläche | `alloc` |

azviz: `Canvas` (Zeichenfläche), `Scene` + `run` (Fenster, `--png`, `--gif`),
`Plane` (Ausschnitt der Ebene), `Camera` + `Batch` (schattierte Dreiecke im
Raum, mit Tiefenprüfung je Bildpunkt).

## Gewebe: Ablauf eines Teilschritts

1. Freie Bewegung unter Schwerkraft und Dämpfung. Die Lage wird
   kompensiert weitergerückt: Was beim Runden verloren geht, bleibt als
   Rest je Teilchen erhalten (siehe „Lage und Rest“).
2. Abstandsbedingungen (XPBD) für Dehnung, Scherung, Biegung, in
   Durchläufen paarweise vor und zurück (siehe „Durchlaufrichtung“).
3. Hindernisse Ebene und Kugel: hinausschieben, Coulomb-Reibung.
4. Abstoßung Teilchen–Dreieck innerhalb der Dicke; unelastisch, d. h. sie
   bremst die Annäherung, schleudert aber nicht weg. Dazu Coulomb-Reibung,
   auch zwischen den Lagen des Tuchs.
5. Sicherung: CCD für jedes nahe Paar Punkt–Dreieck und Kante–Kante. Bei
   einer Berührung wird die Bewegung der vier Teilchen kurz vor dem
   Berührzeitpunkt beendet (80 %, dann 40 % der Zeit, sonst zurück an den
   Anfang). Wiederholen, bis ein Durchlauf nichts mehr findet; nach vier
   Durchläufen nur noch hartes Zurücksetzen, nach acht alle Teilchen zurück.
6. Geschwindigkeit: die der freien Bewegung, geändert um das, was die
   Korrekturen des Teilschritts an der Lage verschoben haben.

### Durchlaufrichtung

Die Verbindungen werden nacheinander korrigiert. Jede einzelne Korrektur
verkleinert den Fehler; ihre Hintereinanderschaltung in stets derselben
Reihenfolge ist als Ganzes aber unsymmetrisch und kann Fehler im Kreis
drehen. Zusammen mit der Geschwindigkeit aus der Lageänderung wächst das an:
Ein feines Tuch (ab etwa 24 Teilchen je Kante), das ruhig in seiner Ruhelage
hing, war nach einer halben Sekunde mehrere Meter je Sekunde schnell und
faltete sich in der eigenen Ebene zusammen oder lief bis zum Überlauf davon.

Ein Durchlauf vorwärts und einer rückwärts sind zusammen symmetrisch und
verkleinern jeden Fehler. Die Durchläufe laufen deshalb paarweise; eine
ungerade Zahl in `ClothParams::iterations` wird aufgerundet. Der Wechsel der
Richtung nur von Teilschritt zu Teilschritt genügt nicht.

Test: `cloth::tests::hanging_cloth_stays_at_rest`.

### Reibung

Coulomb, an Ebene, Kugel, Dreiecken und zwischen Tuchlagen: Die Strecke, um
die ein Teilchen aus dem Hindernis geschoben wird, steht für die
Anpresskraft; das Gleiten längs der Oberfläche wird um höchstens Reibwert
mal diese Strecke gekürzt. Gleitet das Teilchen weniger weit, haftet es. Auf
der schiefen Ebene heißt das: liegen bleiben unter der Grenzsteigung, sonst
rutschen mit `g·(sin α − µ·cos α)` – unabhängig von der Schrittweite.

Die frühere Reibung kürzte das Gleiten um einen festen Anteil je
Teilschritt. Nichts haftete; ein Tuch kroch von der Kugel, und wie schnell,
hing von der Schrittweite ab. An Dreiecken gab es gar keine Reibung.

Tests: `friction_on_a_slope` (gegen die Formel, auf ein Prozent),
`friction_on_triangles`.

### Was das Verfahren nicht kann

Gemessen, kein Fehler, aber gut zu wissen:

- **Dehnung.** Mit wenigen Durchläufen gibt das Tuch unter Last nach: ein
  hängendes Tuch aus 24 Teilchen je Kante um 1,6 %, aus 40 um 4,7 %
  (8 Teilschritte, 2 Durchläufe). Doppelt so viele Teilschritte halbieren
  das etwa; mehr Durchläufe helfen weniger. Beim Aufprall sind kurzzeitig
  30 bis 40 % möglich.
- **Drehung wird gebremst.** Ein frei rotierendes Tuch verliert Drehimpuls,
  bei 3 rad/s und 8 Teilschritten 1,8 % je Sekunde. Der Verlust ist
  proportional zur Schrittweite und unabhängig von Durchläufen und
  Auflösung; Impuls und Schwerpunkt bleiben exakt.
- **Symmetrie.** Ein mittig auf die Kugel fallendes Tuch bleibt nicht
  spiegelbildlich: Die Reihenfolge der Korrekturen bricht die Symmetrie, und
  Falten verstärken kleine Unterschiede. Eine bevorzugte Richtung gibt es
  nicht; der Schwerpunkt bleibt auf wenige Millimeter in der Mitte.

### Lage und Rest

Die Lage eines Teilchens ist auf das Raster ihres Zahltyps gerundet, und
dieses Raster wird gröber, je weiter das Teilchen vom Ursprung entfernt ist
(`f32`: 6·10⁻⁸ m bei 1 m, 10⁻⁶ m bei 10 m). Der Weg eines Teilschritts ist
oft nur wenige tausend Rasterschritte lang. Zwei Folgen, beide systematisch
und nicht zufällig:

- Beim Weiterrücken geht der Rest des Wegs verloren, in jedem Teilschritt
  derselbe. Das Teilchen ist dauerhaft etwas zu schnell oder zu langsam.
- Wird die Geschwindigkeit aus dem gerundeten Weg zurückgerechnet, rastet
  sie auf Vielfache von Rasterschritt durch Teilschritt ein. Nachbarn mit
  verschiedenem Raster bekommen verschiedene Geschwindigkeiten, und die
  Verbindungen gleichen den Scheinfehler mit Stößen in Rastergröße aus.

Ein frei fliegendes Tuch trieb deshalb in jeder Richtung ab, in der es sich
bewegte; eine Richtung ohne Bewegung blieb exakt. Deshalb:

- `shift` rückt die Lage als fehlerfreie Summe weiter und merkt sich den
  Rest (`carry`). Lage und Rest zusammen sind die genaue Lage.
- Die Verbindungen rechnen mit der genauen Lage und verschieben über
  `shift`.
- Die Geschwindigkeit wird nicht aus dem Weg zurückgerechnet, sondern nur
  um die Korrekturen des Teilschritts geändert. Ein Teilchen, das nichts
  korrigiert hat, behält sie bitgenau.

Hindernisse, Abstoßung und Kollisionserkennung lesen nur die gerundete
Lage; der Rest ist kleiner als ein Rasterschritt und für sie ohne Belang.

Test: `cloth::tests::free_flight_does_not_drift` (250 m vom Ursprung, alle
Korrekturen eingeschaltet).

### Warum nichts durchdringt

- Zwei Dreiecksnetze können sich nur durchdringen, indem ein Punkt durch
  ein Dreieck oder eine Kante durch eine Kante tritt. Vorher sinkt der
  Abstand des Paars unter die Toleranz.
- Der Abstand wird über den ganzen Teilschritt verfolgt (`advance`): Er
  ändert sich je Zeit höchstens um die zurückgelegten Wege, also wird
  jeweils um `(Abstand − Toleranz) / Weg` vorgerückt. Das übersieht keine
  Annäherung, auch keine, bei der die Teile nie in einer Ebene liegen. Eine
  frühere Fassung suchte nur die Zeitpunkte, an denen die vier Punkte in
  einer Ebene liegen; sie ließ Annäherungen unter die Toleranz ohne
  Durchtritt zu, und solche Paare blieben danach aneinander stehen.
- Das Ergebnis (`ccd::Contact`) kennt drei Fälle: `Free` (sicher keine
  Berührung), `Touch(t)` (sicher Berührung zur Zeit `t`) und `Undecided(t)`
  (bis `t` sicher frei, danach nach 64 Schritten nicht entschieden –
  schleichende Annäherung). Das Gewebe behandelt `Undecided` wie eine
  Berührung und kürzt die Bewegung; `StepStats::ccd_undecided` zählt mit.
  In den Beispielen kommt der Fall nicht vor.
- Gerechnet wird mit der Relativbewegung (mittlere Verschiebung abgezogen).
  Ein Tuch im freien Flug stößt so nicht scheinbar mit sich selbst zusammen.
- Bewegungen werden nur gekürzt. Der letzte Durchlauf hat alle Paare frei
  gefunden, oder nichts hat sich bewegt.
- Die Kandidatenpaare stammen aus Nachbarschaftslisten mit
  Sicherheitsabstand. Sie werden neu gebaut, sobald ein Teilchen seinen
  Bereich um mehr als den halben Sicherheitsabstand verlässt; ein Paar, das
  nicht in der Liste steht, kann sich nicht berühren.

Tests dazu: `ccd::tests` (darunter 5000 zufällige Durchstöße, keiner
übersehen), `cloth::tests::tunnelling` (300 m/s gegen ein Dreieck; das Tuch
kommt an und bleibt davor; im freien Flug greift die Sicherung nicht ein),
`cloth::tests::folding` (nach jedem Schritt `intersections() == 0`).

### Grenzen

- Die Zusage gilt für Dreiecke und für geradlinige Bewegung im Teilschritt.
  Ebene und Kugel sind nicht Teil der Sicherung.
- Teile innerhalb der Toleranz (ein Zehntel der Dicke) gelten als berührt
  und bewegen sich nicht mehr gegeneinander, bis die Abstoßung sie trennt.
- Abstoßung und Reibung nur Teilchen–Dreieck, nicht Kante–Kante.
- Die Sicherung hält an, statt abprallen zu lassen: kein Stoßgesetz.

## Strömung

`GridFluid` (versetztes Gitter, *stable fluids*): Auftrieb → Transport der
Geschwindigkeit (Semi-Lagrange, wahlweise MacCormack) → Zähigkeit (implizit,
Gauß-Seidel) → Druckkorrektur → Transport des Rauchs. Die Druckgleichung
löst `azmath::solve::preconditioned_conjugate_gradient` mit unvollständiger
Cholesky-Zerlegung (MIC(0)) als Vorkonditionierer. Ränder je Seite: Wand,
Einlass, offen. Hindernisse sind ganze Zellen.

`ParticleFluid` (*position based fluids*): Dichte aus den Nachbarn, Korrektur
der Lagen bis zur Ruhedichte, nur Druck (kein Zug), Abstoßung gegen
Verklumpen, XSPH-Zähigkeit. Wände und Hindernisse sind ruhende Teilchen in
drei Lagen.

Die Wandteilchen haben eine eigene Dichte und einen eigenen Druckfaktor. Das
ist nötig, damit Wasser an einer Wand zur Ruhe kommt: Spiegelt die Wand nur
den Druck ihres Gegenübers, trägt sie kein Druckgefälle mit der Tiefe. Der
Spalte an der Wand fehlt dann etwa die halbe Stütze gegen die Schwerkraft,
sie sinkt mit rund 0,3 m/s ab, und es bleibt eine Umwälzung an jeder Wand
(Test `fluid::tests::pool`). Drei Lagen braucht die innerste Lage, um selbst
die volle Dichte zu sehen; Hindernisse sollten deshalb mindestens drei
Teilchenabstände dick sein.

Grenzen: nur in der Ebene; erster Ordnung in der Zeit; keine
Oberflächenspannung; Teilchenflüssigkeit ohne Kopplung an das Gitter.

## Möglicher Ausbau

- Gewebe: Abstoßung Kante–Kante, Stoßgesetz mit Rückprall, steifere
  Verbindungen bei gleicher Rechenzeit, Dreiecksnetze beliebiger Form (bisher Rechteckgitter), bewegte
  Hindernisse aus Dreiecken.
- Kopplung Gewebe–Strömung.
- Strömung im Raum; freie Oberfläche auf dem Gitter (FLIP/PIC).
- Starre Körper.
