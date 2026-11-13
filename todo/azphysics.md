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

1. Freie Bewegung unter Schwerkraft und Dämpfung.
2. Abstandsbedingungen (XPBD) für Dehnung, Scherung, Biegung.
3. Hindernisse Ebene und Kugel: hinausschieben, Reibung.
4. Abstoßung Teilchen–Dreieck innerhalb der Dicke; unelastisch, d. h. sie
   bremst die Annäherung, schleudert aber nicht weg.
5. Sicherung: CCD für jedes nahe Paar Punkt–Dreieck und Kante–Kante. Bei
   einer Berührung wird die Bewegung der vier Teilchen kurz vor dem
   Berührzeitpunkt beendet (80 %, dann 40 % der Zeit, sonst zurück an den
   Anfang). Wiederholen, bis ein Durchlauf nichts mehr findet; nach vier
   Durchläufen nur noch hartes Zurücksetzen, nach acht alle Teilchen zurück.
6. Geschwindigkeit aus der erreichten Lage.

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
- Abstoßung nur Teilchen–Dreieck, nicht Kante–Kante.
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

- Gewebe: Abstoßung Kante–Kante, Stoßgesetz mit Rückprall, Reibung zwischen
  Lagen, Dreiecksnetze beliebiger Form (bisher Rechteckgitter), bewegte
  Hindernisse aus Dreiecken.
- Kopplung Gewebe–Strömung.
- Strömung im Raum; freie Oberfläche auf dem Gitter (FLIP/PIC).
- Starre Körper.
