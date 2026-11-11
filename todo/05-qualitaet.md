# Qualität: Tests, Prüfmatrix, Doku

## Wer prüft

Die Prüfmatrix läuft mit

```bash
scripts/pruefen.sh
```

und schreibt nach `pruefung/ergebnis.txt`. Jeder Stand auf `master` besteht
sie. Die Beispiele der README werden als Tests mit übersetzt und ausgeführt
(`src/readme.rs`).

## Prüfmatrix

Jede Kombination muss bauen, Tests bestehen und Clippy-frei sein:

| Kombination | Zweck |
|-------------|-------|
| `--no-default-features` | Regel 1: Kern ohne alles |
| `--no-default-features --features libm` | `no_std` mit Backend |
| `--no-default-features --features std` | std-Backend |
| `--no-default-features --features alloc` | Heap ohne Backend |
| `--all-features` | Rangfolge der Backends |

Zusätzlich ein Build für ein Ziel ohne Betriebssystem (z. B.
`thumbv7em-none-eabihf`) mit und ohne `libm`: beweist echtes `no_std`.

## Tests

- **Generische Testfunktionen**: ein Test ist eine Funktion über `T`, die für
  `f32`, `f64`, `Dual<f64>` usw. instanziiert wird. Neue Zahlentypen
  durchlaufen so automatisch alle vorhandenen Tests.
- **Backend-Gleichheit**: dieselben Tests laufen mit jedem Backend; Toleranz
  als Vielfaches von `EPSILON`, nie exakter Vergleich transzendenter Werte.
- **Referenzwerte** nur aus zitierten Quellen (Norm, Lehrbuch, Tabellenwerk),
  Quelle steht im Test.
- **Autodiff** gegen analytische Ableitungen prüfen.
- **Löser** über Rückeinsetzen prüfen (Residuum), nicht nur gegen Zielwerte.
- Tests dürfen `unwrap`/`panic!` nicht verwenden (Workspace-Lints gelten auch
  dort); Fehler mit `assert!` auf `Result` prüfen.

## Doku

- Jedes öffentliche Element hat einen deutschen Doku-Kommentar
  (`missing_docs`).
- Jedes Modul nennt: was es löst, welche Traits es fordert, ob es Heap
  braucht.
- Verfahren nennen Quelle und Voraussetzungen (z. B. „Matrix symmetrisch
  positiv definit“).
