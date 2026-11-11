# Phase 2 – Zahlentypen (Detailplan)

Status: **umgesetzt** (Schritte A–E, je ein Commit). Von der KI nach jedem
Schritt mit Standard-Features und ohne Features geprüft (Clippy und Tests);
die volle Prüfmatrix über `scripts/pruefen.sh` steht noch aus.

## Abweichungen vom Plan bei der Umsetzung

- `Real` setzt jetzt `Commutative` voraus. Ohne das ließ sich nicht
  ausdrücken, dass die Verdopplung einer reellen Algebra ein Körper ist.
- `Dual` und `DualN` liegen gemeinsam in `src/number/dual/` und teilen sich
  die Implementierung über Makros; einen Ordner `dual_n/` gibt es nicht.
- `CayleyDickson` hat den Operator `/` (Rechtsdivision) auf jeder Stufe, also
  auch für Oktonionen; `Field` bleibt auf Stufe 1 und 2 beschränkt.
- `Quaternion` hat zusätzlich `inv()` und `sqrt()`.
- `Complex` hat zusätzlich `mul_i()` und `mul_neg_i()`.
- `Dual`/`DualN` implementieren `Conjugate` nicht; `constant`/`variable`
  brauchen nur `T: Algebra` statt `T: Field`.
- Die Umkehrfunktionen von `Complex` (`asin`, `acos`, ...) nutzen die
  Lehrbuchformeln und sind nahe den Verzweigungsschnitten nicht auf die
  letzte Stelle genau.
- Für Tests gibt es das interne Modul `src/testutil.rs`.

## Ziel

Die Zahlentypen oberhalb der primitiven Gleitkommazahlen: komplexe Zahlen,
Dualzahlen (Autodiff vorwärts), Cayley-Dickson-Algebren und Quaternionen.
Alle implementieren die Skalar-Traits, sodass spätere Algorithmen sie ohne
Sonderfälle verarbeiten.

## Nicht-Ziele dieser Phase

- Keine Vektor-/Matrix-Typen (Phase 3); Rotationen arbeiten auf `[T; 3]`.
- Keine Rotationsmatrizen, keine Euler-Winkel (Phase 3/5).
- Kein Rückwärts-Autodiff (Phase 10), keine Hilfsfunktionen
  `derivative`/`gradient`/`jacobian` (Phase 4).
- Kein eigener `HyperDual`-Typ; zweite Ableitungen über `Dual<Dual<T>>`.
- Keine Formatierung über `Display` hinaus, kein Parsen aus Text.

## Abhängigkeiten und Features

- Keine neue Abhängigkeit, kein neues Feature.
- Alles läuft in `no_std` ohne `alloc`.
- Im Bibliothekscode ist kein `cfg` nötig: Elementarfunktionen hängen nur an
  Trait-Schranken (`T: RealElementary`). Nur Tests mit `f32`/`f64` und
  Elementarfunktionen stehen hinter `cfg(any(feature = "std", feature = "libm"))`.

## Schritt A – Umbau der Skalar-Traits (ändert den Phase-1-Entwurf)

Neue Hierarchie:

```
Algebra ── Ring ── Field ──┬── Elementary ──┐
   │                       └── Real ────────┴── RealElementary
   └── Conjugate            Commutative (Marker)
```

```rust
/// Rechenoperationen mit Null und Eins, ohne Zusage der Assoziativität.
pub trait Algebra:
    Copy + Debug + PartialEq
    + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Neg<Output = Self>
    + AddAssign + SubAssign + MulAssign
{
    const ZERO: Self;
    const ONE: Self;
}

/// Marker: Multiplikation ist assoziativ.
pub trait Ring: Algebra {}

/// Marker: Multiplikation ist kommutativ.
pub trait Commutative: Algebra {}

pub trait Field: Ring + Div<Output = Self> + DivAssign {
    fn from_f64(value: f64) -> Self;   // neu hier, bisher in Real
    fn recip(self) -> Self { .. }
    fn powi(self, n: i32) -> Self { .. }
}

/// Konjugation und Norm; für reelle Typen ist `conj` die Identität.
pub trait Conjugate: Algebra {
    type Real: Real;
    fn conj(self) -> Self;
    fn norm_sqr(self) -> Self::Real;
    fn scale(self, factor: Self::Real) -> Self;
}
```

- `ZERO`/`ONE` wandern von `Ring` nach `Algebra`; `T::ZERO` funktioniert
  unverändert.
- `Real` verliert `from_f64` (kommt über `Field`), behält `to_f64`.
- Primitive: Ganzzahlen `Algebra + Ring + Commutative`; `f32`/`f64`
  zusätzlich `Field`, `Real`, `Conjugate<Real = Self>`.

## Schritt B – `Complex<T>`

```rust
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Complex<T> { pub re: T, pub im: T }
```

| Für | Angebot |
|-----|---------|
| `T: Real` | `new`, `I`, `from_real`; Operatoren `+ - * /` und Zuweisungsformen; `Mul<T>`/`Div<T>` (Skalierung); `Algebra`, `Ring`, `Commutative`, `Field`, `Conjugate<Real = T>` |
| `T: RealElementary` | `norm`, `arg`, `from_polar`, `to_polar`; `Elementary` (alle Funktionen, Hauptzweige, siehe F13) |
| `f32`, `f64` konkret | `T * Complex<T>` (links-skalar; generisch wegen Kohärenzregeln nicht möglich) |
| Umwandlung | `From<T>`, `From<(T, T)>`, von/nach `CayleyDickson<T>` |

Hinweise zur Numerik: Division nach Smith (kein vermeidbarer Überlauf),
`norm` über `hypot`, `sqrt` in der überlaufsicheren Form.

## Schritt C – `Dual<T>` und `DualN<T, N>`

```rust
#[derive(Clone, Copy, Debug, Default)]
pub struct Dual<T> { pub re: T, pub eps: T }

#[derive(Clone, Copy, Debug)]
pub struct DualN<T, const N: usize> { pub re: T, pub eps: [T; N] }
```

| Für | Angebot |
|-----|---------|
| `T: Field` | `constant(x)` (Ableitung 0), `variable(x)` (Ableitung 1); Operatoren; `Algebra`, `Ring`, `Field`; `Commutative`, wenn `T: Commutative` |
| `T: Elementary` | `Elementary` über die Kettenregel |
| `T: Real` | `Real`; `PartialEq`/`PartialOrd` vergleichen **nur `re`** |
| `T: RealElementary` | `RealElementary`; `floor`/`ceil`/`round`/`trunc` haben Ableitung 0 |
| `DualN` | wie `Dual`; `variable(x, index) -> Option<Self>` (kein Panic bei falschem Index) |

- `Dual<Complex<T>>` funktioniert (braucht nur `Elementary`).
- `Dual<Dual<T>>` liefert zweite Ableitungen.
- An Unstetigkeiten (`abs` bei 0, `sqrt` bei 0) gilt die Ableitung der
  IEEE-Rechnung (ggf. unendlich/NaN); wird dokumentiert, nicht abgefangen.

## Schritt D – `CayleyDickson<A>`

```rust
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct CayleyDickson<A> { pub a: A, pub b: A }

pub type Octonion<T> = CayleyDickson<CayleyDickson<CayleyDickson<T>>>;
```

Rechenregeln (mit `*` für Konjugation):

- `(a, b) + (c, d) = (a + c, b + d)`
- `(a, b) · (c, d) = (a·c − d*·b, d·a + b·c*)`
- `(a, b)* = (a*, −b)`
- `norm_sqr = a.norm_sqr() + b.norm_sqr()`
- `inv = conj / norm_sqr`

| Schranke | `CayleyDickson<A>` ist |
|----------|------------------------|
| `A: Conjugate` | `Algebra`, `Conjugate<Real = A::Real>`, Methode `inv()` |
| `A: Ring + Commutative + Conjugate` | zusätzlich `Ring`, `Field` (Division = Rechtsdivision, F12) |
| `A: Real + Conjugate` | zusätzlich `Commutative` |

Ergebnis: Stufe 1 (komplex) ist kommutativer Körper, Stufe 2 (Quaternion)
Schiefkörper, Stufe 3 (Oktonion) nur noch `Algebra` – genau wie in der
Mathematik, und der Compiler erzwingt es.

## Schritt E – `Quaternion<T>`

```rust
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Quaternion<T> { pub w: T, pub x: T, pub y: T, pub z: T }
```

| Für | Angebot |
|-----|---------|
| `T: Real` | `new`, `IDENTITY`, `I`/`J`/`K`, `from_real`, `from_parts(w, [x, y, z])`, `scalar()`, `vector() -> [T; 3]`; Operatoren inkl. `Mul<T>`/`Div<T>`; `Algebra`, `Ring`, `Field`, `Conjugate<Real = T>`; `dot`; `nlerp` braucht Wurzel → siehe unten |
| `T: RealElementary` | `norm`, `normalize() -> Option<Self>`, `exp`, `ln`, `powf(T)`, `from_axis_angle([T; 3], T) -> Option<Self>`, `to_axis_angle() -> ([T; 3], T)`, `rotate([T; 3]) -> [T; 3]`, `slerp`, `nlerp` |
| Umwandlung | von/nach `CayleyDickson<CayleyDickson<T>>` und `CayleyDickson<Complex<T>>`-ähnlicher Zerlegung `(Complex<T>, Complex<T>)` |

- Kein `Commutative`, kein `Elementary` (Potenz zweier Quaternionen ist nicht
  eindeutig); `exp`/`ln`/`powf` sind eigene Methoden.
- `rotate` setzt eine Einheitsquaternion voraus (dokumentiert). Ein eigener
  Typ für Einheitsquaternionen kommt als `SO3` in Phase 5.
- `slerp` nimmt den kürzeren Weg und fällt für fast gleiche Quaternionen auf
  `nlerp` zurück.

## Fehlerbehandlung

Keine Panics. Wo es scheitern kann (Normieren der Null, Achse der Länge 0,
Index außerhalb), ist der Rückgabetyp `Option`. Ein Fehler-Enum ist in dieser
Phase noch nicht nötig.

## Dateien

```
src/scalar/algebra.rs        Algebra, Commutative            (neu)
src/scalar/conjugate.rs      Conjugate                       (neu)
src/scalar/ring.rs  field.rs  real.rs  primitives.rs  tests.rs   (angepasst)
src/scalar.rs                Deklarationen/Re-Exporte        (angepasst)

src/number.rs                Deklarationen/Re-Exporte
src/number/complex.rs        Deklarationen
src/number/complex/repr.rs  ops.rs  traits.rs  elementary.rs  convert.rs  tests.rs
src/number/dual.rs
src/number/dual/repr.rs  ops.rs  traits.rs  elementary.rs  tests.rs
src/number/dual_n.rs
src/number/dual_n/repr.rs  ops.rs  traits.rs  elementary.rs  tests.rs
src/number/cayley_dickson.rs
src/number/cayley_dickson/repr.rs  ops.rs  traits.rs  tests.rs
src/number/quaternion.rs
src/number/quaternion/repr.rs  ops.rs  traits.rs  rotation.rs  convert.rs  tests.rs

src/lib.rs                   `pub mod number;` + Re-Exporte  (angepasst)
```

## Tests

Alle Erwartungswerte sind exakt von Hand herleitbar oder Identitäten.

- **Traits**: generische Körper-Axiome (Assoziativität, Distributivität,
  `x * x.recip() == 1` mit Toleranz) für jeden `Field`-Typ.
- **Complex**: `I * I == -1`; `exp(iπ) ≈ -1`; `sqrt(-1) == i`;
  `ln(exp(z)) ≈ z`; `sin² + cos² ≈ 1`; Polarform hin und zurück; Division
  großer Beträge ohne Überlauf.
- **Dual**: Ableitungen von `x³`, `1/x`, `sqrt`, `exp`, `ln`, `sin`, `cos`,
  `atan2` gegen die bekannte Formel; Produkt- und Kettenregel;
  `Dual<Dual<f64>>`: zweite Ableitung von `x³` ist `6x`; Vergleich ignoriert
  `eps`.
- **DualN**: Gradient von `x²·y` ist `(2xy, x²)`; falscher Index → `None`.
- **CayleyDickson**: Stufe 1 stimmt mit `Complex` überein, Stufe 2 mit
  `Quaternion`; Oktonion: `(e1·e2)·e4 ≠ e1·(e2·e4)` (nicht assoziativ),
  aber `norm(x·y) ≈ norm(x)·norm(y)` und `x·inv(x) ≈ 1`.
- **Quaternion**: `i² = j² = k² = ijk = −1`; Drehung um z um 90° bildet
  x-Achse auf y-Achse ab; zwei Drehungen = Produkt; `slerp` bei 0, ½, 1;
  `to_axis_angle(from_axis_angle(..))`; `ln(exp(q)) ≈ q`; Normieren der
  Null → `None`.
- Jeder Test läuft für `f32` und `f64`.

## Reihenfolge der Umsetzung

1. Schritt A (Trait-Umbau) – eigener Commit, danach Prüflauf.
2. Schritt B `Complex`.
3. Schritt C `Dual`, dann `DualN`.
4. Schritt D `CayleyDickson`.
5. Schritt E `Quaternion`.

Nach jedem Schritt ein Commit; Prüflauf durch den Menschen mindestens nach
Schritt A und am Ende.

## Festlegungen

- `/` ist Rechtsdivision `a * b⁻¹` (entschieden).
- Hauptzweige wie C99 (entschieden).
- `Dual`/`DualN`/`Complex`-Felder sind öffentlich (`re`, `eps`, `im`).
