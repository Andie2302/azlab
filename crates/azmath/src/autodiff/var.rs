//! Die Variable des Rückwärtsmodus.

use super::tape::{DETACHED, Record, TapeError};
use crate::linalg::Vector;
use crate::scalar::{
    Algebra, Commutative, Conjugate, Elementary, Field, Real, RealElementary, Ring,
};
use core::cmp::Ordering;
use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};

/// Zahl, deren Rechenschritte auf einem Band aufgezeichnet werden.
///
/// `Var` implementiert alle Skalar-Traits. Jeder generische Algorithmus
/// läuft deshalb unverändert mit `Var`, und [`gradient`](Self::gradient)
/// liefert danach die Ableitungen des Ergebnisses nach beliebig vielen
/// Eingängen in einem einzigen Rückwärtsdurchlauf. Das lohnt sich gegenüber
/// den Dualzahlen, wenn es viele Eingänge und wenige Ausgänge gibt.
///
/// Vergleiche (`==`, `<`, …) sehen nur den Wert. Konstanten (`ZERO`, `ONE`,
/// `from_f64`, [`constant`](Self::constant)) belegen keinen Platz auf dem
/// Band.
pub struct Var<'t, T, R> {
    value: T,
    index: usize,
    tape: Option<&'t R>,
}

impl<T: Copy, R> Clone for Var<'_, T, R> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<T: Copy, R> Copy for Var<'_, T, R> {}

impl<T: fmt::Debug, R> fmt::Debug for Var<'_, T, R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("Var").field(&self.value).finish()
    }
}

impl<T: PartialEq, R> PartialEq for Var<'_, T, R> {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl<T: PartialOrd, R> PartialOrd for Var<'_, T, R> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.value.partial_cmp(&other.value)
    }
}

impl<'t, T: Algebra, R: Record<T>> Var<'t, T, R> {
    /// Neue unabhängige Variable auf dem Band `tape`.
    pub fn new(tape: &'t R, value: T) -> Self {
        let index = tape.push([DETACHED; 2], [T::ZERO; 2]);
        Self { value, index, tape: Some(tape) }
    }

    /// Konstante: nimmt an der Rechnung teil, hat aber keine Ableitung.
    pub const fn constant(value: T) -> Self {
        Self { value, index: DETACHED, tape: None }
    }

    /// Der Wert.
    pub fn value(self) -> T {
        self.value
    }

    /// Ergebnis eines Schritts mit einem Vorgänger.
    fn unary(self, value: T, partial: T) -> Self {
        match self.tape {
            Some(tape) if self.index != DETACHED => {
                Self { value, index: tape.push([self.index, DETACHED], [partial, T::ZERO]), tape: Some(tape) }
            }
            _ => Self { value, ..self },
        }
    }

    /// Ergebnis eines Schritts mit zwei Vorgängern.
    fn binary(self, other: Self, value: T, partial_self: T, partial_other: T) -> Self {
        match self.tape.or(other.tape) {
            Some(tape) => Self {
                value,
                index: tape.push([self.index, other.index], [partial_self, partial_other]),
                tape: Some(tape),
            },
            None => Self::constant(value),
        }
    }

    /// Ableitungen dieses Werts nach den Variablen `wrt`, in einem
    /// Rückwärtsdurchlauf.
    ///
    /// Fehler, wenn das Band zu klein war oder eine Variable nicht zu ihm
    /// gehört.
    pub fn gradient<const K: usize>(self, wrt: [Self; K]) -> Result<[T; K], TapeError> {
        let mut out = [T::ZERO; K];
        let Some(tape) = self.tape.or_else(|| wrt.iter().find_map(|v| v.tape)) else {
            // Nur Konstanten: alle Ableitungen null.
            return Ok(out);
        };
        tape.backward(self.index, &wrt.map(|v| v.index), &mut out)?;
        Ok(out)
    }

    /// Ableitungen nach allen Komponenten eines Vektors von Variablen.
    pub fn gradient_vector<const N: usize>(
        self,
        wrt: Vector<Self, N>,
    ) -> Result<Vector<T, N>, TapeError> {
        self.gradient(wrt.0).map(Vector)
    }
}

impl<T: Algebra, R: Record<T>> Add for Var<'_, T, R> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        self.binary(rhs, self.value + rhs.value, T::ONE, T::ONE)
    }
}

impl<T: Algebra, R: Record<T>> Sub for Var<'_, T, R> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self.binary(rhs, self.value - rhs.value, T::ONE, -T::ONE)
    }
}

impl<T: Algebra, R: Record<T>> Neg for Var<'_, T, R> {
    type Output = Self;

    fn neg(self) -> Self {
        self.unary(-self.value, -T::ONE)
    }
}

impl<T: Algebra, R: Record<T>> Mul for Var<'_, T, R> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        self.binary(rhs, self.value * rhs.value, rhs.value, self.value)
    }
}

impl<T: Field, R: Record<T>> Div for Var<'_, T, R> {
    type Output = Self;

    // Quotientenregel: Die Ableitungen enthalten Produkte und Differenzen.
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: Self) -> Self {
        let quotient = self.value / rhs.value;
        self.binary(rhs, quotient, rhs.value.recip(), -quotient / rhs.value)
    }
}

impl_op_assign!(['t, T, R] Var<'t, T, R>, AddAssign, add_assign, +, where T: Algebra, R: Record<T>);
impl_op_assign!(['t, T, R] Var<'t, T, R>, SubAssign, sub_assign, -, where T: Algebra, R: Record<T>);
impl_op_assign!(['t, T, R] Var<'t, T, R>, MulAssign, mul_assign, *, where T: Algebra, R: Record<T>);
impl_op_assign!(['t, T, R] Var<'t, T, R>, DivAssign, div_assign, /, where T: Field, R: Record<T>);

impl<T: Algebra, R: Record<T>> Algebra for Var<'_, T, R> {
    const ZERO: Self = Self::constant(T::ZERO);
    const ONE: Self = Self::constant(T::ONE);
}

impl<T: Ring, R: Record<T>> Ring for Var<'_, T, R> {}

impl<T: Commutative, R: Record<T>> Commutative for Var<'_, T, R> {}

impl<T: Field, R: Record<T>> Field for Var<'_, T, R> {
    fn from_f64(value: f64) -> Self {
        Self::constant(T::from_f64(value))
    }
}

impl<T: Real, R: Record<T>> Conjugate for Var<'_, T, R> {
    type Real = Self;

    fn conj(self) -> Self {
        self
    }

    fn norm_sqr(self) -> Self {
        self * self
    }

    fn real_part(self) -> Self {
        self
    }

    fn abs_sum(self) -> Self {
        <Self as Real>::abs(self)
    }

    fn scale(self, factor: Self) -> Self {
        self * factor
    }
}

/// Alle Abfragen und Vergleiche beziehen sich nur auf den Wert.
impl<T: Real, R: Record<T>> Real for Var<'_, T, R> {
    const PI: Self = Self::constant(T::PI);
    const E: Self = Self::constant(T::E);
    const EPSILON: Self = Self::constant(T::EPSILON);
    const INFINITY: Self = Self::constant(T::INFINITY);
    const NEG_INFINITY: Self = Self::constant(T::NEG_INFINITY);
    const NAN: Self = Self::constant(T::NAN);

    fn to_f64(self) -> f64 {
        self.value.to_f64()
    }

    fn abs(self) -> Self {
        if self.value < T::ZERO { -self } else { self }
    }

    fn min(self, other: Self) -> Self {
        if other.value < self.value || self.value.is_nan() { other } else { self }
    }

    fn max(self, other: Self) -> Self {
        if other.value > self.value || self.value.is_nan() { other } else { self }
    }

    fn copysign(self, sign: Self) -> Self {
        let value = self.value.copysign(sign.value);
        if value == self.value { self } else { self.unary(value, -T::ONE) }
    }

    fn next_up(self) -> Self {
        self.unary(self.value.next_up(), T::ONE)
    }

    fn next_down(self) -> Self {
        self.unary(self.value.next_down(), T::ONE)
    }

    fn is_finite(self) -> bool {
        self.value.is_finite()
    }

    fn is_nan(self) -> bool {
        self.value.is_nan()
    }
}

impl<T: Elementary + Real, R: Record<T>> Elementary for Var<'_, T, R> {
    fn sqrt(self) -> Self {
        let value = self.value.sqrt();
        self.unary(value, (value + value).recip())
    }

    fn exp(self) -> Self {
        let value = self.value.exp();
        self.unary(value, value)
    }

    fn exp2(self) -> Self {
        let value = self.value.exp2();
        self.unary(value, value * T::from_f64(core::f64::consts::LN_2))
    }

    fn exp_m1(self) -> Self {
        self.unary(self.value.exp_m1(), self.value.exp())
    }

    fn ln(self) -> Self {
        self.unary(self.value.ln(), self.value.recip())
    }

    fn log2(self) -> Self {
        let scale = T::from_f64(core::f64::consts::LN_2);
        self.unary(self.value.log2(), (self.value * scale).recip())
    }

    fn log10(self) -> Self {
        let scale = T::from_f64(core::f64::consts::LN_10);
        self.unary(self.value.log10(), (self.value * scale).recip())
    }

    fn ln_1p(self) -> Self {
        self.unary(self.value.ln_1p(), (T::ONE + self.value).recip())
    }

    /// Ist der Exponent eine Konstante, wird `ln` der Basis nicht gebraucht;
    /// so bleibt z. B. `x²` auch für negative `x` differenzierbar.
    fn powf(self, n: Self) -> Self {
        let value = self.value.powf(n.value);
        let by_base = n.value * self.value.powf(n.value - T::ONE);
        let by_exponent = if n.index == DETACHED { T::ZERO } else { value * self.value.ln() };
        self.binary(n, value, by_base, by_exponent)
    }

    fn sin(self) -> Self {
        self.unary(self.value.sin(), self.value.cos())
    }

    fn cos(self) -> Self {
        self.unary(self.value.cos(), -self.value.sin())
    }

    fn tan(self) -> Self {
        let value = self.value.tan();
        self.unary(value, T::ONE + value * value)
    }

    fn asin(self) -> Self {
        self.unary(self.value.asin(), (T::ONE - self.value * self.value).sqrt().recip())
    }

    fn acos(self) -> Self {
        self.unary(self.value.acos(), -(T::ONE - self.value * self.value).sqrt().recip())
    }

    fn atan(self) -> Self {
        self.unary(self.value.atan(), (T::ONE + self.value * self.value).recip())
    }

    fn sinh(self) -> Self {
        self.unary(self.value.sinh(), self.value.cosh())
    }

    fn cosh(self) -> Self {
        self.unary(self.value.cosh(), self.value.sinh())
    }

    fn tanh(self) -> Self {
        let value = self.value.tanh();
        self.unary(value, T::ONE - value * value)
    }

    fn asinh(self) -> Self {
        self.unary(self.value.asinh(), (self.value * self.value + T::ONE).sqrt().recip())
    }

    fn acosh(self) -> Self {
        self.unary(self.value.acosh(), (self.value * self.value - T::ONE).sqrt().recip())
    }

    fn atanh(self) -> Self {
        self.unary(self.value.atanh(), (T::ONE - self.value * self.value).recip())
    }
}

impl<T: RealElementary, R: Record<T>> RealElementary for Var<'_, T, R> {
    fn cbrt(self) -> Self {
        let value = self.value.cbrt();
        self.unary(value, (T::from_f64(3.0) * value * value).recip())
    }

    fn atan2(self, x: Self) -> Self {
        let (y0, x0) = (self.value, x.value);
        let denominator = x0 * x0 + y0 * y0;
        self.binary(x, y0.atan2(x0), x0 / denominator, -y0 / denominator)
    }

    fn hypot(self, other: Self) -> Self {
        let value = self.value.hypot(other.value);
        self.binary(other, value, self.value / value, other.value / value)
    }

    /// Stückweise konstant: Ableitung null.
    fn floor(self) -> Self {
        Self::constant(self.value.floor())
    }

    /// Stückweise konstant: Ableitung null.
    fn ceil(self) -> Self {
        Self::constant(self.value.ceil())
    }

    /// Stückweise konstant: Ableitung null.
    fn round(self) -> Self {
        Self::constant(self.value.round())
    }

    /// Stückweise konstant: Ableitung null.
    fn trunc(self) -> Self {
        Self::constant(self.value.trunc())
    }

    fn mul_add(self, a: Self, b: Self) -> Self {
        self * a + b
    }
}
