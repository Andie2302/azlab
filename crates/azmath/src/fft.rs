//! Diskrete Fourier-Transformation.
//!
//! - [`fft`] und [`ifft`]: schnelle Transformation an Ort und Stelle für
//!   Längen, die Zweierpotenzen sind.
//! - [`dft`] und [`idft`]: direkte Summe für beliebige Längen, mit
//!   quadratischem Aufwand.
//!
//! Vereinbarung: Die Hintransformation ist
//! `X[k] = Σₙ x[n]·e^(−2πi·k·n/N)` ohne Vorfaktor, die Rücktransformation
//! hat den Vorfaktor `1/N`. Hin und zurück ergibt also das Ausgangssignal.

use crate::linalg::LinalgError;
use crate::number::Complex;
use crate::scalar::RealElementary;

/// Gemeinsamer Kern von [`fft`] und [`ifft`]; `sign` ist −1 bzw. +1.
fn transform<T: RealElementary>(data: &mut [Complex<T>], sign: T) -> Result<(), LinalgError> {
    let length = data.len();
    if length == 0 {
        return Ok(());
    }
    if !length.is_power_of_two() {
        return Err(LinalgError::DimensionMismatch);
    }
    // Bitumkehr-Vertauschung.
    let bits = length.trailing_zeros();
    for index in 0..length {
        let reversed = if bits == 0 { 0 } else { index.reverse_bits() >> (usize::BITS - bits) };
        if reversed > index {
            data.swap(index, reversed);
        }
    }
    // Schmetterlinge mit wachsender Blocklänge.
    let two_pi = T::PI * T::from_f64(2.0);
    let mut block = 2;
    while block <= length {
        let half = block / 2;
        let angle = sign * two_pi / T::from_f64(block as f64);
        for start in (0..length).step_by(block) {
            for offset in 0..half {
                // Drehfaktor direkt aus dem Winkel: kein Aufsummieren von
                // Rundungsfehlern.
                let twiddle = Complex::from_polar(T::ONE, angle * T::from_f64(offset as f64));
                let upper = data[start + offset];
                let lower = data[start + offset + half] * twiddle;
                data[start + offset] = upper + lower;
                data[start + offset + half] = upper - lower;
            }
        }
        block *= 2;
    }
    Ok(())
}

/// Schnelle Fourier-Transformation an Ort und Stelle.
///
/// Fehler, wenn die Länge keine Zweierpotenz ist; dann [`dft`] nehmen.
pub fn fft<T: RealElementary>(data: &mut [Complex<T>]) -> Result<(), LinalgError> {
    transform(data, -T::ONE)
}

/// Schnelle Rücktransformation an Ort und Stelle (mit Vorfaktor `1/N`).
///
/// Fehler, wenn die Länge keine Zweierpotenz ist; dann [`idft`] nehmen.
pub fn ifft<T: RealElementary>(data: &mut [Complex<T>]) -> Result<(), LinalgError> {
    transform(data, T::ONE)?;
    let scale = T::from_f64(data.len() as f64).recip();
    for value in data.iter_mut() {
        *value = *value * scale;
    }
    Ok(())
}

/// Gemeinsamer Kern von [`dft`] und [`idft`].
fn direct<T: RealElementary>(
    input: &[Complex<T>],
    output: &mut [Complex<T>],
    sign: T,
) -> Result<(), LinalgError> {
    let length = input.len();
    if output.len() != length {
        return Err(LinalgError::DimensionMismatch);
    }
    let step = sign * T::PI * T::from_f64(2.0) / T::from_f64(length.max(1) as f64);
    let zero = Complex::new(T::ZERO, T::ZERO);
    for (k, out) in output.iter_mut().enumerate() {
        *out = input.iter().enumerate().fold(zero, |sum, (n, &value)| {
            // k·n modulo N hält den Winkel klein und damit genau.
            let turns = (k * n) % length;
            sum + value * Complex::from_polar(T::ONE, step * T::from_f64(turns as f64))
        });
    }
    Ok(())
}

/// Diskrete Fourier-Transformation als direkte Summe, für jede Länge.
///
/// Fehler, wenn `output` nicht so lang ist wie `input`.
pub fn dft<T: RealElementary>(
    input: &[Complex<T>],
    output: &mut [Complex<T>],
) -> Result<(), LinalgError> {
    direct(input, output, -T::ONE)
}

/// Rücktransformation als direkte Summe (mit Vorfaktor `1/N`), für jede
/// Länge.
///
/// Fehler, wenn `output` nicht so lang ist wie `input`.
pub fn idft<T: RealElementary>(
    input: &[Complex<T>],
    output: &mut [Complex<T>],
) -> Result<(), LinalgError> {
    direct(input, output, T::ONE)?;
    let scale = T::from_f64(input.len().max(1) as f64).recip();
    for value in output.iter_mut() {
        *value = *value * scale;
    }
    Ok(())
}

#[cfg(all(test, any(feature = "std", feature = "libm")))]
mod tests {
    use super::{dft, fft, idft, ifft};
    use crate::linalg::LinalgError;
    use crate::number::Complex;
    use crate::scalar::{Conjugate, RealElementary};

    fn near<T: RealElementary>(a: &[Complex<T>], b: &[Complex<T>]) -> bool {
        let limit = T::EPSILON.sqrt();
        a.len() == b.len() && a.iter().zip(b).all(|(&x, &y)| (x - y).norm() <= limit)
    }

    fn check<T: RealElementary>() {
        let f = T::from_f64;
        let real = |v: f64| Complex::from_real(f(v));
        let zero = real(0.0);

        // Einheitsimpuls → konstantes Spektrum; Konstante → Impuls bei 0.
        let mut impulse = [zero; 8];
        impulse[0] = real(1.0);
        assert!(fft(&mut impulse).is_ok() && near(&impulse, &[real(1.0); 8]));
        let mut constant = [real(1.0); 8];
        let mut expected = [zero; 8];
        expected[0] = real(8.0);
        assert!(fft(&mut constant).is_ok() && near(&constant, &expected));

        // Reine Schwingung mit drei Perioden: Spitze genau bei k = 3.
        let mut wave: [Complex<T>; 16] = core::array::from_fn(|n| {
            Complex::from_polar(T::ONE, T::PI * f(2.0) * f(3.0) * f(n as f64) / f(16.0))
        });
        let mut expected = [zero; 16];
        expected[3] = real(16.0);
        assert!(fft(&mut wave).is_ok() && near(&wave, &expected));

        // Hin und zurück, Vergleich mit der direkten Summe, Parseval.
        let signal: [Complex<T>; 8] = core::array::from_fn(|n| {
            Complex::new(f(n as f64) - f(3.0), f(0.5) * f((n * n) as f64) - f(2.0))
        });
        let mut spectrum = signal;
        let mut direct = [zero; 8];
        assert!(fft(&mut spectrum).is_ok() && dft(&signal, &mut direct).is_ok());
        assert!(near(&spectrum, &direct));
        let energy = |values: &[Complex<T>]| values.iter().fold(T::ZERO, |sum, z| sum + z.norm_sqr());
        assert!((energy(&spectrum) - f(8.0) * energy(&signal)).abs() <= T::EPSILON.sqrt() * energy(&spectrum));
        let mut back = spectrum;
        assert!(ifft(&mut back).is_ok() && near(&back, &signal));

        // Beliebige Länge über die direkte Summe.
        let odd: [Complex<T>; 5] = core::array::from_fn(|n| Complex::new(f(n as f64), f(1.0)));
        let mut spectrum = [zero; 5];
        let mut back = [zero; 5];
        assert!(dft(&odd, &mut spectrum).is_ok() && idft(&spectrum, &mut back).is_ok());
        assert!(near(&back, &odd));
        // Erster Koeffizient ist die Summe.
        assert!((spectrum[0] - Complex::new(f(10.0), f(5.0))).norm() <= T::EPSILON.sqrt());

        // Sonderfälle und Fehler.
        let mut single = [real(3.0)];
        assert!(fft(&mut single).is_ok() && single[0] == real(3.0));
        let mut empty: [Complex<T>; 0] = [];
        assert!(fft(&mut empty).is_ok() && ifft(&mut empty).is_ok());
        assert!(dft(&empty, &mut []).is_ok() && idft(&empty, &mut []).is_ok());
        let mut three = [zero; 3];
        assert_eq!(fft(&mut three), Err(LinalgError::DimensionMismatch));
        assert_eq!(ifft(&mut three), Err(LinalgError::DimensionMismatch));
        assert_eq!(dft(&odd, &mut three), Err(LinalgError::DimensionMismatch));
        assert_eq!(idft(&odd, &mut three), Err(LinalgError::DimensionMismatch));
    }

    #[test]
    fn fft_f32() {
        check::<f32>();
    }

    #[test]
    fn fft_f64() {
        check::<f64>();
    }
}
