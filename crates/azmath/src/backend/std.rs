//! Backend `std`: nutzt die Gleitkomma-Methoden der Standardbibliothek.

impl_backend!(Elementary for f32 {
    sqrt() => f32::sqrt;
    exp() => f32::exp;
    exp2() => f32::exp2;
    exp_m1() => f32::exp_m1;
    ln() => f32::ln;
    log2() => f32::log2;
    log10() => f32::log10;
    ln_1p() => f32::ln_1p;
    powf(n) => f32::powf;
    sin() => f32::sin;
    cos() => f32::cos;
    tan() => f32::tan;
    asin() => f32::asin;
    acos() => f32::acos;
    atan() => f32::atan;
    sinh() => f32::sinh;
    cosh() => f32::cosh;
    tanh() => f32::tanh;
    asinh() => f32::asinh;
    acosh() => f32::acosh;
    atanh() => f32::atanh;
});

impl_backend!(RealElementary for f32 {
    cbrt() => f32::cbrt;
    atan2(x) => f32::atan2;
    hypot(other) => f32::hypot;
    floor() => f32::floor;
    ceil() => f32::ceil;
    round() => f32::round;
    trunc() => f32::trunc;
    mul_add(a, b) => f32::mul_add;
});

impl_backend!(Elementary for f64 {
    sqrt() => f64::sqrt;
    exp() => f64::exp;
    exp2() => f64::exp2;
    exp_m1() => f64::exp_m1;
    ln() => f64::ln;
    log2() => f64::log2;
    log10() => f64::log10;
    ln_1p() => f64::ln_1p;
    powf(n) => f64::powf;
    sin() => f64::sin;
    cos() => f64::cos;
    tan() => f64::tan;
    asin() => f64::asin;
    acos() => f64::acos;
    atan() => f64::atan;
    sinh() => f64::sinh;
    cosh() => f64::cosh;
    tanh() => f64::tanh;
    asinh() => f64::asinh;
    acosh() => f64::acosh;
    atanh() => f64::atanh;
});

impl_backend!(RealElementary for f64 {
    cbrt() => f64::cbrt;
    atan2(x) => f64::atan2;
    hypot(other) => f64::hypot;
    floor() => f64::floor;
    ceil() => f64::ceil;
    round() => f64::round;
    trunc() => f64::trunc;
    mul_add(a, b) => f64::mul_add;
});
