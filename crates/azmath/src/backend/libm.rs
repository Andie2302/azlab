//! Backend `libm`: reine Rust-Implementierung, läuft in `no_std` ohne `alloc`.

impl_backend!(Elementary for f32 {
    sqrt() => libm::sqrtf;
    exp() => libm::expf;
    exp2() => libm::exp2f;
    exp_m1() => libm::expm1f;
    ln() => libm::logf;
    log2() => libm::log2f;
    log10() => libm::log10f;
    ln_1p() => libm::log1pf;
    powf(n) => libm::powf;
    sin() => libm::sinf;
    cos() => libm::cosf;
    tan() => libm::tanf;
    asin() => libm::asinf;
    acos() => libm::acosf;
    atan() => libm::atanf;
    sinh() => libm::sinhf;
    cosh() => libm::coshf;
    tanh() => libm::tanhf;
    asinh() => libm::asinhf;
    acosh() => libm::acoshf;
    atanh() => libm::atanhf;
});

impl_backend!(RealElementary for f32 {
    cbrt() => libm::cbrtf;
    atan2(x) => libm::atan2f;
    hypot(other) => libm::hypotf;
    floor() => libm::floorf;
    ceil() => libm::ceilf;
    round() => libm::roundf;
    trunc() => libm::truncf;
    mul_add(a, b) => libm::fmaf;
});

impl_backend!(Elementary for f64 {
    sqrt() => libm::sqrt;
    exp() => libm::exp;
    exp2() => libm::exp2;
    exp_m1() => libm::expm1;
    ln() => libm::log;
    log2() => libm::log2;
    log10() => libm::log10;
    ln_1p() => libm::log1p;
    powf(n) => libm::pow;
    sin() => libm::sin;
    cos() => libm::cos;
    tan() => libm::tan;
    asin() => libm::asin;
    acos() => libm::acos;
    atan() => libm::atan;
    sinh() => libm::sinh;
    cosh() => libm::cosh;
    tanh() => libm::tanh;
    asinh() => libm::asinh;
    acosh() => libm::acosh;
    atanh() => libm::atanh;
});

impl_backend!(RealElementary for f64 {
    cbrt() => libm::cbrt;
    atan2(x) => libm::atan2;
    hypot(other) => libm::hypot;
    floor() => libm::floor;
    ceil() => libm::ceil;
    round() => libm::round;
    trunc() => libm::trunc;
    mul_add(a, b) => libm::fma;
});
