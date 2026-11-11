//! Interne Makros für die Operator-Implementierungen der Zahlentypen.

/// Leitet einen zusammengesetzten Zuweisungsoperator (`+=`, ...) aus dem
/// zugehörigen zweistelligen Operator ab.
macro_rules! impl_op_assign {
    ([$($gen:tt)*] $ty:ty, $tr:ident, $method:ident, $op:tt) => {
        impl<$($gen)*> core::ops::$tr for $ty {
            #[inline]
            fn $method(&mut self, rhs: Self) {
                *self = *self $op rhs;
            }
        }
    };
}
