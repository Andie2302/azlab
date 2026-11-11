//! Interne Makros für die Operator-Implementierungen der Zahlentypen.

/// Leitet einen zusammengesetzten Zuweisungsoperator (`+=`, ...) aus dem
/// zugehörigen zweistelligen Operator ab. Optional mit `where`-Schranken.
macro_rules! impl_op_assign {
    ([$($gen:tt)*] $ty:ty, $tr:ident, $method:ident, $op:tt) => {
        impl_op_assign!([$($gen)*] $ty, $tr, $method, $op, where);
    };
    ([$($gen:tt)*] $ty:ty, $tr:ident, $method:ident, $op:tt, where $($bounds:tt)*) => {
        impl<$($gen)*> core::ops::$tr for $ty
        where
            $($bounds)*
        {
            #[inline]
            fn $method(&mut self, rhs: Self) {
                *self = *self $op rhs;
            }
        }
    };
}
