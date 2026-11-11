//! Makro, mit dem jedes Backend seine Funktionen den Traits zuordnet.

/// Implementiert einen Skalar-Trait für einen Typ aus einer Zuordnungstabelle
/// `methode(weitere argumente) => backend_funktion;`.
macro_rules! impl_backend {
    ($tr:ident for $t:ty { $( $method:ident($($arg:ident),*) => $func:path; )* }) => {
        impl $crate::scalar::$tr for $t {
            $(
                #[inline]
                fn $method(self $(, $arg: Self)*) -> Self {
                    $func(self $(, $arg)*)
                }
            )*
        }
    };
}
