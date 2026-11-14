//! Unabhängige Arbeit in Stücken, mit dem Feature `parallel` auf mehrere
//! Prozessorkerne verteilt.
//!
//! Das Ergebnis ist in beiden Fällen dasselbe, Bit für Bit: Die Stücke
//! werden in ihrer ursprünglichen Reihenfolge aneinandergehängt.

use alloc::vec::Vec;
use core::ops::Range;

/// Zahl der Elemente je Stück. Groß genug, dass sich das Verteilen lohnt,
/// klein genug, dass alle Kerne gleichmäßig zu tun haben.
#[cfg(feature = "parallel")]
const PIECE: usize = 512;

/// Ruft `work` für die Elemente `0..count` auf, stückweise; `work` hängt an
/// seine Liste, was es findet. Zurück kommen alle Funde in der Reihenfolge
/// der Elemente.
///
/// `work` muss für einen Bereich dasselbe liefern wie nacheinander für
/// dessen Teile, darf also nichts von einem Element zum nächsten mitnehmen.
#[cfg(feature = "parallel")]
pub(super) fn gather<R: Send>(count: usize, work: impl Fn(Range<usize>, &mut Vec<R>) + Sync) -> Vec<R> {
    use rayon::prelude::*;
    if count <= PIECE {
        let mut found = Vec::new();
        work(0..count, &mut found);
        return found;
    }
    let pieces: Vec<Vec<R>> = (0..count.div_ceil(PIECE))
        .into_par_iter()
        .map(|piece| {
            let mut found = Vec::new();
            work(piece * PIECE..((piece + 1) * PIECE).min(count), &mut found);
            found
        })
        .collect();
    let mut found = Vec::with_capacity(pieces.iter().map(Vec::len).sum());
    for piece in pieces {
        found.extend(piece);
    }
    found
}

/// Ohne das Feature `parallel`: alles in einem Stück.
#[cfg(not(feature = "parallel"))]
pub(super) fn gather<R>(count: usize, work: impl Fn(Range<usize>, &mut Vec<R>)) -> Vec<R> {
    let mut found = Vec::new();
    work(0..count, &mut found);
    found
}

#[cfg(test)]
mod tests {
    use super::gather;
    use alloc::vec::Vec;

    /// Die Funde kommen vollständig und in der Reihenfolge der Elemente
    /// zurück, gleich ob in einem Stück oder auf Kerne verteilt.
    #[test]
    fn keeps_order_and_finds_everything() {
        for count in [0, 1, 511, 512, 513, 10_000] {
            // Jedes dritte Element liefert zwei Funde, die übrigen keinen.
            let found = gather(count, |elements, found| {
                for element in elements {
                    if element % 3 == 0 {
                        found.push((element, 0));
                        found.push((element, 1));
                    }
                }
            });
            let expected: Vec<(usize, u8)> =
                (0..count).filter(|element| element % 3 == 0).flat_map(|element| [(element, 0), (element, 1)]).collect();
            assert_eq!(found, expected);
        }
    }
}
