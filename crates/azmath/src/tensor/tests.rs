//! Tests für Tensoren: Struktur, Fehlerfälle und Vergleiche mit Vektoren
//! und Matrizen fester Größe. Keine vorab ausgerechneten Erwartungswerte.

use super::{TensorView, TensorViewMut, contract, dot, outer};
use crate::linalg::{LinalgError, Matrix, Vector};
use crate::scalar::Real;

fn check<T: Real>() {
    let f = T::from_f64;
    let data: [T; 12] = core::array::from_fn(|i| f(i as f64) - f(4.5));

    // Sicht der Stufe 3 mit Form [2, 3, 2]: letzter Index läuft am schnellsten.
    assert!(TensorView::new(&data, [2, 3, 2]).is_ok_and(|t| {
        t.shape() == [2, 3, 2]
            && t.len() == 12
            && !t.is_empty()
            && t[[1, 2, 0]] == data[(3 + 2) * 2]
            && t.get([0, 1, 1]) == Some(&data[2 + 1])
            && t.get([2, 0, 0]).is_none()
            && t.get([0, 3, 0]).is_none()
            && t.get([0, 0, 2]).is_none()
            && t.data().len() == 12
    }));
    assert!(matches!(
        TensorView::new(&data, [2, 3, 3]),
        Err(LinalgError::DimensionMismatch)
    ));
    // Stufe 0: ein einzelner Wert.
    assert!(TensorView::new(&data[..1], []).is_ok_and(|t| t[[]] == data[0] && t.len() == 1));
    // Leerer Tensor.
    assert!(TensorView::<T, 2>::new(&[], [0, 5]).is_ok_and(|t| t.is_empty()));

    // Umformen behält die Daten, auch über Stufen hinweg.
    assert!(TensorView::new(&data, [2, 3, 2]).is_ok_and(|t| {
        t.reshape([4, 3]).is_ok_and(|m| m[[1, 1]] == data[3 + 1] && m.shape() == [4, 3])
            && t.reshape([12]).is_ok_and(|v| v[[7]] == data[7])
            && t.reshape([5, 2]).is_err()
    }));

    // Schreibende Sicht.
    let mut buffer = data;
    assert!(matches!(
        TensorViewMut::new(&mut buffer, [5, 5]),
        Err(LinalgError::DimensionMismatch)
    ));
    assert!(TensorViewMut::new(&mut buffer, [3, 4]).is_ok_and(|mut t| {
        t[[0, 0]] = T::ONE;
        if let Some(slot) = t.get_mut([2, 3]) {
            *slot = T::ZERO;
        }
        t.shape() == [3, 4]
            && t.len() == 12
            && !t.is_empty()
            && t[[0, 0]] == T::ONE
            && t.get([2, 3]) == Some(&T::ZERO)
            && t.get([3, 0]).is_none()
            && t.get_mut([0, 4]).is_none()
            && t.as_view()[[0, 1]] == data[1]
    }));
    assert_eq!((buffer[0], buffer[11]), (T::ONE, T::ZERO));

    // Elementweise Operationen: x + x = 2x.
    let mut doubled = data;
    let mut summed = data;
    assert!(TensorViewMut::new(&mut doubled, [3, 4]).is_ok_and(|mut t| {
        t.map_in_place(|x| x * f(2.0));
        true
    }));
    assert!(TensorView::new(&data, [3, 4]).is_ok_and(|other| {
        TensorViewMut::new(&mut summed, [3, 4]).is_ok_and(|mut t| {
            t.zip_in_place(other, |x, y| x + y).is_ok()
                && other.reshape([4, 3]).is_ok_and(|wrong| t.zip_in_place(wrong, |x, _| x).is_err())
        })
    }));
    assert_eq!(doubled, summed);

    // Vergleich mit Vektoren und Matrizen fester Größe.
    let u = Vector([f(1.0), f(-2.0), f(0.5)]);
    let w = Vector([f(0.3), f(0.7), f(-1.1)]);
    let m = Matrix([[f(1.0), f(0.5), f(-0.3)], [f(0.2), f(1.5), f(0.4)]]);
    let n = Matrix([[f(4.0), f(1.0)], [f(-1.0), f(5.0)], [f(0.3), f(-0.7)]]);
    let m_flat = [m[(0, 0)], m[(0, 1)], m[(0, 2)], m[(1, 0)], m[(1, 1)], m[(1, 2)]];
    let n_flat = [n[(0, 0)], n[(0, 1)], n[(1, 0)], n[(1, 1)], n[(2, 0)], n[(2, 1)]];
    let (tu, tw) = (TensorView::new(&u.0, [3]), TensorView::new(&w.0, [3]));
    let (tm, tn) = (TensorView::new(&m_flat, [2, 3]), TensorView::new(&n_flat, [3, 2]));

    let ok = match (tu, tw, tm, tn) {
        (Ok(tu), Ok(tw), Ok(tm), Ok(tn)) => {
            // Äußeres Produkt zweier Vektoren.
            let mut product = [T::ZERO; 9];
            let expected = u.outer(w);
            let outer_ok = outer(tu, tw, &mut product).is_ok()
                && TensorView::new(&product, [3, 3]).is_ok_and(|t| {
                    (0..3).all(|i| (0..3).all(|j| t[[i, j]] == expected[(i, j)]))
                });

            // Skalarprodukt, Matrix · Vektor, Matrix · Matrix.
            let mut scalar = [T::ZERO; 1];
            let scalar_ok = dot(tu, tw, &mut scalar).is_ok() && scalar[0] == u.dot(w);
            let mut image = [T::ZERO; 2];
            let vector_ok = dot(tm, tu, &mut image).is_ok() && image == (m * u).0;
            let mut square = [T::ZERO; 4];
            let mn = m * n;
            let matrix_ok = dot(tm, tn, &mut square).is_ok()
                && square == [mn[(0, 0)], mn[(0, 1)], mn[(1, 0)], mn[(1, 1)]];

            // Kontraktion: Spur, und Spur des äußeren Produkts = Skalarprodukt.
            let mut trace = [T::ZERO; 1];
            let trace_ok = TensorView::new(&square, [2, 2])
                .is_ok_and(|t| contract(t, 0, 1, &mut trace).is_ok())
                && trace[0] == mn.trace();
            let mut outer_trace = [T::ZERO; 1];
            let identity_ok = TensorView::new(&product, [3, 3])
                .is_ok_and(|t| contract(t, 1, 0, &mut outer_trace).is_ok())
                && crate::testutil::close(outer_trace[0], u.dot(w), 8.0);

            // Kontraktion bei Stufe 3 lässt den übrigen Index stehen.
            let mut rest = [T::ZERO; 3];
            let cube: [T; 12] = data;
            let rank3_ok = TensorView::new(&cube, [2, 3, 2]).is_ok_and(|t| {
                contract(t, 0, 2, &mut rest).is_ok()
                    && (0..3).all(|j| rest[j] == t[[0, j, 0]] + t[[1, j, 1]])
            });

            // Fehlerfälle.
            let mut wrong = [T::ZERO; 5];
            let errors_ok = outer(tu, tw, &mut wrong).is_err()
                && dot(tu, tw, &mut wrong).is_err()
                && dot(tm, tm, &mut square).is_err()
                && contract(tm, 0, 1, &mut scalar).is_err()
                && contract(tn, 0, 0, &mut scalar).is_err()
                && contract(tn, 0, 2, &mut scalar).is_err()
                && TensorView::new(&square, [2, 2])
                    .is_ok_and(|t| contract(t, 0, 1, &mut wrong).is_err());

            outer_ok
                && scalar_ok
                && vector_ok
                && matrix_ok
                && trace_ok
                && identity_ok
                && rank3_ok
                && errors_ok
        }
        _ => false,
    };
    assert!(ok);
}

#[test]
fn tensor_f32() {
    check::<f32>();
}

#[test]
fn tensor_f64() {
    check::<f64>();
}

#[cfg(feature = "alloc")]
mod owned {
    use super::*;
    use crate::tensor::DTensor;
    use alloc::vec;
    use alloc::vec::Vec;

    fn check<T: Real>() {
        let f = T::from_f64;
        let a: Vec<T> = (0..6).map(|i| f(f64::from(i)) - f(2.5)).collect();
        let b: Vec<T> = (0..6).map(|i| f(0.5) * f(f64::from(i)) + f(1.0)).collect();
        let ta = DTensor::from_vec(a.clone(), [2, 3]);
        let tb = DTensor::from_vec(b.clone(), [3, 2]);

        let ok = match (ta, tb) {
            (Ok(ta), Ok(tb)) => {
                // Dieselben Ergebnisse wie die Operationen auf Sichten.
                let mut flat = vec![T::ZERO; 36];
                let outer_ok = outer(ta.view(), tb.view(), &mut flat).is_ok()
                    && DTensor::<T, 4>::outer(ta.view(), tb.view())
                        .is_ok_and(|t| t.shape() == [2, 3, 3, 2] && t.into_vec() == flat);
                let mut product = vec![T::ZERO; 4];
                let dot_ok = dot(ta.view(), tb.view(), &mut product).is_ok()
                    && DTensor::<T, 2>::dot(ta.view(), tb.view())
                        .is_ok_and(|t| t.shape() == [2, 2] && t.clone().into_vec() == product);
                let mut trace = vec![T::ZERO; 1];
                let contract_ok = DTensor::<T, 2>::dot(ta.view(), tb.view()).is_ok_and(|square| {
                    contract(square.view(), 0, 1, &mut trace).is_ok()
                        && DTensor::<T, 0>::contract(square.view(), 0, 1)
                            .is_ok_and(|t| t.shape().is_empty() && t.into_vec() == trace)
                });

                // Falsche Zielstufe.
                let rank_ok = DTensor::<T, 3>::outer(ta.view(), tb.view()).is_err()
                    && DTensor::<T, 3>::dot(ta.view(), tb.view()).is_err()
                    && DTensor::<T, 1>::contract(ta.view(), 0, 1).is_err();

                // Sichten auf den eigenen Speicher.
                let mut copy = ta.clone();
                copy.view_mut()[[0, 0]] = T::ONE;
                let view_ok = copy.view()[[0, 0]] == T::ONE
                    && copy.view()[[1, 2]] == ta.view()[[1, 2]]
                    && copy != ta
                    && ta.shape() == [2, 3];

                outer_ok && dot_ok && contract_ok && rank_ok && view_ok
            }
            _ => false,
        };
        assert!(ok);

        assert!(DTensor::<T, 2>::zeros([2, 3]).is_ok_and(|t| t.into_vec() == vec![T::ZERO; 6]));
        assert!(DTensor::<T, 2>::zeros([usize::MAX, 2]).is_err());
        assert!(DTensor::from_vec(a, [4, 2]).is_err());
        assert!(DTensor::from_vec(b, [6]).is_ok());
    }

    #[test]
    fn owned_f32() {
        check::<f32>();
    }

    #[test]
    fn owned_f64() {
        check::<f64>();
    }
}
