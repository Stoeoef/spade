use crate::{HasPosition, SpadeNum};

impl<S> HasPosition for mint::Point2<S>
where
    S: SpadeNum,
{
    type Scalar = S;

    fn position(&self) -> crate::Point2<Self::Scalar> {
        (*self).into()
    }
}

impl<S> From<mint::Point2<S>> for crate::Point2<S> {
    fn from(value: mint::Point2<S>) -> Self {
        crate::Point2::new(value.x, value.y)
    }
}

impl<S> From<crate::Point2<S>> for mint::Point2<S> {
    fn from(value: crate::Point2<S>) -> Self {
        mint::Point2 {
            x: value.x,
            y: value.y,
        }
    }
}
