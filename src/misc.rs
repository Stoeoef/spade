use cgmath::{BaseFloat, BaseNum, Vector2};

// A call to l.min(r) does not seem to be inlined, thus we define it ourselves
// This does improve performance significantly, especially for larger node sizes
#[inline]
pub fn fmin<'a, S: BaseFloat>(l: &'a S, r: &'a S) -> &'a S {
    if l < r {
        l
    } else {
        r
    }
}

#[inline]
pub fn fmax<'a, S: BaseFloat>(l: &'a S, r: &'a S) -> &'a S {
    if l > r {
        l
    } else {
        r
    }
}

#[inline]
pub fn clamp<S: BaseNum>(lower: S, upper: S, value: S) -> S {
    upper.partial_min(lower.partial_max(value))
}

pub fn length2<S: BaseNum>(vec: &Vector2<S>) -> S {
    vec.x * vec.x + vec.y * vec.y
}
