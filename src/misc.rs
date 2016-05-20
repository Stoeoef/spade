use cgmath::BaseFloat;

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
pub fn fclamp<'a, S: BaseFloat>(lower: &'a S, upper: &'a S, value: &'a S) -> &'a S {
    fmin(upper, fmax(lower, value))
}
