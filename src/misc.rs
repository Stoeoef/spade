// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::traits::{SpadeNum};

// A call to l.min(r) does not seem to be inlined, thus we define it ourselves
// This does improve performance significantly.
#[inline]
pub fn min_inline<S: SpadeNum>(a: S, b: S) -> S {
    if a < b {
        a
    } else {
        b
    }
}

#[inline]
pub fn max_inline<S: SpadeNum>(a: S, b: S) -> S {
    if a > b {
        a
    } else {
        b
    }
}
