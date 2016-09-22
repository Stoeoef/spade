// Copyright 2016 The Spade Developers. For a full listing of the authors,
// refer to the Cargo.toml file at the top-level directory of this distribution.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

use traits::{SpadeNum};

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
