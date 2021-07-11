use doc_comment::doc_comment;

/// Describes on which side of a line a a point lies.
///
/// Created by [DirectedEdgeHandle::side_query](handles/type.DirectedEdgeHandle.html#method.side_query)
#[derive(Debug, Clone, Copy)]
pub struct LineSideInfo {
    signed_side: f64,
}

impl PartialEq for LineSideInfo {
    fn eq(&self, other: &LineSideInfo) -> bool {
        if self.is_on_line() || other.is_on_line() {
            self.is_on_line() && other.is_on_line()
        } else {
            self.is_on_right_side() == other.is_on_right_side()
        }
    }
}

macro_rules! gen_edge_side_method {
    ($name:ident, $method:expr, $first_side:expr, $second_side:expr) => {
        gen_edge_side_method!(
            $name,
            $method,
            $first_side,
            $second_side,
            "This method returns `false` if the point lies exactly on the line."
        );
    };
    ($name:ident, $method:expr, $first_side:expr, $second_side:expr, $on_line:expr) => {
        doc_comment! {
            concat!("Returns `true` if a point lies on the ", $first_side, " side of a line.",
            "\n\n",
            "For left handed coordinate systems, this method returns if a point lies",
            " on the ", $second_side, " side of a line.",
            "\n\n",
            $on_line
        ),
            #[inline]
            pub fn $name(self) -> bool {
                $method(self.signed_side)
            }
        }
    };
}

impl LineSideInfo {
    #[inline]
    pub(crate) fn from_determinant(s: f64) -> LineSideInfo {
        LineSideInfo { signed_side: s }
    }

    gen_edge_side_method!(is_on_left_side, (|det| det > 0.0), "left", "right");

    gen_edge_side_method!(is_on_right_side, (|det| det < 0.0), "right", "left");

    gen_edge_side_method!(
        is_on_left_side_or_on_line,
        (|det| det >= 0.0),
        "left",
        "right",
        ""
    );

    gen_edge_side_method!(
        is_on_right_side_or_on_line,
        (|det| det >= 0.0),
        "right",
        "left",
        ""
    );

    /// Returns `true` if a point lies exactly on this line.
    #[inline]
    pub fn is_on_line(self) -> bool {
        self.signed_side.abs() == 0.0
    }

    /// Returns the opposite of this `LineSideInfo`.
    pub fn reversed(self) -> LineSideInfo {
        LineSideInfo {
            signed_side: -self.signed_side,
        }
    }
}
