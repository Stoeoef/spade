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

use traits::{SpadeNum, HasPosition, TwoDimensional};
use delaunay::DelaunayTriangulation;
use primitives::{SimpleEdge, EdgeSideInfo};
use bigvec::{BigVec2, AdaptiveInt};
use exactpred::{orient2d, incircle};

/// Determines how a delaunay triangulation performs its basic geometry computation.
/// 
/// Every delaunay triangulation is based on two basic geometry operations: orientation tests
/// (on which side of a line lies a point?) and in-circle tests (is a point contained in the
/// circumference of a triangle?). These questions can be answered *approximately* or *precise*,
/// their calculation can take overflow issues into account (if handling integer coordinates).
///
/// Since each application has different needs, a `DelaunayKernel` will define how these geometric
/// queries are calculated for each triangulation. It is recommended to use one of the predefined
/// kernels that fits your needs.
pub trait DelaunayKernel<D: SpadeNum>: ::std::marker::Sized {
    /// Returns true if pd is contained in the circumference of the triangle spanned by pa, pb, pc.
    /// pa, pb, pc have to be ordered clockwise.
    fn contained_in_circumference<V: TwoDimensional<Scalar=D>>(pa: &V, pb: &V, pc: &V, pd: &V) -> bool {
        let adx = pa[0].clone() - pd[0].clone();
        let ady = pa[1].clone() - pd[1].clone();
        let bdx = pb[0].clone() - pd[0].clone();
        let bdy = pb[1].clone() - pd[1].clone();
        let cdx = pc[0].clone() - pd[0].clone();
        let cdy = pc[1].clone() - pd[1].clone();

        let abdet = adx.clone() * bdy.clone() - bdx.clone() * ady.clone();
        let bcdet = bdx.clone() * cdy.clone() - cdx.clone() * bdy.clone();
        let cadet = cdx.clone() * ady.clone() - adx.clone() * cdy.clone();
        let alift = adx.clone() * adx.clone() + ady.clone() * ady.clone();
        let blift = bdx.clone() * bdx.clone() + bdy.clone() * bdy.clone();
        let clift = cdx.clone() * cdx.clone() + cdy.clone() * cdy.clone();

        let det = alift.clone() * bcdet.clone() 
            + blift.clone() * cadet.clone()
            + clift.clone() * abdet.clone();
        det.is_negative()
    }

    /// Returns an `EdgeSideInfo` yielding on which side of a line a point lies.
    fn side_query<V: TwoDimensional<Scalar=D>>(edge: &SimpleEdge<V>, position: &V) -> EdgeSideInfo<D> {
        edge.side_query(position)
    }

    /// Another formulation of `side_query`, will return `true` if `v0`, `v1` and `v2` are ordered
    /// counterclockwise.
    fn is_ordered_ccw<V: TwoDimensional<Scalar=D>>(v0: &V, v1: &V, v2: &V) -> bool {
        let edge = SimpleEdge::new(v0.clone(), v1.clone());
        Self::side_query(&edge, v2).is_on_left_side_or_on_line()
    }

    /// Creates a new delaunay triangulation using this kernel.
    fn new_triangulation<T: HasPosition>() -> DelaunayTriangulation<T, Self>
                         where T::Vector: TwoDimensional<Scalar=D>
    {
        DelaunayTriangulation::new()
    }

    /// Returns if a point lies on the infinite edge going through
    /// `edge.from` and `edge.to`.
    fn point_on_edge<V: TwoDimensional<Scalar=D>>(
        edge: &SimpleEdge<V>, position: &V) -> bool {
        Self::side_query(edge, position).is_on_line()
            && edge.is_projection_on_edge(position)
    }
}

/// Offers fast, possibly approximative and possibly inaccurate geometric calculations.
///
/// Use this kernel if you are working with small integral coordinates (e.g. `Vector2<i64>`
/// in the range &#177;100000) or if you are willing to risk crashes when using `Vector2<f64>`
/// as coordinate. Offers best performance. Run this in debug mode to be notified early if
/// something goes amiss (e.g. over / underflow, contradictory geometric calculations) since
/// `DelaunayTriangulation` uses a few `debug_assert`s to keep its internal structure healthy.
///
/// If you run into said problems, consider using `FloatKernel` or `AdaptiveIntKernel`.
pub struct TrivialKernel { }

impl <N: SpadeNum> DelaunayKernel<N> for TrivialKernel { }

/// Delaunay kernel for integral coordinates with a larger value range.
///
/// Integer calculations do not suffer from precision loss (since no divisions have to be used),
/// yet they are prone to over- and underflow errors.
/// This kernel will use `AdaptiveInt`s for its internal computation. These integers will
/// transparently require more bits if they encounter an under- or overflow.
pub struct AdaptiveIntKernel { }

impl DelaunayKernel<i64> for AdaptiveIntKernel {
    fn contained_in_circumference<V: TwoDimensional<Scalar=i64>>(v1: &V, v2: &V, v3: &V, p: &V) -> bool {
        let to_bigvec = |v: &V| BigVec2::new(
            AdaptiveInt::from_i64(&v[0]), AdaptiveInt::from_i64(&v[1]));
        // Cast input to adaptive ints to prevent overflows
        let v1: BigVec2<_> = to_bigvec(v1);
        let v2: BigVec2<_> = to_bigvec(v2);
        let v3: BigVec2<_> = to_bigvec(v3);
        let p: BigVec2<_> = to_bigvec(p);

        TrivialKernel::contained_in_circumference(&v1, &v2, &v3, &p)
    }
}

/// Offers a fast, precise kernel working with `f64` coordinates.
///
/// Performing a delaunay triangulation is often a tradeoff between accuracy and speed:
/// a triangulation working on native `f64`-operations will fail in rare cases due to rounding
/// errors, while switching to precise floats (like `BigRationals` from the `num` crate) reduces
/// performance by several orders of magnitude.
/// This kernel works with adaptive precision: if a calculation is inaccurate,
/// it will increase its precision until the result is accurate enough. Since most calculations
/// are accurate enough in their simplest form, only the overhead of checking the precision is
///  usually encountered. For more information, refer to 
/// [this link](https://www.cs.cmu.edu/~quake/robust.html) describing the techniques behind
/// the adaptive approach.
pub struct FloatKernel { }

impl DelaunayKernel<f64> for FloatKernel {
    fn contained_in_circumference<V: TwoDimensional<Scalar=f64>>(v1: &V, v2: &V, v3: &V, p: &V) -> bool {
        incircle(v1, v2, v3, p) < 0.0
    }

    fn side_query<V: TwoDimensional<Scalar=f64>>(edge: &SimpleEdge<V>, position: &V) -> EdgeSideInfo<f64> {
        let det = orient2d(&edge.from, &edge.to, &position);
        EdgeSideInfo::from_determinant(det)
    }
}

#[cfg(test)]
mod test {
    use super::{TrivialKernel, DelaunayKernel};
    use nalgebra as na;

    #[test]
    fn test_contained_in_circumference() {
        let (a1, a2, a3) = (1f32, 2f32, 3f32);
        let offset = na::Vector2::new(0.5, 0.7);
        let v1 = na::Vector2::new(a1.sin(), a1.cos()) * 2. + offset;
        let v2 = na::Vector2::new(a2.sin(), a2.cos()) * 2. + offset;
        let v3 = na::Vector2::new(a3.sin(), a3.cos()) * 2. + offset;
        assert!(TrivialKernel::contained_in_circumference(&v1, &v2, &v3, &offset));
        let shrunk = (v1 - offset) * 0.9 + offset;
        assert!(TrivialKernel::contained_in_circumference(&v1, &v2, &v3, &shrunk));
        let expanded = (v1 - offset) * 1.1 + offset;
        assert!(!TrivialKernel::contained_in_circumference(&v1, &v2, &v3, &expanded));
        assert!(!TrivialKernel::contained_in_circumference(
            &v1, &v2, &v3, &(na::Vector2::new(2.0, 2.0) + offset)));
        assert!(TrivialKernel::contained_in_circumference(
            &na::Vector2::new(0f32, -1f32), &na::Vector2::new(0f32, 0f32),
            &na::Vector2::new(-1f32, 0f32), &na::Vector2::new(0f32, 1f32)));
    }
}
