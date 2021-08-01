use std::convert::TryInto;

use super::super::DCEL;
use super::public_handles::{InnerOuterMarker, PossiblyOuterTag};

#[cfg(feature = "serde")]
use serde_crate::{Deserialize, Serialize};

pub trait DelaunayElementType: Sized + Default {
    fn num_elements<V, DE, UE, F>(dcel: &DCEL<V, DE, UE, F>) -> usize;
}
/// Internal type definition that is only exposed for documentation purposes.
///
/// Rust will currently not generate documentation for type definitions depending
/// `pub(crate)` types, see [#32077](https://github.com/rust-lang/rust/issues/32077).
///
/// Do not use these types. Their removal will not be considered a breaking change.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(crate = "serde_crate")
)]
pub struct FixedHandleImpl<Type, InnerOuter: InnerOuterMarker> {
    index: u32,
    ty: Type,
    inner_outer: InnerOuter,
}

impl<Type: Default, InnerOuter: InnerOuterMarker> optional::Noned
    for FixedHandleImpl<Type, InnerOuter>
{
    fn is_none(&self) -> bool {
        optional::wrap(self.index).is_none()
    }

    fn get_none() -> Self {
        Self {
            index: optional::Optioned::<_>::none().unpack(),
            ty: Default::default(),
            inner_outer: Default::default(),
        }
    }
}

impl<Type: Default, InnerOuter: InnerOuterMarker> optional::OptEq
    for FixedHandleImpl<Type, InnerOuter>
{
    fn opt_eq(&self, other: &Self) -> bool {
        self.index().opt_eq(&other.index())
    }
}

impl<Type: Default, InnerOuter: InnerOuterMarker> optional::OptOrd
    for FixedHandleImpl<Type, InnerOuter>
{
    fn opt_cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.index().opt_cmp(&other.index())
    }
}

impl<Type, InnerOuter: InnerOuterMarker> std::fmt::Debug for FixedHandleImpl<Type, InnerOuter> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FixedHandle")
            .field("index", &self.index)
            .finish()
    }
}

pub const fn new_fixed_face_handle(index: usize) -> FixedHandleImpl<FaceTag, PossiblyOuterTag> {
    FixedHandleImpl {
        index: index as u32,
        ty: FaceTag,
        inner_outer: PossiblyOuterTag,
    }
}

impl<Type: Default, InnerOuter: InnerOuterMarker> FixedHandleImpl<Type, InnerOuter> {
    pub(crate) fn new(index: usize) -> Self {
        Self::new_internal(
            index
                .try_into()
                .expect("Index too big - at most 2^32 elements supported"),
        )
    }

    pub(crate) fn max() -> Self {
        Self::new_internal(u32::MAX)
    }

    pub(crate) fn index(&self) -> usize {
        self.index as usize
    }

    pub(crate) fn adjust_inner_outer<TargetInnerOuter: InnerOuterMarker>(
        &self,
    ) -> FixedHandleImpl<Type, TargetInnerOuter> {
        FixedHandleImpl::<_, _>::new_internal(self.index)
    }

    fn new_internal(index: u32) -> Self {
        Self {
            index,
            ty: Type::default(),
            inner_outer: InnerOuter::default(),
        }
    }
}

/// Internal type definition that is only exposed for documentation purposes.
///
/// Rust will currently not generate documentation for type definitions depending
/// `pub(crate)` types, see [#32077](https://github.com/rust-lang/rust/issues/32077).
///
/// Do not use these types. Their removal will not be considered a breaking change.
pub struct DynamicHandleImpl<'a, V, DE, UE, F, Type, InnerOuter: InnerOuterMarker> {
    pub(super) dcel: &'a DCEL<V, DE, UE, F>,
    pub(super) handle: FixedHandleImpl<Type, InnerOuter>,
}

impl<'a, V, DE, UE, F, Type: Default, InnerOuter: InnerOuterMarker>
    DynamicHandleImpl<'a, V, DE, UE, F, Type, InnerOuter>
{
    #[inline]
    pub(crate) fn new(
        dcel: &'a DCEL<V, DE, UE, F>,
        handle: FixedHandleImpl<Type, InnerOuter>,
    ) -> Self {
        Self { dcel, handle }
    }

    pub(in super::super) fn adjust_inner_outer<TargetInnerOuter: InnerOuterMarker>(
        &self,
    ) -> DynamicHandleImpl<'a, V, DE, UE, F, Type, TargetInnerOuter> {
        DynamicHandleImpl::new(self.dcel, self.handle.adjust_inner_outer())
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(crate = "serde_crate")
)]
pub struct VertexTag;
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(crate = "serde_crate")
)]
pub struct DirectedEdgeTag;
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(crate = "serde_crate")
)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
pub struct UndirectedEdgeTag;
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(crate = "serde_crate")
)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
pub struct FaceTag;
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
pub struct DirectedVoronoiEdgeTag;
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(crate = "serde_crate")
)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
pub struct UndirectedVoronoiEdgeTag;
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
pub struct VoronoiVertexTag;
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Default, Hash)]
pub struct VoronoiFaceTag;

impl DelaunayElementType for VertexTag {
    fn num_elements<V, DE, UE, F>(dcel: &DCEL<V, DE, UE, F>) -> usize {
        dcel.num_vertices()
    }
}

impl DelaunayElementType for DirectedEdgeTag {
    fn num_elements<V, DE, UE, F>(dcel: &DCEL<V, DE, UE, F>) -> usize {
        dcel.num_directed_edges()
    }
}

impl DelaunayElementType for UndirectedEdgeTag {
    fn num_elements<V, DE, UE, F>(dcel: &DCEL<V, DE, UE, F>) -> usize {
        dcel.num_undirected_edges()
    }
}

impl DelaunayElementType for FaceTag {
    fn num_elements<V, DE, UE, F>(dcel: &DCEL<V, DE, UE, F>) -> usize {
        dcel.num_faces()
    }
}

impl DelaunayElementType for VoronoiFaceTag {
    fn num_elements<V, DE, UE, F>(dcel: &DCEL<V, DE, UE, F>) -> usize {
        dcel.num_vertices()
    }
}

impl DelaunayElementType for DirectedVoronoiEdgeTag {
    fn num_elements<V, DE, UE, F>(dcel: &DCEL<V, DE, UE, F>) -> usize {
        dcel.num_directed_edges()
    }
}

impl DelaunayElementType for UndirectedVoronoiEdgeTag {
    fn num_elements<V, DE, UE, F>(dcel: &DCEL<V, DE, UE, F>) -> usize {
        dcel.num_undirected_edges()
    }
}

#[cfg(test)]
pub mod test {
    use optional::Optioned;

    use crate::handles::FixedDirectedEdgeHandle;

    #[test]
    fn optioned_handle_is_smaller_than_option_type() {
        let optioned_size = std::mem::size_of::<Optioned<FixedDirectedEdgeHandle>>();
        let option_size = std::mem::size_of::<Option<FixedDirectedEdgeHandle>>();

        assert!(optioned_size < option_size);
    }
}
