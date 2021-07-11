use std::marker::PhantomData;

use super::super::handle_defs::{DelaunayElementType, DynamicHandleImpl, FixedHandleImpl};
use super::super::InnerOuterMarker;

use crate::delaunay_core::DCEL;

pub struct FixedHandleIterator<Type, InnerOuter> {
    range: std::ops::Range<usize>,
    ty: PhantomData<Type>,
    inner_outer: PhantomData<InnerOuter>,
}

impl<Type: Default, InnerOuter: InnerOuterMarker> FixedHandleIterator<Type, InnerOuter> {
    pub(crate) fn new(to: usize) -> Self {
        Self {
            range: 0..to,
            ty: Default::default(),
            inner_outer: Default::default(),
        }
    }
}

impl<Type: Default, InnerOuter: InnerOuterMarker> Iterator
    for FixedHandleIterator<Type, InnerOuter>
{
    type Item = FixedHandleImpl<Type, InnerOuter>;

    fn next(&mut self) -> Option<Self::Item> {
        self.range.next().map(FixedHandleImpl::new)
    }
}

pub struct DynamicHandleIterator<'a, V, DE, UE, F, Type, InnerOuter> {
    fixed_iterator: FixedHandleIterator<Type, InnerOuter>,
    dcel: &'a DCEL<V, DE, UE, F>,
    inner_outer: std::marker::PhantomData<InnerOuter>,
}

impl<'a, V, DE, UE, F, Type, InnerOuter> DynamicHandleIterator<'a, V, DE, UE, F, Type, InnerOuter>
where
    Type: DelaunayElementType + Default,
    InnerOuter: InnerOuterMarker,
{
    pub(crate) fn new(dcel: &'a DCEL<V, DE, UE, F>) -> Self {
        DynamicHandleIterator {
            fixed_iterator: FixedHandleIterator::new(Type::num_elements(dcel)),
            dcel,
            inner_outer: std::marker::PhantomData::default(),
        }
    }
}

impl<'a, V, DE, UE, F, Type, InnerOuter> Iterator
    for DynamicHandleIterator<'a, V, DE, UE, F, Type, InnerOuter>
where
    Type: DelaunayElementType,
    InnerOuter: InnerOuterMarker,
{
    type Item = DynamicHandleImpl<'a, V, DE, UE, F, Type, InnerOuter>;

    fn next(&mut self) -> Option<Self::Item> {
        self.fixed_iterator
            .next()
            .map(|handle| DynamicHandleImpl::new(self.dcel, handle.adjust_inner_outer()))
    }
}
