// -*- c++ -*-
/*
 * Copyright (c) 2010-2012, Jim Bosch
 * All rights reserved.
 *
 * ndarray is distributed under a simple BSD-like license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/ndarray/ndarray
 */
#ifndef NDARRAY_DETAIL_BinaryOp_h_INCLUDED
#define NDARRAY_DETAIL_BinaryOp_h_INCLUDED

/** 
 *  @file ndarray/detail/BinaryOp.h
 *
 *  @brief Lazy binary expression templates.
 */

#include "ndarray/ExpressionBase.h"
#include "ndarray/vectorize.h"
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

namespace ndarray {
namespace detail {

/**
 *  @internal @brief An iterator for binary expression templates.
 *
 *  @ingroup ndarrayInternalGroup
 *
 *  Acts as a combination "zip" and "transform" iterator.
 */
template <typename Operand1, typename Operand2, typename BinaryFunction>
class BinaryOpIterator : public boost::iterator_adaptor<
    BinaryOpIterator<Operand1,Operand2,BinaryFunction>,
    boost::zip_iterator<
        boost::tuple<
            typename ExpressionTraits<Operand1>::Iterator,
            typename ExpressionTraits<Operand2>::Iterator
            >
        >,
    typename ExpressionTraits< BinaryOpExpression<Operand1,Operand2,BinaryFunction> >::Value,
    boost::use_default,
    typename ExpressionTraits< BinaryOpExpression<Operand1,Operand2,BinaryFunction> >::Reference
    > {
    typedef BinaryOpExpression<Operand1,Operand2,BinaryFunction> Operation;
public:
    typedef typename ExpressionTraits<Operand1>::Iterator BaseIterator1;
    typedef typename ExpressionTraits<Operand2>::Iterator BaseIterator2;
    typedef typename ExpressionTraits<Operation>::Value Value;
    typedef typename ExpressionTraits<Operation>::Reference Reference;

    BinaryOpIterator() : BinaryOpIterator::iterator_adaptor_(), _functor() {}

    BinaryOpIterator(
        BaseIterator1 const & baseIter1, 
        BaseIterator2 const & baseIter2, 
        BinaryFunction const & functor
    ) :
        BinaryOpIterator::iterator_adaptor_(boost::make_tuple(baseIter1,baseIter2)),
        _functor(functor) {}

    BinaryOpIterator(BinaryOpIterator const & other) :
        BinaryOpIterator::iterator_adaptor_(other), _functor(other._functor) {}

private:
    friend class boost::iterator_core_access;

    Reference dereference() const {
        return vectorize(
            _functor, 
            this->base_reference()->template get<0>(),
            this->base_reference()->template get<1>()
        );
    }
    
    BinaryFunction _functor;
};

/**
 *  @internal @brief A binary expression template.
 *
 *  @ingroup ndarrayInternalGroup
 *
 *  Represents the lazy evaluation of a binary expression.
 */
template <typename Operand1, typename Operand2, typename BinaryFunction, int N>
class BinaryOpExpression : public ExpressionBase< BinaryOpExpression<Operand1,Operand2,BinaryFunction,N> > {
    typedef BinaryOpExpression<Operand1,Operand2,BinaryFunction,N> Self;
public:
    typedef typename ExpressionTraits<Self>::Element Element;
    typedef typename ExpressionTraits<Self>::ND ND;
    typedef typename ExpressionTraits<Self>::Iterator Iterator;
    typedef typename ExpressionTraits<Self>::Value Value;
    typedef typename ExpressionTraits<Self>::Reference Reference;
    typedef Vector<Size,N> Index;
    
    BinaryOpExpression(
        Operand1 const & operand1,
        Operand2 const & operand2,
        BinaryFunction const & functor
    ) :
        _operand1(operand1), _operand2(operand2), _functor(functor) {
        NDARRAY_ASSERT(_operand1.getShape() == _operand2.getShape());
    }

    Reference operator[](Size n) const {
        return Reference(_operand1[n],_operand2[n],_functor);
    }

    Iterator begin() const {
        return Iterator(_operand1.begin(),_operand2.begin(),_functor);
    }

    Iterator end() const {
        return Iterator(_operand1.end(),_operand2.end(),_functor);
    }

    template <int P> Size getSize() const {
        return _operand1.template getSize<P>();
    }

    Index getShape() const {
        return _operand1.getShape();
    }

    Operand1 _operand1;
    Operand2 _operand2;
    BinaryFunction _functor;
};

} // namespace detail
} // namespace ndarray

#endif // !NDARRAY_DETAIL_BinaryOp_h_INCLUDED
