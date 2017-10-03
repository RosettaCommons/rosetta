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
#ifndef NDARRAY_DETAIL_UnaryOp_h_INCLUDED
#define NDARRAY_DETAIL_UnaryOp_h_INCLUDED

/** 
 *  @file ndarray/detail/UnaryOp.h
 *
 *  @brief Lazy unary expression templates.
 */

#include "ndarray/ExpressionBase.h"
#include "ndarray/vectorize.h"
#include <boost/iterator/iterator_adaptor.hpp>

namespace ndarray {
namespace detail {

/**
 *  @internal @brief An iterator for unary expression templates.
 *
 *  @ingroup ndarrayInternalGroup
 *
 *  Acts as a "transform" iterator.
 */
template <typename Operand, typename UnaryFunction>
class UnaryOpIterator : public boost::iterator_adaptor<
    UnaryOpIterator<Operand,UnaryFunction>,
    typename ExpressionTraits<Operand>::Iterator,
    typename ExpressionTraits< UnaryOpExpression<Operand,UnaryFunction> >::Value,
    boost::use_default,
    typename ExpressionTraits< UnaryOpExpression<Operand,UnaryFunction> >::Reference
    > {
    typedef UnaryOpExpression<Operand,UnaryFunction> Operation;
public:
    typedef typename ExpressionTraits<Operand>::Iterator BaseIterator;
    typedef typename ExpressionTraits<Operation>::Value Value;
    typedef typename ExpressionTraits<Operation>::Reference Reference;

    UnaryOpIterator() : UnaryOpIterator::iterator_adaptor_(), _functor() {}

    UnaryOpIterator(BaseIterator const & baseIter, UnaryFunction const & functor) :
        UnaryOpIterator::iterator_adaptor_(baseIter), _functor(functor) {}

    UnaryOpIterator(UnaryOpIterator const & other) :
        UnaryOpIterator::iterator_adaptor_(other), _functor(other._functor) {}

private:
    friend class boost::iterator_core_access;

    Reference dereference() const {
        return vectorize(_functor,*this->base_reference());
    }
    
    UnaryFunction _functor;
};

/**
 *  @internal @brief A unary expression template.
 *
 *  @ingroup ndarrayInternalGroup
 *
 *  Represents the lazy evaluation of a unary expression.
 */
template <typename Operand, typename UnaryFunction, int N>
class UnaryOpExpression : public ExpressionBase< UnaryOpExpression<Operand,UnaryFunction,N> > {
    typedef UnaryOpExpression<Operand,UnaryFunction,N> Self;
public:
    typedef typename ExpressionTraits<Self>::Element Element;
    typedef typename ExpressionTraits<Self>::ND ND;
    typedef typename ExpressionTraits<Self>::Iterator Iterator;
    typedef typename ExpressionTraits<Self>::Value Value;
    typedef typename ExpressionTraits<Self>::Reference Reference;
    typedef Vector<Size,N> Index;
    
    UnaryOpExpression(Operand const & operand, UnaryFunction const & functor) :
        _operand(operand), _functor(functor) {}

    Reference operator[](Size n) const {
        return Reference(_operand[n],_functor);
    }

    Iterator begin() const {
        return Iterator(_operand.begin(),_functor);
    }

    Iterator end() const {
        return Iterator(_operand.end(),_functor);
    }

    template <int P> Size getSize() const {
        return _operand.template getSize<P>();
    }

    Index getShape() const {
        return _operand.getShape();
    }

    Operand _operand;
    UnaryFunction _functor;
};

} // namespace detail
} // namespace ndarray

#endif // !NDARRAY_DETAIL_UnaryOp_h_INCLUDED
