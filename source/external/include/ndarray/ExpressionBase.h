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
#ifndef NDARRAY_ExpressionBase_h_INCLUDED
#define NDARRAY_ExpressionBase_h_INCLUDED

/** 
 *  @file ndarray/ExpressionBase.h
 *
 *  @brief Definitions for ExpressionBase.
 */

#include "ndarray/ExpressionTraits.h"
#include "ndarray/Vector.h"

namespace ndarray {

/**
 *  @class ExpressionBase
 *  @brief CRTP base class for all multidimensional expressions.
 *
 *  @ingroup MainGroup
 *
 *  ExpressionBase is a CRTP base class for both true array objects (subclasses of
 *  ArrayBase) and lazy array expressions, which are created by most arithmetic,
 *  bitwise, and logical operations on arrays.  These lazy expressions have most
 *  of the features of a true array to a const data type.
 *
 *  ExpressionBase also provides implementations for a few STL compatibility and
 *  convenience member functions.
 */
template <typename Derived>
class ExpressionBase {
public:
    /// @brief Data type of expression elements.
    typedef typename ExpressionTraits<Derived>::Element Element;
    /// @brief Number of dimensions (boost::mpl::int_).
    typedef typename ExpressionTraits<Derived>::ND ND;
    /// @brief Nested expression or element iterator.
    typedef typename ExpressionTraits<Derived>::Iterator Iterator;
    /// @brief Nested expression or element reference.
    typedef typename ExpressionTraits<Derived>::Reference Reference;
    /// @brief Nested expression or element value type.
    typedef typename ExpressionTraits<Derived>::Value Value;
    /// @brief Vector type for N-dimensional indices.
    typedef Vector<Size,ND::value> Index;
    /// @brief CRTP derived type.
    typedef Derived Self;

    /// @brief Return a single nested expression or element.
    Reference operator[](Size n) const { return getSelf().operator[](n); }

    /// @brief Return the first nested expression or element.
    Reference front() const { return this->operator[](0); }

    /// @brief Return the last nested expression or element.
    Reference back() const { return this->operator[](this->template getSize<0>()-1); }

    /// @brief Return an Iterator to the beginning of the expression.
    Iterator begin() const { return getSelf().begin(); }

    /// @brief Return an Iterator to one past the end of the expression.
    Iterator end() const { return getSelf().end(); }

    /// @brief Return the size of a specific dimension.
    template <int P> Size getSize() const { return getSelf().template getSize<P>(); }

    /// @brief Return a Vector of the sizes of all dimensions.
    Index getShape() const { return getSelf().getShape(); }

    /// @brief Return the total number of elements in the expression.
    Size getNumElements() const { return getSelf().getNumElements(); }

    /* ------------------------- STL Compatibility -------------------------- */

    typedef Value value_type;
    typedef Iterator iterator;
    typedef Iterator const_iterator;
    typedef Reference reference;
    typedef Reference const_reference;
    typedef Iterator pointer;
    typedef Offset difference_type;
    typedef Size size_type;

    /// @brief Return the size of the first dimension.
    size_type size() const { return this->template getSize<0>(); }

    /// @brief Return true if the first dimension has no elements.
    bool empty() const { return this->template getSize<0>() == 0; }

    /* ---------------------------------------------------------------------- */

protected:
    Self & getSelf() { return *static_cast<Self *>(this); }
    Self const & getSelf() const { return *static_cast<Self const *>(this); }
};

} // namespace ndarray

#endif // !NDARRAY_ExpressionBase_h_INCLUDED
