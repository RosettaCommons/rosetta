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
#ifndef NDARRAY_vectorize_h_INCLUDED
#define NDARRAY_vectorize_h_INCLUDED

/** 
 *  \file ndarray/vectorize.h @brief Code to apply arbitrary scalar functors to arrays.
 */

#include "ndarray_fwd.h"
#include "ndarray/detail/UnaryOp.h"

#include <boost/mpl/and.hpp>
#include <boost/utility/enable_if.hpp>

namespace ndarray {
namespace result_of {

template <typename T1, typename T2, typename T3=void>
struct vectorize {
    typedef T1 BinaryFunction;
    typedef T2 Argument1;
    typedef T3 Argument2;

    typename boost::mpl::if_<
        boost::mpl::and_<
            typename ExpressionTraits<Argument1>::IsScalar,
            typename ExpressionTraits<Argument2>::IsScalar
        >,
        typename BinaryFunction::result_type,
        detail::BinaryOpExpression<Argument1,Argument2,BinaryFunction>
    >::type type;
    
};

template <typename T1, typename T2>
struct vectorize<T1,T2,void> {
    typedef T1 UnaryFunction;
    typedef T2 Argument;

    typedef typename boost::mpl::if_<
        typename ExpressionTraits<Argument>::IsScalar,
        typename UnaryFunction::result_type,
        detail::UnaryOpExpression<Argument,UnaryFunction>
    >::type type;
};

} // namespace result_of

/// @addtogroup ndarrayMainGroup
/// @{

/** 
 *  @brief Apply a non-mutating unary function object to a scalar.
 *
 *  This overload exists to allow recursive usage of the Array-argument vectorize functions.
 */
template <typename Scalar, typename UnaryFunction>
#ifndef DOXYGEN
typename boost::enable_if<typename ExpressionTraits<Scalar>::IsScalar,
                          typename UnaryFunction::result_type>::type
#else
typename UnaryFunction::result_type
#endif
vectorize(
    UnaryFunction const & functor,
    Scalar const & scalar
) {
    return functor(scalar);
}

/** 
 *  @brief Apply a non-mutating unary function object to each element of a multidimensional Expression.
 *
 *  Evaluation is lazy.
 */
template <typename Derived, typename UnaryFunction>
detail::UnaryOpExpression<Derived,UnaryFunction>
vectorize(
    UnaryFunction const & functor,
    ExpressionBase<Derived> const & operand
) {
    return detail::UnaryOpExpression<Derived,UnaryFunction>(
        static_cast<Derived const &>(operand),
        functor
    );
}

/** 
 *  @brief Apply a non-mutating binary function object to a pair of scalars.
 *
 *  This overload exists to allow recursive usage of the Array-argument vectorize functions.
 */
template <typename Scalar1, typename Scalar2, typename BinaryFunction>
#ifndef DOXYGEN
typename boost::enable_if_c<
    (ExpressionTraits<Scalar1>::IsScalar::value
     && ExpressionTraits<Scalar2>::IsScalar::value),
    typename BinaryFunction::result_type
    >::type
#else
typename BinaryFunction::result_type
#endif
vectorize(
    BinaryFunction const & functor,
    Scalar1 const & scalar1,
    Scalar2 const & scalar2
) {
    return functor(scalar1,scalar2);
}

/** 
 *  @brief Apply a non-mutating binary function object pairwise to 
 *  the elements of two multidimensional Expressions.
 *
 *  Evaluation is lazy.
 */
template <typename Derived1, typename Derived2, typename BinaryFunction>
detail::BinaryOpExpression<Derived1,Derived2,BinaryFunction>
vectorize(
    BinaryFunction const & functor,
    ExpressionBase<Derived1> const & operand1,
    ExpressionBase<Derived2> const & operand2
) {
    return detail::BinaryOpExpression<Derived1,Derived2,BinaryFunction>(
        static_cast<Derived1 const &>(operand1),
        static_cast<Derived2 const &>(operand2),
        functor
    );
}

/// @}

} // namespace ndarray

#endif // !NDARRAY_vectorize_h_INCLUDED
