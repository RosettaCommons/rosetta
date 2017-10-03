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
#ifndef NDARRAY_DETAIL_ArrayAccess_h_INCLUDED
#define NDARRAY_DETAIL_ArrayAccess_h_INCLUDED

/** 
 *  @file ndarray/detail/ArrayAccess.h
 *
 *  @brief Definitions for ArrayAccess
 */

#include "ndarray/ExpressionTraits.h"

namespace ndarray {
namespace detail {

template <typename Array_>
class ArrayAccess {
public:
    typedef typename ExpressionTraits< Array_ >::Element Element;
    typedef typename ExpressionTraits< Array_ >::Core Core;
    typedef typename ExpressionTraits< Array_ >::CorePtr CorePtr;

    static CorePtr const & getCore(Array_ const & array) {
        return array._core;
    }

    static Array_ construct(Element * data, CorePtr const & core) {
        return Array_(data, core);
    }

};

} // namespace detail
} // namespace ndarray

#endif // !NDARRAY_DETAIL_ArrayAccess_h_INCLUDED
