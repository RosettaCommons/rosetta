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
#ifndef NDARRAY_ArrayTraits_h_INCLUDED
#define NDARRAY_ArrayTraits_h_INCLUDED

/** 
 *  @file ndarray/ArrayTraits.h
 *
 *  @brief Traits for Array.
 */

#include "ndarray_fwd.h"
#include "ndarray/ExpressionTraits.h"
#include "ndarray/detail/Core.h"
#include <boost/mpl/int.hpp>
#include <boost/mpl/bool.hpp>

namespace ndarray {
namespace detail {

template <int N, typename T2, int C2, typename T1, int C1>
struct Convertible : public boost::mpl::bool_<
    (((C2>=C1 && C1>=0) || (C2<=C1 && C1<=0) || (N == 1 && C2 == -C1)) 
     && boost::is_convertible<T2*,T1*>::value)
> {};

} // namespace detail


/**
 *  @brief Dimension-specialized traits shared by Array and ArrayRef.
 *
 *  @ingroup MainGroup
 */
template <typename T, int N, int C>
struct ArrayTraits {
    typedef T Element;
    typedef boost::mpl::int_<N> ND;
    typedef boost::mpl::int_<C> RMC;
    typedef detail::NestedIterator<T,N,C> Iterator;
    typedef ArrayRef<T,N-1,(N==C)?(N-1):((C>0)?C:0)> Reference;
    typedef Array<T,N-1,(N==C)?(N-1):((C>0)?C:0)> Value;
    typedef detail::Core<N> Core;
    typedef typename Core::ConstPtr CorePtr;

    static Reference makeReference(Element * data, CorePtr const & core) {
        return Reference(data, core);
    }
    static Iterator makeIterator(Element * data, CorePtr const & core, Offset stride) {
        return Iterator(Reference(data, core), stride);
    }
    static void fill(Iterator iter, Iterator const & end, Element value) {
        // We can't use std::fill here because NestedIterator is not formally an STL ForwardIterator;
        // it has random access traversal, but it does not dereference to an addressable type (see
        // http://www.boost.org/doc/libs/1_55_0/libs/iterator/doc/new-iter-concepts.html#motivation)
        // Most C++ standard libraries have a fill implementation that will accept NestedIterator
        // anyway, but Clang's libc++ is more strictly compliant and does not.
        for (; iter != end; ++iter) {
            *iter = value;
        }
    }
};

template <typename T>
struct ArrayTraits<T,1,0> {
    typedef T Element;
    typedef boost::mpl::int_<1> ND;
    typedef boost::mpl::int_<0> RMC;
    typedef detail::StridedIterator<Element> Iterator;
    typedef Element & Reference;
    typedef Element Value;
    typedef detail::Core<1> Core;
    typedef typename Core::ConstPtr CorePtr;

    static Reference makeReference(Element * data, CorePtr const & core) {
        return *data;
    }
    static Iterator makeIterator(Element * data, CorePtr const & core, Offset stride) {
        return Iterator(data, stride);
    }
    static void fill(Iterator iter, Iterator const & end, Element value) {
        std::fill(iter, end, value);
    }
};

template <typename T>
struct ArrayTraits<T,1,1> {
    typedef T Element;
    typedef boost::mpl::int_<1> ND;
    typedef boost::mpl::int_<1> RMC;
    typedef Element * Iterator;
    typedef Element & Reference;
    typedef Element Value;
    typedef detail::Core<1> Core;
    typedef typename Core::ConstPtr CorePtr;

    static Reference makeReference(Element * data, CorePtr const & core) {
        return *data;
    }
    static Iterator makeIterator(Element * data, CorePtr const & core, Offset stride) {
        return data;
    }
    static void fill(Iterator iter, Iterator const & end, Element value) {
        std::fill(iter, end, value);
    }
};

template <typename T>
struct ArrayTraits<T,1,-1> {
    typedef T Element;
    typedef boost::mpl::int_<1> ND;
    typedef boost::mpl::int_<-1> RMC;
    typedef Element * Iterator;
    typedef Element & Reference;
    typedef Element Value;
    typedef detail::Core<1> Core;
    typedef typename Core::ConstPtr CorePtr;

    static Reference makeReference(Element * data, CorePtr const & core) {
        return *data;
    }
    static Iterator makeIterator(Element * data, CorePtr const & core, Offset stride) {
        return data;
    }
    static void fill(Iterator iter, Iterator const & end, Element value) {
        std::fill(iter, end, value);
    }
};

template <typename T, int N, int C>
struct ExpressionTraits< Array<T,N,C> > : public ArrayTraits<T,N,C> {
    typedef Array<T,N,C> Self;
    typedef boost::mpl::false_ IsScalar;
};

template <typename T, int N, int C>
struct ExpressionTraits< ArrayRef<T,N,C> > : public ArrayTraits<T,N,C> {
    typedef ArrayRef<T,N,C> Self;
    typedef boost::mpl::false_ IsScalar;
};

} // namespace ndarray

#endif // !NDARRAY_ArrayTraits_h_INCLUDED
