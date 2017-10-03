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
#ifndef NDARRAY_BP_ArrayRef_h_INCLUDED
#define NDARRAY_BP_ArrayRef_h_INCLUDED

#include "ndarray/bp/Array.h"

namespace ndarray {

template <typename T, int N, int C>
class ToBoostPython< ArrayRef<T,N,C> > {
public:

    typedef boost::numpy::ndarray result_type;

    static boost::numpy::ndarray apply(ArrayRef<T,N,C> const & array) {
        return ToBoostPython< Array<T,N,C> >::apply(array);
    }

};

template <typename T, int N, int C>
class FromBoostPython< ArrayRef<T,N,C> > {
public:

    explicit FromBoostPython(boost::python::object const & input) : _impl(input) {}

    bool convertible() { return _impl.convertible(); }

    ArrayRef<T,N,C> operator()() { return ArrayRef<T,N,C>(_impl()); }

private:
    FromBoostPython< Array<T,N,C> > _impl;
};

} // namespace ndarray

namespace boost { namespace numpy {

template <typename T, int N, int C>
numpy::ndarray array(::ndarray::ArrayRef<T,N,C> const & arg) {
    return ::ndarray::ToBoostPython< ::ndarray::ArrayRef<T,N,C> >::apply(arg);
}

}} // namespace boost::numpy

#endif // !NDARRAY_BP_ArrayRef_h_INCLUDED
