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
#ifndef NDARRAY_BP_Vector_h_INCLUDED
#define NDARRAY_BP_Vector_h_INCLUDED

#include "boost/python.hpp"
#include "ndarray/Vector.h"
#include "ndarray/bp_fwd.h"

namespace ndarray {

template <typename T, int N>
class ToBoostPython< Vector<T,N> > {
public:

    typedef boost::python::tuple result_type;

    static boost::python::tuple apply(Vector<T,N> const & x) {
        boost::python::handle<> t(PyTuple_New(N));
        for (int n=0; n<N; ++n) {
            boost::python::object item(x[n]);
            Py_INCREF(item.ptr());
            PyTuple_SET_ITEM(t.get(), n, item.ptr());
        }
        return boost::python::tuple(t);
    }

};

template <typename T, int N>
class FromBoostPython< Vector<T,N> > {
public:

    explicit FromBoostPython(boost::python::object const & input_) : input(input_) {}

    bool convertible() {
        try {
            boost::python::tuple t(input);
            if (len(t) != N) return false;
            input = t;
        } catch (boost::python::error_already_set) {
            boost::python::handle_exception();
            PyErr_Clear();
            return false;
        }
        return true;
    }

    Vector<T,N> operator()() {
        boost::python::tuple t = boost::python::extract<boost::python::tuple>(input);
        if (len(t) != N) {
            PyErr_SetString(PyExc_ValueError, "Incorrect size for ndarray::Vector.");
            boost::python::throw_error_already_set();
        }
        Vector<T,N> r;
        for (int n=0; n<N; ++n) {
            r[n] = boost::python::extract<T>(t[n]);
        }
        return r;
    }

    boost::python::object input;
};

} // namespace ndarray

#endif // !NDARRAY_BP_Vector_h_INCLUDED
