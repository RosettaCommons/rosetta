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
#ifndef NDARRAY_BP_Array_h_INCLUDED
#define NDARRAY_BP_Array_h_INCLUDED

#include "boost/numpy.hpp"
#include "ndarray.h"
#include "ndarray/bp_fwd.h"
#include <vector>

namespace ndarray {
namespace detail {

inline void destroyManagerCObject(void * p) {
#if PY_MAJOR_VERSION == 2
    Manager::Ptr * b = reinterpret_cast<Manager::Ptr*>(p);
#else
    Manager::Ptr * b = reinterpret_cast<Manager::Ptr*>(PyCapsule_GetPointer(reinterpret_cast<PyObject*>(p), 0));
#endif
    delete b;
}

inline boost::python::object makePyObject(Manager::Ptr const & x) {
    boost::intrusive_ptr< ExternalManager<boost::python::object> > y
        = boost::dynamic_pointer_cast< ExternalManager<boost::python::object> >(x);
    if (y) {
        return y->getOwner();
    }
#if PY_MAJOR_VERSION == 2
    boost::python::handle<> h(::PyCObject_FromVoidPtr(new Manager::Ptr(x), &destroyManagerCObject));
#else
    boost::python::handle<> h(::PyCapsule_New(new Manager::Ptr(x), 0, (PyCapsule_Destructor)&destroyManagerCObject));
#endif
    return boost::python::object(h);
}

} // namespace detail

template <typename T, int N, int C>
class ToBoostPython< Array<T,N,C> > {
public:

    typedef boost::numpy::ndarray result_type;

    static boost::numpy::ndarray apply(Array<T,N,C> const & array) {
        boost::numpy::dtype dtype
            = boost::numpy::dtype::get_builtin<typename boost::remove_const<T>::type>();
        boost::python::object owner = detail::makePyObject(array.getManager());
        Py_ssize_t itemsize = dtype.get_itemsize();
        ndarray::Vector<Size,N> shape_elements = array.getShape();
        ndarray::Vector<Offset,N> strides_elements = array.getStrides();
        std::vector<Py_intptr_t> shape_bytes(N);
        std::vector<Py_intptr_t> strides_bytes(N);
        for (int n=0; n<N; ++n) {
            shape_bytes[n] = shape_elements[n];
            strides_bytes[n] = strides_elements[n] * itemsize;
        }
        return boost::numpy::from_data(array.getData(), dtype, shape_bytes, strides_bytes, owner);
    }

};

template <typename T, int N, int C>
class FromBoostPython< Array<T,N,C> > {
public:

    explicit FromBoostPython(boost::python::object const & input_) : input(input_) {}

    bool convertible() {
        if (input.is_none()) return true;
        try {
            boost::numpy::ndarray array = boost::python::extract<boost::numpy::ndarray>(input);
            boost::numpy::dtype dtype
                = boost::numpy::dtype::get_builtin<typename boost::remove_const<T>::type>();
            boost::numpy::ndarray::bitflag flags = array.get_flags();
            if (dtype != array.get_dtype()) return false;
            if (N != array.get_nd()) return false;
            if (!boost::is_const<T>::value && !(flags & boost::numpy::ndarray::WRITEABLE)) return false;
            if (C > 0) {
                Offset requiredStride = sizeof(T);
                for (int i = 0; i < C; ++i) {
                    if ((array.shape(N-i-1) > 1) && (array.strides(N-i-1) != requiredStride)) {
                        return false;
                    }
                    requiredStride *= array.shape(N-i-1);
                }
            } else if (C < 0) {
                Offset requiredStride = sizeof(T);
                for (int i = 0; i < -C; ++i) {
                    if ((array.shape(i) > 1) && (array.strides(i) != requiredStride)) {
                        return false;
                    }
                    requiredStride *= array.shape(i);
                }
            }
        } catch (boost::python::error_already_set) {
            boost::python::handle_exception();
            PyErr_Clear();
            return false;
        }
        return true;
    }

    Array<T,N,C> operator()() {
        if (input.is_none()) return Array<T,N,C>();
        boost::numpy::ndarray array = boost::python::extract<boost::numpy::ndarray>(input);
        boost::numpy::dtype dtype = array.get_dtype();
        Py_ssize_t itemsize = dtype.get_itemsize();
        for (int i = 0; i < N; ++i) {
            if ((array.shape(i) > 1) && (array.strides(i) % itemsize != 0)) {
                PyErr_SetString(
                    PyExc_TypeError,
                    "Cannot convert array to C++: strides must be an integer multiple of the element size"
                );
                boost::python::throw_error_already_set();
            }
        }
        boost::python::object obj_owner = array.get_base();
        if (obj_owner.is_none()) {
            obj_owner = array;
        }
        Vector<Size,N> shape;
        Vector<Offset,N> strides;
        for (int i=0; i<N; ++i) {
            shape[i] = array.shape(i);
            if (shape[i] > 1) {
                strides[i] = array.strides(i) / itemsize;
            } else {
                strides[i] = 1;
            }
        }
        Array<T,N,C> r = ndarray::external(
            reinterpret_cast<T*>(array.get_data()), shape, strides, obj_owner
        );
        return r;
    }

    boost::python::object input;
};

} // namespace ndarray

namespace boost { namespace numpy {

template <typename T, int N, int C>
numpy::ndarray array(::ndarray::Array<T,N,C> const & arg) {
    return ::ndarray::ToBoostPython< ::ndarray::Array<T,N,C> >::apply(arg);
}

}} // namespace boost::numpy

#endif // !NDARRAY_BP_Array_h_INCLUDED
