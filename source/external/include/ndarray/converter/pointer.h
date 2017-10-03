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
#pragma once

/**
 *  @file ndarray/converter/pointer.h
 *  @brief PyPtr and PyCapsule support functions.
 */

#include "Python.h"
#include "ndarray.h"

#ifndef DOXYGEN
//namespace boost {
inline void intrusive_ptr_add_ref(PyObject * obj) { Py_INCREF(obj); }
inline void intrusive_ptr_release(PyObject * obj) { Py_DECREF(obj); }
//}
#endif

namespace ndarray {

/**
 *  @brief A reference-counting smart pointer for PyObject.
 */
typedef boost::intrusive_ptr<PyObject> PyPtr;

}


namespace ndarray {
namespace detail {

/**
 *  @internal @ingroup ndarrayPythonInternalGroup
 *  @brief A destructor for a Python CObject that owns a shared_ptr.
 */
inline void destroyCapsule(PyObject * p) {
    void * m = PyCapsule_GetPointer(p, "ndarray.Manager");
    ndarray::Manager::Ptr * b = reinterpret_cast<ndarray::Manager::Ptr*>(m);
    delete b;
}

}}
