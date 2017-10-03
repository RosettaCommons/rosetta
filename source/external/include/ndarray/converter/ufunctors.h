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
#ifndef NDARRAY_CONVERTER_ufunctors_h_INCLUDED
#define NDARRAY_CONVERTER_ufunctors_h_INCLUDED

/**
 *  \file ndarray/converter/ufunctors.h
 *  @brief Python wrappers to create numpy ufunc objects from C++ function objects.
 */

#include "ndarray/converter/numpy.h"

namespace ndarray {

template <typename TUnaryFunctor,
          typename TArgument=typename TUnaryFunctor::argument_type,
          typename TResult=typename TUnaryFunctor::result_type>
struct PyUnaryUFunctor {

    static PyObject* _call_(TUnaryFunctor const& self, PyObject* input, PyObject* output) {
        PyObject* input_array = PyArray_FROM_OTF(input,detail::NumpyTraits<TArgument>::getCode(),
                                                 NPY_ALIGNED);
        if (input_array == NULL) return NULL;
        PyObject* output_array;
        if (output == NULL || output == Py_None) {
            output_array = PyArray_SimpleNew(PyArray_NDIM(input_array),
                                             PyArray_DIMS(input_array),
                                             detail::NumpyTraits<TResult>::getCode());
        } else {
            output_array = PyArray_FROM_OTF(output,detail::NumpyTraits<TResult>::getCode(),
                                            NPY_ALIGNED | NPY_WRITEABLE | NPY_UPDATEIFCOPY);
        }
        if (output_array == NULL) {
            Py_DECREF(input_array);
            return NULL;
        }
        PyObject* iter = PyArray_MultiIterNew(2,input_array,output_array);
        if (iter == NULL) {
            Py_DECREF(input_array);
            Py_DECREF(output_array);
            return NULL;
        }
        int size = ((PyArrayMultiIterObject*)(iter))->size;
        while (size--) {
            TArgument* arg = (TArgument*)PyArray_MultiIter_DATA(iter,0);
            TResult* res = (TResult*)PyArray_MultiIter_DATA(iter,1);
            *res = self(*arg);
            PyArray_MultiIter_NEXT(iter);
        }
        Py_DECREF(input_array);
        Py_DECREF(iter);
        return PyArray_Return(reinterpret_cast<PyArrayObject*>(output_array));
    }

};


template <typename TBinaryFunctor,
          typename TArgument1=typename TBinaryFunctor::first_argument_type,
          typename TArgument2=typename TBinaryFunctor::second_argument_type,
          typename TResult=typename TBinaryFunctor::result_type>
struct PyBinaryUFunctor {

    static PyObject* _call_(TBinaryFunctor const& self, PyObject* input1, PyObject* input2,
                              PyObject* output) {
        PyObject* input1_array = PyArray_FROM_OTF(input1,detail::NumpyTraits<TArgument1>::getCode(),
                                                  NPY_ALIGNED);
        PyObject* input2_array = PyArray_FROM_OTF(input2,detail::NumpyTraits<TArgument1>::getCode(),
                                                  NPY_ALIGNED);
        if (input1_array == NULL || input2_array == NULL) {
            Py_XDECREF(input1_array);
            Py_XDECREF(input2_array);
            return NULL;
        }
        PyObject* output_array;
        if (output == NULL || output == Py_None) {
            PyObject* tmp = PyArray_MultiIterNew(2,input1_array,input2_array);
            if (tmp == NULL) {
                Py_XDECREF(input1_array);
                Py_XDECREF(input2_array);
                return NULL;
            }
            PyArrayMultiIterObject* tmp_iter = (PyArrayMultiIterObject*)tmp;
            output_array = PyArray_SimpleNew(tmp_iter->nd,tmp_iter->dimensions,
                                             detail::NumpyTraits<TResult>::getCode());
            Py_DECREF(tmp);
        } else {
            output_array = PyArray_FROM_OTF(output,detail::NumpyTraits<TResult>::getCode(),
                                            NPY_ALIGNED | NPY_WRITEABLE | NPY_UPDATEIFCOPY);
        }
        if (output_array == NULL) {
            Py_DECREF(input1_array);
            Py_DECREF(input2_array);
            return NULL;
        }
        PyObject* iter = PyArray_MultiIterNew(3,input1_array,input2_array,output_array);
        if (iter == NULL) {
            Py_DECREF(input1_array);
            Py_DECREF(input2_array);
            Py_DECREF(output_array);
            return NULL;
        }
        int size = ((PyArrayMultiIterObject*)(iter))->size;
        while (size--) {
            TArgument1* arg1 = (TArgument1*)PyArray_MultiIter_DATA(iter,0);
            TArgument2* arg2 = (TArgument2*)PyArray_MultiIter_DATA(iter,1);
            TResult* res = (TResult*)PyArray_MultiIter_DATA(iter,2);
            *res = self(*arg1,*arg2);
            PyArray_MultiIter_NEXT(iter);
        }
        Py_DECREF(input1_array);
        Py_DECREF(input2_array);
        Py_DECREF(iter);
        return PyArray_Return(reinterpret_cast<PyArrayObject*>(output_array));
    }

};

} // namespace ndarray

#endif // !NDARRAY_CONVERTER_ufunctors_h_INCLUDED
