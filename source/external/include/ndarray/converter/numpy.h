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
#ifndef NDARRAY_CONVERTER_numpy_h_INCLUDED
#define NDARRAY_CONVERTER_numpy_h_INCLUDED

/**
 *  @file ndarray/converter/numpy.h
 *  @brief Python C-API conversions between ndarray and numpy.
 */

#include "Python.h"
#include "ndarray.h"
#include "ndarray/converter/pointer.h"
#include "ndarray/converter/PyConverter.h"

namespace ndarray {
namespace detail {

/**
 *  @internal @ingroup ndarrayPythonInternalGroup
 *  @brief Traits class that specifies Numpy typecodes for numeric types.
 */
template <typename T> struct NumpyTraits {
    static int getCode();
};

/// \cond SPECIALIZATIONS

template <> struct NumpyTraits<bool> {
    static int getCode() {
	if (sizeof(bool)==sizeof(npy_bool)) return NPY_BOOL;
	if (sizeof(bool)==1) return NPY_UBYTE;
	if (sizeof(bool)==2 && sizeof(short)==2) return NPY_USHORT;
	if (sizeof(bool)==4 && sizeof(int)==4) return NPY_UINT;
	assert(false);
	return 0;
    }
};

template <> struct NumpyTraits<npy_ubyte> { static int getCode() { return NPY_UBYTE; } };
template <> struct NumpyTraits<npy_byte> { static int getCode() { return NPY_BYTE; } };
template <> struct NumpyTraits<npy_ushort> { static int getCode() { return NPY_USHORT; } };
template <> struct NumpyTraits<npy_short> { static int getCode() { return NPY_SHORT; } };
// NPY_INT is on Windows a virtual flag that is never actually used. It must be
// checked against the platform dtype.
// see http://mail.scipy.org/pipermail/numpy-discussion/2010-June/051057.html.
template <> struct NumpyTraits<npy_uint> {
    static int getCode() {
#ifdef _MSC_VER
        switch(sizeof(int)) {
        case 1: return NPY_UBYTE;
        case 2: return NPY_USHORT;
        case 4: return NPY_ULONG;
        case 8: return NPY_ULONGLONG;
        // no datatype here...
        default: throw std::exception();
        }
#else
        return NPY_UINT;
#endif
    }
};
template <> struct NumpyTraits<npy_int> {
    static int getCode() {
#ifdef _MSC_VER
        switch(sizeof(int)) {
        case 1: return NPY_BYTE;
        case 2: return NPY_SHORT;
        case 4: return NPY_LONG;
        case 8: return NPY_LONGLONG;
        // no datatype here...
        default: throw std::exception();
        }
#else
    return NPY_INT;
#endif
    }
};
template <> struct NumpyTraits<npy_ulong> { static int getCode() { return NPY_ULONG; } };
template <> struct NumpyTraits<npy_long> { static int getCode() { return NPY_LONG; } };
template <> struct NumpyTraits<npy_ulonglong> { static int getCode() { return NPY_ULONGLONG; } };
template <> struct NumpyTraits<npy_longlong> { static int getCode() { return NPY_LONGLONG; } };
template <> struct NumpyTraits<npy_float> { static int getCode() { return NPY_FLOAT; } };
template <> struct NumpyTraits<npy_double> { static int getCode() { return NPY_DOUBLE; } };
#if (npy_double != npy_longdouble)
  template <> struct NumpyTraits<npy_longdouble> { static int getCode() { return NPY_LONGDOUBLE; } };
#endif
template <> struct NumpyTraits<npy_cfloat> { static int getCode() { return NPY_CFLOAT; } };
template <> struct NumpyTraits<npy_cdouble> { static int getCode() { return NPY_CDOUBLE; } };
template <> struct NumpyTraits<npy_clongdouble> { static int getCode() { return NPY_CLONGDOUBLE; } };

template <> struct NumpyTraits<std::complex<float> > {
    static int getCode() { assert(sizeof(std::complex<float>)==sizeof(npy_cfloat)); return NPY_CFLOAT; }
};

template <> struct NumpyTraits<std::complex<double> > {
    static int getCode() { assert(sizeof(std::complex<double>)==sizeof(npy_cdouble)); return NPY_CDOUBLE; }
};

template <> struct NumpyTraits<std::complex<long double> > {
    static int getCode() {
	assert(sizeof(std::complex<long double>)==sizeof(npy_clongdouble));
	return NPY_CLONGDOUBLE;
    }
};

/// \endcond

} // namespace ndarray::detail

/**
 *  @ingroup ndarrayPythonGroup
 *  @brief A traits class providing Python conversion functions for Array.
 *
 *  This specialization, for Array, adds addititional optional arguments
 *  to the toPython() conversion member function.
 */
template <typename T, int N, int C>
struct PyConverter< Array<T,N,C> > : public detail::PyConverterBase< Array<T,N,C> > {
    typedef typename Array<T,N,C>::Element Element;
    typedef typename boost::remove_const<Element>::type NonConst;

    /**
     *  @brief Check if a Python object is convertible to T
     *  and optionally begin the conversion by replacing the
     *  input with an intermediate.
     *
     *  \return true if a conversion may be possible, and
     *  false if it is not (with a Python exception set).
     */
    static bool fromPythonStage1(
        PyPtr & p /**< On input, a Python object to be converted.
                   *   On output, a Python object to be passed to
                   *   fromPythonStage2().
                   */
    ) {
        if (!PyArray_Check(p.get())) {
            PyErr_SetString(PyExc_TypeError, "numpy.ndarray argument required");
            return false;
        }
        int actualType = PyArray_TYPE(p.get());
        int requiredType = detail::NumpyTraits<NonConst>::getCode();
        if (actualType != requiredType) {
            PyErr_SetString(PyExc_ValueError, ("numpy.ndarray argument has incorrect data type"));
            return false;
        }
        if (PyArray_NDIM(p.get()) != N) {
            PyErr_SetString(PyExc_ValueError, "numpy.ndarray argument has incorrect number of dimensions");
            return false;
        }
        bool writeable = !boost::is_const<Element>::value;
        if (writeable && !(PyArray_FLAGS(p.get()) & NPY_WRITEABLE)) {
            PyErr_SetString(PyExc_TypeError, "numpy.ndarray argument must be writeable");
            return false;
        }
        if (C > 0) {
            Offset requiredStride = sizeof(Element);
            for (int i = 0; i < C; ++i) {
                Offset actualStride = PyArray_STRIDE(p.get(), N-i-1);
                if (actualStride != requiredStride) {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "numpy.ndarray does not have enough row-major contiguous dimensions"
                    );
                    return false;
                }
                requiredStride *= PyArray_DIM(p.get(), N-i-1);
            }
        } else if (C < 0) {
            int requiredStride = sizeof(Element);
            for (int i = 0; i < -C; ++i) {
                Offset actualStride = PyArray_STRIDE(p.get(), i);
                if (actualStride != requiredStride) {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "numpy.ndarray does not have enough column-major contiguous dimensions"
                    );
                    return false;
                }
                requiredStride *= PyArray_DIM(p.get(), i);
            }
        }
        return true;
    }

    /**
     *  @brief Complete a Python to C++ conversion begun with fromPythonStage1().
     *
     *  The copy will be shallow if possible and deep if necessary to meet the data type
     *  and contiguousness requirements.  If a non-const array is required, the copy will
     *  always be shallow; if this is not possible, ValueError will be raised.
     *
     *  The output Array's shared_ptr owner attribute will own a reference to the numpy
     *  array that ultimately owns the data (either the original or the copy).
     *
     *  \return true on success, false on failure (with a Python exception set).
     */
    static bool fromPythonStage2(
        PyPtr const & input,  ///< Result of fromPythonStage1().
        Array<T,N,C> & output ///< Reference to existing output C++ object.
    ) {
        if (!(PyArray_FLAGS(input.get()) & NPY_ALIGNED)) {
            PyErr_SetString(PyExc_TypeError, "unaligned arrays cannot be converted to C++");
            return false;
        }
        Offset itemsize = sizeof(Element);
        for (int i = 0; i < N; ++i) {
            if ((PyArray_DIM(input.get(), i) > 1) && (PyArray_STRIDE(input.get(), i) % itemsize != 0)) {
                PyErr_SetString(
                    PyExc_TypeError,
                    "Cannot convert array to C++: strides must be an integer multiple of the element size"
                );
                return false;
            }
        }
        Vector<Size,N> shape;
        Vector<Offset,N> strides;
        std::copy(PyArray_DIMS(input.get()), PyArray_DIMS(input.get()) + N, shape.begin());
        std::copy(PyArray_STRIDES(input.get()), PyArray_STRIDES(input.get()) + N , strides.begin());
        for (int i = 0; i < N; ++i) strides[i] /= sizeof(Element);
        output = external(
            reinterpret_cast<Element*>(PyArray_DATA(input.get())),
            shape, strides, input
        );
        return true;
    }

    /**
     *  @brief Create a numpy.ndarray from an ndarray::Array.
     *
     *  The Array will be shallow-copied with reference counting if either
     *  m.getManager() is not empty or the optional owner argument is supplied;
     *  otherwise a deep copy will be made.
     *
     *  \return a new Python object, or NULL on failure (with
     *  a Python exception set).
     */
    static PyObject* toPython(
        Array<T,N,C> const & m, ///< The input Array object.
        PyObject* owner=NULL    /**< A Python object that owns the memory in the Array.
                                 *   If NULL, one will be constructed from m.getManager(). */
    ) {
        int flags = NPY_ALIGNED;
        if (C==N) flags |= NPY_C_CONTIGUOUS;
        bool writeable = !boost::is_const<Element>::value;
        if (writeable) flags |= NPY_WRITEABLE;
        npy_intp outShape[N];
        npy_intp outStrides[N];
        Vector<Size,N> inShape = m.getShape();
        Vector<Offset,N> inStrides = m.getStrides();
        std::copy(inShape.begin(), inShape.end(), outShape);
        for (int i = 0; i < N; ++i) outStrides[i] = inStrides[i] * sizeof(Element);
        PyPtr array(
            PyArray_New(
                &PyArray_Type, N, outShape, detail::NumpyTraits<NonConst>::getCode(),
                outStrides, const_cast<NonConst*>(m.getData()), sizeof(Element), flags, NULL
            ),
            false
        );
        if (!array) return NULL;
        if (!m.getManager() && owner == NULL) {
            flags = NPY_CARRAY_RO | NPY_ENSURECOPY | NPY_C_CONTIGUOUS;
            if (writeable) flags |= NPY_WRITEABLE;
            PyPtr r = PyArray_FROM_OF(array.get(),flags);
            if (!r) return NULL;
            array.swap(r);
        } else {
            if (owner != NULL) {
                Py_INCREF(owner);
            } else {
                owner = PyCapsule_New(
                    new Manager::Ptr(m.getManager()),
                    "ndarray.Manager",
                    detail::destroyCapsule
                );
            }
            reinterpret_cast<PyArrayObject*>(array.get())->base = owner;
        }
        Py_INCREF(array.get());
        return PyArray_Return(reinterpret_cast<PyArrayObject*>(array.get()));
    }

    static PyTypeObject const * getPyType() { return &PyArray_Type; }
};

} // namespace ndarray

#endif // !NDARRAY_CONVERTER_numpy_h_INCLUDED
