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
#ifndef NDARRAY_CONVERTER_PyConverter_h_INCLUDED
#define NDARRAY_CONVERTER_PyConverter_h_INCLUDED

/**
 *  @file ndarray/converter/PyConverter.h
 *  @brief Python C-API conversions for standard numeric types.
 */
#include <Python.h>

#include <boost/intrusive_ptr.hpp>
#include <complex>

#include<ndarray/converter/pointer.h>

namespace ndarray {

template <typename T> struct PyConverter;

namespace detail {

/**
 *  @internal @ingroup ndarrayPythonInternalGroup
 *  @brief Non-specialized base class for PyConverter.
 */
template <typename T>
struct PyConverterBase {

    /**
     *  @brief Check if a Python object might be convertible to T.
     *
     *  \return true if the conversion may be successful, and
     *  false if it definitely is not.  Will not raise a Python
     *  exception.
     *
     *  This is mostly useful for wrapper generators like CONVERTER or
     *  Boost.Python, which could use matches() to check if an
     *  Python arguments match a particular signature for an
     *  overloaded C++ function.
     *
     *  \sa PyConverter<T>::fromPythonStage1()
     */
    static bool matches(
        PyObject * arg // input Python object (borrowed)
    ) {
        PyPtr p(arg,true);
        if (!PyConverter<T>::fromPythonStage1(p)) {
            PyErr_Clear();
            return false;
        }
        return true;
    }

    /**
     *  @brief Convert a Python object to a C++ object.
     *
     *  fromPython() is appropriate for use with the PyArg_ParseX
     *  Python C-API functions as a "O&" converter function.
     *
     *  \return true if the conversion was successful, and false
     *  otherwise (with a Python exception set).
     *
     *  \sa PyConverter<T>::fromPythonStage1()
     *  \sa PyConverter<T>::fromPythonStage2()
     **/
    static int fromPython(
        PyObject * arg,  // input Python object (borrowed)
        T * output       // pointer to an existing C++
    ) {
        PyPtr p(arg,true);
        if (!PyConverter<T>::fromPythonStage1(p)) return false;
        return PyConverter<T>::fromPythonStage2(arg,*output);
    }

};

/**
 *  @internal @ingroup ndarrayPythonInternalGroup
 *
 *  Intermediate base class for PyConverter implementations for non-bool integer types,
 *  providing a common stage 1 test.
 */
template <typename T>
struct PyIntConverterBase : public PyConverterBase<T> {

    static bool fromPythonStage1(PyPtr & input) {
        if (
#if PY_MAJOR_VERSION <= 2
            !PyInt_Check(input.get()) &&
#endif
            !PyLong_Check(input.get())
        ) {
            PyPtr s(PyObject_Repr(input.get()));
            if (!s) return false;
            char * cs = PyBytes_AsString(s.get());
            if (!cs) return false;
            PyErr_Format(PyExc_TypeError, "'%s' is not a valid C++ integer value.", cs);
        }
        return true;
    }

};

/**
 *  @internal @ingroup ndarrayPythonInternalGroup
 *
 *  Implementation of PyConverter for non-bool integer types, specialized on
 *  whether the type is unsigned and whether it is smaller than, equal to, or
 *  greater than long (which is the maximum for Python's 'int') in size.
 */
template <
    typename T,
    bool isUnsigned=boost::is_unsigned<T>::value,
    int cmpToLong=(sizeof(T) < sizeof(long)) ? -1 : ((sizeof(T)==sizeof(long)) ? 0 : 1)
    >
struct PyIntConverter;

/**
 *  @internal @ingroup ndarrayPythonInternalGroup
 *
 *  Specialized implementation of PyConverter for signed integers smaller than long in size.
 */
template <typename T>
struct PyIntConverter<T,false,-1> : public PyIntConverterBase<T> {

    static bool fromPythonStage2(PyPtr const & input, T & output) {
        NDARRAY_ASSERT(input);
        output = PyLong_AsLong(input.get());
        if (PyErr_Occurred()) return false; // could get OverflowError here.
        return true;
    }

    static PyObject * toPython(T input) {
        return PyLong_FromLong(input);
    }

    static PyTypeObject const * getPyType() { return &PyLong_Type; }

};

/**
 *  @internal @ingroup ndarrayPythonInternalGroup
 *
 *  Specialized implementation of PyConverter for unsigned integers smaller than long in size.
 */
template <typename T>
struct PyIntConverter<T,true,-1> : public PyIntConverterBase<T> {

    static bool fromPythonStage2(PyPtr const & input, T & output) {
        NDARRAY_ASSERT(input);
        output = PyLong_AsLong(input.get());
        if (PyErr_Occurred()) return false; // could get OverflowError here.
        return true;
    }

    static PyObject * toPython(T input) {
        return PyLong_FromLong(input);
    }

    static PyTypeObject const * getPyType() { return &PyLong_Type; }

};

/**
 *  @internal @ingroup ndarrayPythonInternalGroup
 *
 *  Specialized implementation of PyConverter for signed integers equal to long in size.
 */
template <typename T>
struct PyIntConverter<T,false,0> : public PyIntConverterBase<T> {

    static bool fromPythonStage2(PyPtr const & input, T & output) {
        NDARRAY_ASSERT(input);
        output = PyLong_AsLong(input.get());
        if (PyErr_Occurred()) return false; // could get OverflowError here.
        return true;
    }

    static PyObject * toPython(T input) {
        return PyLong_FromLong(input);
    }

    static PyTypeObject const * getPyType() { return &PyLong_Type; }

};

/**
 *  @internal @ingroup ndarrayPythonInternalGroup
 *
 *  Specialized implementation of PyConverter for unsigned integers equal to long in size.
 */
template <typename T>
struct PyIntConverter<T,true,0> : public PyIntConverterBase<T> {

    static bool fromPythonStage2(PyPtr const & input, T & output) {
        NDARRAY_ASSERT(input);
        output = PyLong_AsUnsignedLong(input.get());
        if (PyErr_Occurred()) return false; // could get OverflowError here.
        return true;
    }

    static PyObject * toPython(T input) {
        return PyLong_FromUnsignedLong(input);
    }

    static PyTypeObject const * getPyType() { return &PyLong_Type; }

};

/**
 *  @internal @ingroup ndarrayPythonInternalGroup
 *
 *  Specialized implementation of PyConverter for signed integers greater than long in size.
 */
template <typename T>
struct PyIntConverter<T,false,1> : public PyIntConverterBase<T> {

    static bool fromPythonStage2(PyPtr const & input, T & output) {
        NDARRAY_ASSERT(input);
        output = PyLong_AsLongLong(input.get());
        if (PyErr_Occurred()) return false; // could get OverflowError here.
        return true;
    }

    static PyObject * toPython(T input) {
        return PyLong_FromLongLong(input);
    }

    static PyTypeObject const * getPyType() { return &PyLong_Type; }

};

/**
 *  @internal @ingroup ndarrayPythonInternalGroup
 *
 *  Specialized implementation of PyConverter for unsigned integers greater than long in size.
 */
template <typename T>
struct PyIntConverter<T,true,1> : public PyIntConverterBase<T> {

    static bool fromPythonStage2(PyPtr const & input, T & output) {
        NDARRAY_ASSERT(input);
#ifndef _MSC_VER
        output = static_cast<T>(PyLong_AsUnsignedLongLong(input.get()));
#else
        // Christoph Lassner reports that the above segfaults on Windows; until we have a better idea
        // what triggers that, it's safer just to use PyLong_AsLongLong and just not support the full
        // range of unsigned long long values.
        PY_LONG_LONG tmp = PyLong_AsLongLong(input.get());
        if (tmp < 0) {
            PyErr_SetString(PyExc_OverflowError, "Negative integer in conversion to unsigned");
        } else {
            output = static_cast<T>(tmp);
        }
#endif
        if (PyErr_Occurred()) return false; // could get OverflowError here.
        return true;
    }

    static PyObject * toPython(T input) {
        return PyLong_FromUnsignedLongLong(input);
    }

    static PyTypeObject const * getPyType() { return &PyLong_Type; }

};

} // namespace ndarray::detail

/**
 *  @ingroup ndarrayPythonGroup
 *  @brief A class providing Python conversion functions for T.
 *
 *  Undocumented specializations exist for bool, (unsigned) short, (unsigned) int, (unsigned) long,
 *  (unsigned long long), float, double, std::complex<float>, std::complex<double>, and std::string.
 */
template <typename T>
struct PyConverter : public detail::PyConverterBase<T> {

    /**
     *  @brief Convert a C++ object to a new Python object.
     *
     *  \return A new Python object, or NULL on failure (with
     *  a Python exception set).
     */
    static PyObject * toPython(
        T const & input ///< Input C++ object.
    );

    /**
     *  @brief Return the Python TypeObject that corresponds to
     *  the object the toPython() function returns.
     */
    static PyTypeObject const * getPyType();

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
    );

    /**
     *  @brief Complete a Python to C++ conversion begun
     *  with fromPythonStage1().
     *
     *  \return true if the conversion was successful,
     *  and false otherwise (with a Python exception set).
     */
    static bool fromPythonStage2(
        PyPtr const & p, ///< A Python object processed by fromPythonStage1().
        T & output       ///< The output C++ object.
    );

};

/// \cond SPECIALIZATIONS

template <>
struct PyConverter<bool> : public detail::PyConverterBase<bool> {

    static bool fromPythonStage1(PyPtr & input) {
        if (!PyBool_Check(input.get())) {
            PyPtr s(PyObject_Repr(input.get()));
            if (!s) return false;
            char * cs = PyBytes_AsString(s.get());
            if (!cs) return false;
            PyErr_Format(PyExc_TypeError,"'%s' is not a valid C++ bool value.",cs);
        }
        return true;
    }

    static bool fromPythonStage2(PyPtr const & input, bool & output) {
        NDARRAY_ASSERT(input);
        output = (input.get() == Py_True);
        return true;
    }

    static PyObject * toPython(bool input) {
        if (input) Py_RETURN_TRUE;
        Py_RETURN_FALSE;
    }

    static PyTypeObject const * getPyType() { return &PyBool_Type; }
};

template <> struct PyConverter<short> : public detail::PyIntConverter<short> {};
template <> struct PyConverter<unsigned short> : public detail::PyIntConverter<unsigned short> {};
template <> struct PyConverter<int> : public detail::PyIntConverter<int> {};
template <> struct PyConverter<unsigned int> : public detail::PyIntConverter<unsigned int> {};
template <> struct PyConverter<long> : public detail::PyIntConverter<long> {};
template <> struct PyConverter<unsigned long> : public detail::PyIntConverter<unsigned long> {};
template <> struct PyConverter<long long> : public detail::PyIntConverter<long long> {};
template <> struct PyConverter<unsigned long long> : public detail::PyIntConverter<unsigned long long> {};

template <>
struct PyConverter<float> : public detail::PyConverterBase<float> {

    static bool fromPythonStage1(PyPtr & input) {
        if (!PyFloat_Check(input.get())) {
            PyPtr s(PyObject_Repr(input.get()));
            if (!s) return false;
            char * cs = PyBytes_AsString(s.get());
            if (!cs) return false;
            PyErr_Format(PyExc_TypeError,"'%s' is not a valid C++ float value.",cs);
        }
        return true;
    }

    static bool fromPythonStage2(PyPtr const & input, float & output) {
        NDARRAY_ASSERT(input);
        output = PyFloat_AsDouble(input.get());
        return true;
    }

    static PyObject * toPython(float input) {
        return PyFloat_FromDouble(input);
    }

    static PyTypeObject const * getPyType() { return &PyFloat_Type; }
};

template <>
struct PyConverter<double> : public detail::PyConverterBase<double> {

    static bool fromPythonStage1(PyPtr & input) {
        if (!PyFloat_Check(input.get())) {
            PyPtr s(PyObject_Repr(input.get()));
            if (!s) return false;
            char * cs = PyBytes_AsString(s.get());
            if (!cs) return false;
            PyErr_Format(PyExc_TypeError,"'%s' is not a valid C++ double value.",cs);
        }
        return true;
    }

    static bool fromPythonStage2(PyPtr const & input, double & output) {
        NDARRAY_ASSERT(input);
        output = PyFloat_AsDouble(input.get());
        return true;
    }

    static PyObject * toPython(double input) {
        return PyFloat_FromDouble(input);
    }

    static PyTypeObject const * getPyType() { return &PyFloat_Type; }
};

template <typename U>
struct PyConverter< std::complex<U> > : public detail::PyConverterBase< std::complex<U> > {

    static bool fromPythonStage1(PyPtr & input) {
        if (!PyComplex_Check(input.get())) {
            PyPtr s(PyObject_Repr(input.get()));
            if (!s) return false;
            char * cs = PyBytes_AsString(s.get());
            if (!cs) return false;
            PyErr_Format(PyExc_TypeError,"'%s' is not a valid C++ complex value.",cs);
        }
        return true;
    }

    static bool fromPythonStage2(PyPtr const & input, std::complex<U> & output) {
        NDARRAY_ASSERT(input);
        output.real() = PyComplex_RealAsDouble(input.get());
        output.imag() = PyComplex_ImagAsDouble(input.get());
        return true;
    }

    static PyObject * toPython(std::complex<U> const & input) {
        return PyComplex_FromDoubles(input.real(),input.imag());
    }

    static PyTypeObject const * getPyType() { return &PyComplex_Type; }
};

template <>
struct PyConverter< std::string > : public detail::PyConverterBase<std::string> {

    static bool fromPythonStage1(PyPtr & input) {
        if (!PyBytes_Check(input.get())) {
            PyPtr s(PyObject_Repr(input.get()));
            if (!s) return false;
            char * cs = PyBytes_AsString(s.get());
            if (!cs) return false;
            PyErr_Format(PyExc_TypeError,"'%s' is not a valid C++ string value.",cs);
        }
        return true;
    }

    static bool fromPythonStage2(PyPtr const & input, std::string & output) {
        NDARRAY_ASSERT(input);
        char * buf = 0;
        Py_ssize_t size = 0;
        if (PyBytes_AsStringAndSize(input.get(),&buf,&size) == -1) return false;
        output = std::string(buf,size);
        return true;
    }

    static PyObject * toPython(std::string const & input) {
        return PyBytes_FromStringAndSize(input.data(),input.size());
    }

    static PyTypeObject const * getPyType() { return &PyBytes_Type; }
};

/// \endcond

} // namespace ndarray

#endif // !NDARRAY_CONVERTER_PyConverter_h_INCLUDED
