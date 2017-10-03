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
#ifndef NDARRAY_CONVERTER_eigen_h_INCLUDED
#define NDARRAY_CONVERTER_eigen_h_INCLUDED

/**
 *  @file ndarray/converter/eigen.h
 *  @brief Python C-API conversions for Eigen matrices.
 */

#include "ndarray/eigen.h"
#include "ndarray/converter.h"
#include "ndarray/converter/eigen.h"

namespace ndarray {

/**
 *  @ingroup ndarrayPythonGroup
 *  @brief Specialization of PyConverter for EigenView.
 */
template <typename T, int N, int C, typename XprKind_, int Rows_, int Cols_>
struct PyConverter< EigenView<T,N,C,XprKind_,Rows_,Cols_> > {

    static bool fromPythonStage1(PyPtr & p) {
        // add or remove dimensions with size one so we have the right number of dimensions
        if (PyArray_Check(p.get())) {
            if ((Rows_ == 1 || Cols_ == 1) && N == 2) {
                npy_intp shape[2] = { -1, -1 };
                if (Rows_ == 1) {
                    shape[0] = 1;
                } else {
                    shape[1] = 1;
                }
                PyArray_Dims dims = { shape, 2 };
                PyPtr r(PyArray_Newshape(reinterpret_cast<PyArrayObject*>(p.get()), &dims, NPY_ANYORDER));
                if (!r) return false;
                p.swap(r);
            } else if (N == 1) {
                PyPtr r(PyArray_Squeeze(reinterpret_cast<PyArrayObject*>(p.get())));
                if (!r) return false;
                p.swap(r);
            }
        } // else let the Array converter raise the exception
        if (!PyConverter< Array<T,N,C> >::fromPythonStage1(p)) return false;
        // check whether the size is correct if it's static
        if (N == 2) {
            if (Rows_ != Eigen::Dynamic && PyArray_DIM(p.get(), 0) != Rows_) {
                PyErr_SetString(PyExc_ValueError, "incorrect number of rows for matrix");
                return false;
            }
            if (Cols_ != Eigen::Dynamic && PyArray_DIM(p.get(), 1) != Cols_) {
                PyErr_SetString(PyExc_ValueError, "incorrect number of columns for matrix");
                return false;
            }
        } else {
            int requiredSize = Rows_ * Cols_;
            if (requiredSize != Eigen::Dynamic && PyArray_SIZE(p.get()) != requiredSize) {
                PyErr_SetString(PyExc_ValueError, "incorrect number of elements for vector");
                return false;
            }
        }
        return true;
    }

    static bool fromPythonStage2(
        PyPtr const & input,
        EigenView<T,N,C,XprKind_,Rows_,Cols_> & output
    ) {
        Array<T,N,C> array;
        if (!PyConverter< Array<T,N,C> >::fromPythonStage2(input, array)) return false;
        output.reset(array);
        return true;
    }

    static PyObject * toPython(EigenView<T,N,C,XprKind_,Rows_,Cols_> const & m, PyObject * owner=NULL) {
        PyPtr r(PyConverter< Array<T,N,C> >::toPython(m.shallow(), owner));
        if (!r) return NULL;
        PyPtr p(PyArray_Squeeze(reinterpret_cast<PyArrayObject*>(r.get())));
        Py_XINCREF(p.get());
        return p.get();
    }

};

namespace detail {

/**
 *  @internal @ingroup ndarrayPythonInternalGroup
 *  @brief Implementations for PyConverter for Eigen objects.
 */
template <typename Matrix>
class EigenPyConverter : public detail::PyConverterBase<Matrix> {
    typedef typename SelectEigenView<Matrix>::Type OutputView;
    typedef typename SelectEigenView<Matrix const, false>::Type InputView;
public:

    static PyObject * toPython(Matrix const & input, PyObject * owner = NULL) {
        return PyConverter<OutputView>::toPython(ndarray::copy(input), owner);
    }

    static PyTypeObject const * getPyType() {
        return &PyArray_Type;
    }

    static bool fromPythonStage1(PyPtr & p) {
        return PyConverter<InputView>::fromPythonStage1(p);
    }

    static bool fromPythonStage2(PyPtr const & p, Matrix & output) {
        InputView v;
        if (!PyConverter<InputView>::fromPythonStage2(p, v)) return false;
        output = v;
        return true;
    }

};

} // namespace detail

/**
 *  @ingroup ndarrayPythonGroup
 *  @brief Specialization of PyConverter for Eigen::Matrix.
 */
template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct PyConverter< Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols> >
    : public detail::EigenPyConverter< Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols> >
{};

/**
 *  @ingroup ndarrayPythonGroup
 *  @brief Specialization of PyConverter for Eigen::Array.
 */
template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct PyConverter< Eigen::Array<Scalar,Rows,Cols,Options,MaxRows,MaxCols> >
    : public detail::EigenPyConverter< Eigen::Array<Scalar,Rows,Cols,Options,MaxRows,MaxCols> >
{};

} // namespace ndarray

#endif // !NDARRAY_CONVERTER_eigen_h_INCLUDED
