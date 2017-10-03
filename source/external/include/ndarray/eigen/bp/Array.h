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
#ifndef NDARRAY_EIGEN_BP_Array_h_INCLUDED
#define NDARRAY_EIGEN_BP_Array_h_INCLUDED

#include "ndarray/eigen/bp/EigenView.h"

namespace ndarray {

template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int  MaxCols>
class ToBoostPython< Eigen::Array<Scalar,Rows,Cols,Options,MaxRows,MaxCols> > {
public:

    typedef Eigen::Array<Scalar,Rows,Cols,Options,MaxRows,MaxCols> Input;

    typedef boost::numpy::ndarray result_type;

    typedef typename SelectEigenView<Input>::Type View;

    static boost::numpy::ndarray apply(Input const & input) {
        return ToBoostPython< View >::apply(ndarray::copy(input));
    }

};

template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int  MaxCols>
class FromBoostPython< Eigen::Array<Scalar,Rows,Cols,Options,MaxRows,MaxCols> > {
public:

    typedef Eigen::Array<Scalar,Rows,Cols,Options,MaxRows,MaxCols> Output;

    // Use noncontiguous EigenView ('false' below) because we don't care about the input strides.
    typedef typename SelectEigenView<Output,false>::Type View;

    explicit FromBoostPython(boost::python::object const & input) : _impl(input) {}

    bool convertible() { return _impl.convertible(); }

    Output operator()() { return Output(_impl()); }

private:
    FromBoostPython< View > _impl;
};

} // namespace ndarray

#endif // !NDARRAY_EIGEN_BP_Matrix_h_INCLUDED
