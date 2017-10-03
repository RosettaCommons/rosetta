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
#ifndef NDARRAY_EIGEN_BP_EigenView_h_INCLUDED
#define NDARRAY_EIGEN_BP_EigenView_h_INCLUDED

#include "boost/numpy.hpp"
#include "ndarray/bp/Array.h"
#include "ndarray/eigen.h"

namespace ndarray {

template <typename T, int N, int C, typename XprKind_, int Rows_, int Cols_>
class ToBoostPython< EigenView<T,N,C,XprKind_,Rows_,Cols_> > {
public:

    typedef boost::numpy::ndarray result_type;

    static boost::numpy::ndarray apply(EigenView<T,N,C,XprKind_,Rows_,Cols_> const & input) {
        boost::numpy::ndarray result = ToBoostPython< Array<T,N,C> >::apply(input.shallow());
        if (Rows_ == 1 || Cols_ == 1) return result.squeeze();
        return result;
    }

};

template <typename T, int N, int C, typename XprKind_, int Rows_, int Cols_>
class FromBoostPython< EigenView<T,N,C,XprKind_,Rows_,Cols_> > {
public:

    BOOST_STATIC_ASSERT( N == 1 || N == 2 );

    BOOST_STATIC_ASSERT( N != 1 || Rows_ == 1 || Cols_ == 1 );

    explicit FromBoostPython(boost::python::object const & input) : _impl(input) {}

    bool convertible() {
        try {
            boost::numpy::ndarray array = boost::python::extract<boost::numpy::ndarray>(_impl.input);
            if (N == 2) {
                if (Rows_ == 1) {
                    array = array.reshape(boost::python::make_tuple(1, -1));
                } else if (Cols_ == 1) {
                    array = array.reshape(boost::python::make_tuple(-1, 1));
                }
                if (Rows_ != Eigen::Dynamic && array.shape(0) != Rows_) return false;
                if (Cols_ != Eigen::Dynamic && array.shape(1) != Cols_) return false;
            } else if (N == 1) {
                array = array.squeeze();
                int requiredSize = Rows_ * Cols_;
                if (requiredSize != Eigen::Dynamic && array.shape(0) != requiredSize) return false;
            }
            _impl.input = array;
            if (!_impl.convertible()) return false;
        } catch (boost::python::error_already_set) {
            boost::python::handle_exception();
            PyErr_Clear();
            return false;
        }
        return true;
    }

    EigenView<T,N,C,XprKind_,Rows_,Cols_> operator()() {
        return EigenView<T,N,C,XprKind_,Rows_,Cols_>(_impl());
    }

private:
    FromBoostPython< Array<T,N,C> > _impl;
};

} // namespace ndarray

#endif // !NDARRAY_EIGEN_BP_EigenView_h_INCLUDED
