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
#ifndef NDARRAY_eigen_fwd_h_INCLUDED
#define NDARRAY_eigen_fwd_h_INCLUDED

/**
 *  @file ndarray/eigen_fwd.h
 *  @brief Forward declarations for ndarray/eigen interface.
 */

/** 
 * \defgroup ndarrayEigenGroup Eigen
 * Interoperability with the Eigen 3 linear algebra library.
 */

namespace Eigen {

struct MatrixXpr;

} // namespace Eigen

namespace ndarray {
namespace detail {

template <int N, int C, int Rows, int Cols> struct EigenStrideTraits;

} // namespace detail

template <
    typename T, int N, int C, 
    typename XprKind_=Eigen::MatrixXpr,
    int Rows_=-1,
    int Cols_=((N == 1) ? 1 : -1)
    >
class EigenView;

} // namespace ndarray

#endif // !NDARRAY_eigen_fwd_h_INCLUDED
