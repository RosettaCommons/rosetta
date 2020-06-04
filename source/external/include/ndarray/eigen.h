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
#ifndef NDARRAY_eigen_h_INCLUDED
#define NDARRAY_eigen_h_INCLUDED

/**
 *  @file ndarray/eigen.h
 *  @brief Functions that return an Eigen Map non-reference-counted view into an ndarray::Array.
 *
 *  \note This file is not included by the main "ndarray.h" header file.
 */

#if defined __GNUC__ && __GNUC__>=6
 #pragma GCC diagnostic ignored "-Wignored-attributes"
 #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "Eigen/Core"
#include "ndarray.h"
#include "ndarray/eigen_fwd.h"
#include "ndarray/EigenView.h"

namespace ndarray {
namespace detail {

template <typename T, int N, int C, typename XprKind>
struct SelectEigenPlainBase;

template <typename T, int N, int C>
struct SelectEigenPlainBase<T, N, C, Eigen::MatrixXpr> {
    typedef typename boost::remove_const<T>::type Scalar;
    typedef Eigen::Matrix<
        Scalar, Eigen::Dynamic, N == 1 ? 1 : Eigen::Dynamic,
        Eigen::AutoAlign|(C >= 0 && N > 1 ? Eigen::RowMajor : Eigen::ColMajor)
    > Type;
};

template <typename T, int N, int C>
struct SelectEigenPlainBase<T, N, C, Eigen::ArrayXpr> {
    typedef typename boost::remove_const<T>::type Scalar;
    typedef Eigen::Array<
        Scalar, Eigen::Dynamic, N == 1 ? 1 : Eigen::Dynamic,
        Eigen::AutoAlign|(C >= 0 && N > 1 ? Eigen::RowMajor : Eigen::ColMajor)
    > Type;
};

template <typename T, int N, int C, typename XprKind, bool AddConst=boost::is_const<T>::value>
struct SelectEigenPlain;

template <typename T, int N, int C, typename XprKind>
struct SelectEigenPlain<T, N, C, XprKind, true> {
    typedef typename SelectEigenPlainBase<T, N, C, XprKind>::Type const Type;
};

template <typename T, int N, int C, typename XprKind>
struct SelectEigenPlain<T, N, C, XprKind, false> {
    typedef typename SelectEigenPlainBase<T, N, C, XprKind>::Type Type;
};

template <typename T, int N, int C, typename XprKind>
struct SelectEigenMap; // unspecialized template parameters means not supported.

// 1-d, not contiguous
template <typename T, typename XprKind>
struct SelectEigenMap<T, 1, 0, XprKind> {

    typedef Eigen::Map<
        typename SelectEigenPlain<T, 1, 0, XprKind>::Type,
        Eigen::Unaligned,
        Eigen::InnerStride<>
    > Type;

    static Type apply(ndarray::Array<T, 1, 0> const & array) {
        return Type(array.getData(), array.template getSize<0>(),
                    Eigen::InnerStride<>(array.template getStride<0>()));
    }

};

// 1-d, row-major contiguous
template <typename T, typename XprKind>
struct SelectEigenMap<T, 1, 1, XprKind> {

    typedef Eigen::Map<
        typename SelectEigenPlain<T, 1, 1, XprKind>::Type,
        Eigen::Unaligned
    > Type;

    static Type apply(ndarray::Array<T, 1, 1> const & array) {
        return Type(array.getData(), array.template getSize<0>());
    }

};

// 1-d, column-major contiguous
template <typename T, typename XprKind>
struct SelectEigenMap<T, 1, -1, XprKind> {

    typedef Eigen::Map<
        typename SelectEigenPlain<T, 1, -1, XprKind>::Type,
        Eigen::Unaligned
    > Type;

    static Type apply(ndarray::Array<T, 1, -1> const & array) {
        return Type(array.getData(), array.template getSize<0>());
    }

};

// 2-d, not contiguous
template <typename T, typename XprKind>
struct SelectEigenMap<T, 2, 0, XprKind> {

    typedef Eigen::Map<
        typename SelectEigenPlain<T, 2, 0, XprKind>::Type,
        Eigen::Unaligned,
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>
    > Type;

    static Type apply(ndarray::Array<T, 2, 0> const & array) {
        return Type(array.getData(), array.template getSize<0>(), array.template getSize<1>(),
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(
                        array.template getStride<0>(), array.template getStride<1>()
                    ));
    }

};

// 2-d, row-major, contiguous within a row
template <typename T, typename XprKind>
struct SelectEigenMap<T, 2, 1, XprKind> {

    typedef Eigen::Map<
        typename SelectEigenPlain<T, 2, 1, XprKind>::Type,
        Eigen::Unaligned,
        Eigen::OuterStride<>
    > Type;

    static Type apply(ndarray::Array<T, 2, 1> const & array) {
        return Type(array.getData(), array.template getSize<0>(), array.template getSize<1>(),
                    Eigen::OuterStride<>(array.template getStride<0>()));
    }

};

// 2-d, row-major, fully contiguous
template <typename T, typename XprKind>
struct SelectEigenMap<T, 2, 2, XprKind> {

    typedef Eigen::Map<
        typename SelectEigenPlain<T, 2, 2, XprKind>::Type,
        Eigen::Unaligned
    > Type;

    static Type apply(ndarray::Array<T, 2, 2> const & array) {
        return Type(array.getData(), array.template getSize<0>(), array.template getSize<1>());
    }

};

// 2-d, column-major, contiguous within a column
template <typename T, typename XprKind>
struct SelectEigenMap<T, 2, -1, XprKind> {

    typedef Eigen::Map<
        typename SelectEigenPlain<T, 2, -1, XprKind>::Type,
        Eigen::Unaligned,
        Eigen::OuterStride<>
    > Type;

    static Type apply(ndarray::Array<T, 2, -1> const & array) {
        return Type(array.getData(), array.template getSize<0>(), array.template getSize<1>(),
                    Eigen::OuterStride<>(array.template getStride<1>()));
    }

};

// 2-d, column-major, fully contiguous
template <typename T, typename XprKind>
struct SelectEigenMap<T, 2, -2, XprKind> {

    typedef Eigen::Map<
        typename SelectEigenPlain<T, 2, -2, XprKind>::Type,
        Eigen::Unaligned
    > Type;

    static Type apply(ndarray::Array<T, 2, -2> const & array) {
        return Type(array.getData(), array.template getSize<0>(), array.template getSize<1>());
    }

};

} // namespace detail


//@{
/**
 *  Use asEigenArray(array) or asEigen<Eigen::ArrayXpr>(array)
 *  to obtain an Eigen::Array-like view of an ndarray array.
 *  Use asEigenMatrix(array) or asEigen<Eigen::MatrixXpr>(array)
 *  to obtain an Eigen::Matrix-like view of an ndarray array.
 *  The returned view is an Eigen::Map whose memory is owned by the ndarray array.
 *
 *  Be careful when calling these functions on a temporary ndarray array,
 *  as you must keep the temporary alive until you are done with the returned Eigen::Map.
 *  For example the following will have undefined behavior:
 *
 *     auto eigenMap = asEigenMatrix(makeNdArray(...));
 *     // at this point the temporary ndarray array is gone and eigenMap is broken
 *     processMatrix(eigenMap);
 *
 *  This version is safe:
 *
 *     processMatrix(asEigenMatrix(makeNdArray(...)));
 *
 *  This version is also safe and allows you to do additional processing of the Eigen map:
 *
 *     auto arr = makeNdArray(...);
 *     auto eigenMap = asEigenMatrix(arr);
 *     processMatrix(eigenMap);
 *     // ... do more with eigenMap
 */
template <typename XprKind, typename T, int N, int C>
typename detail::SelectEigenMap<T, N, C, XprKind>::Type
asEigen(Array<T, N, C> const & a) {
    return detail::SelectEigenMap<T, N, C, XprKind>::apply(a);
}

template <typename XprKind, typename T, int N, int C>
typename detail::SelectEigenMap<T, N, C, XprKind>::Type
asEigen(ArrayRef<T, N, C> const & a) {
    return detail::SelectEigenMap<T, N, C, XprKind>::apply(a);
}

template <typename T, int N, int C>
typename detail::SelectEigenMap<T, N, C, Eigen::ArrayXpr>::Type
asEigenArray(Array<T, N, C> const & a) {
    return asEigen<Eigen::ArrayXpr, T, N, C>(a);
}

template <typename T, int N, int C>
typename detail::SelectEigenMap<T, N, C, Eigen::ArrayXpr>::Type
asEigenArray(ArrayRef<T, N, C> const & a) {
    return asEigen<Eigen::ArrayXpr, T, N, C>(a);
}

template <typename T, int N, int C>
typename detail::SelectEigenMap<T, N, C, Eigen::MatrixXpr>::Type
asEigenMatrix(Array<T, N, C> const & a) {
    return asEigen<Eigen::MatrixXpr, T, N, C>(a);
}

template <typename T, int N, int C>
typename detail::SelectEigenMap<T, N, C, Eigen::MatrixXpr>::Type
asEigenMatrix(ArrayRef<T, N, C> const & a) {
    return asEigen<Eigen::MatrixXpr, T, N, C>(a);
}
//@}

} // namespace ndarray

#endif // !NDARRAY_eigen_h_INCLUDED
