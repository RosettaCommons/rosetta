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
#ifndef NDARRAY_views_h_INCLUDED
#define NDARRAY_views_h_INCLUDED

/** 
 *  \file ndarray/views.h @brief Public interface for arbitrary views into arrays.
 */

#include <boost/fusion/include/push_back.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/include/mpl.hpp>

#include "ndarray_fwd.h"

namespace ndarray {
namespace index {

/** 
 *  @brief Simple structure defining a noncontiguous range of indices.
 */
struct Slice {
    Size start;
    Size stop;
    Offset step;

    Slice(Size start_, Size stop_, Offset step_) : start(start_), stop(stop_), step(step_) {}

    Size computeSize() const { return (step > 1) ? (stop - start + 1) / step : stop - start; }
};

/**
 *  @brief Simple structure defining a contiguous range of indices.
 */
struct Range {
    Size start;
    Size stop;

    Range(Size start_, Size stop_) : start(start_), stop(stop_) {}
};

/**
 *  @brief Empty structure marking a view of an entire dimension.
 */
struct Full {};

/**
 *  @brief Structure marking a single element of a dimension.
 */
struct Scalar {
    Size n;

    explicit Scalar(Size n_) : n(n_) {}
};

} // namespace index

/** 
 *  @brief A template meta-sequence that defines an arbitrary view into an unspecified array. 
 *
 *  A View is constructed from a call to the global view() function
 *  and subsequent chained calls to operator().
 */
template <typename Seq_ = boost::fusion::vector<> >
struct View {
    typedef Seq_ Sequence; ///< A boost::fusion sequence type
    Sequence _seq; ///< A boost::fusion sequence of index objects.

    explicit View(Sequence seq) : _seq(seq) {}

    template <typename OtherSequence>
    explicit View(OtherSequence const & other) : _seq(other) {}

    template <typename OtherSequence>
    View(View<OtherSequence> const & other) : _seq(other._seq) {}

    /// @brief The View that results from chaining an full dimension index <b><tt>()</tt></b> to this.
    typedef View<typename boost::fusion::result_of::push_back<Sequence const,index::Full>::type> Full;

    /// @brief The View that results from chaining a range <b><tt>(start,stop)</tt></b> to this.
    typedef View<typename boost::fusion::result_of::push_back<Sequence const,index::Range>::type> Range;

    /// @brief The View that results from chaining a slice <b><tt>(start,stop,step)</tt></b> to this.
    typedef View<typename boost::fusion::result_of::push_back<Sequence const,index::Slice>::type> Slice;

    /// @brief The View that results from chaining a scalar <b><tt>(n)</tt></b> to this.
    typedef View<typename boost::fusion::result_of::push_back<Sequence const,index::Scalar>::type> Scalar;

    /// @brief Chain the full next dimension to this.
    Full operator()() const { return Full(boost::fusion::push_back(_seq, index::Full())); }
    
    /// @brief Chain a contiguous range of the next dimension to this.
    Range operator()(Size start, Size stop) const {
        return Range(boost::fusion::push_back(_seq, index::Range(start, stop)));
    }

    /// @brief Chain a noncontiguous slice of the next dimension to this.
    Slice operator()(Size start, Size stop, Offset step) const {
        return Slice(boost::fusion::push_back(_seq, index::Slice(start, stop, step)));
    }

    /// @brief Chain a single element of the next dimension to this.
    Scalar operator()(Size n) const {
        return Scalar(boost::fusion::push_back(_seq, index::Scalar(n)));
    }
};

/// @addtogroup ndarrayMainGroup
/// @{

/** @brief Start a view definition that includes the entire first dimension. */
inline View< boost::fusion::vector1<index::Full> > view() {
    return View< boost::fusion::vector1<index::Full> >(
        boost::fusion::make_vector(index::Full())
    );
}

/** @brief Start a view definition that selects a contiguous range in the first dimension. */
inline View< boost::fusion::vector1<index::Range> > view(Size start, Size stop) {
    return View< boost::fusion::vector1<index::Range> >(
        boost::fusion::make_vector(index::Range(start, stop))
    );
}

/** @brief Start a view definition that selects a noncontiguous slice of the first dimension. */
inline View< boost::fusion::vector1<index::Slice> > view(Size start, Size stop, Offset step) {
    return View< boost::fusion::vector1<index::Slice> >(
        boost::fusion::make_vector(index::Slice(start, stop, step))
    );
}

/** @brief Start a view definition that selects single element from the first dimension. */
inline View< boost::fusion::vector1<index::Scalar> > view(Size n) {
    return View< boost::fusion::vector1<index::Scalar> >(
        boost::fusion::make_vector(index::Scalar(n))
    );
}

/// @}

} // namespace ndarray

#endif // !NDARRAY_views_h_INCLUDED
