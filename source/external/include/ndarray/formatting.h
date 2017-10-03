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
#ifndef NDARRAY_formatting_h_INCLUDED
#define NDARRAY_formatting_h_INCLUDED

/** 
 *  @file ndarray/formatting.h
 *
 *  @brief iostream output support for Expression.
 */

#include "ndarray/ExpressionBase.h"

#include <iostream>

namespace ndarray {
namespace detail {
template <typename Derived, int N = Derived::ND::value> class Formatter;
} // namespace detail

/**
 *  @class FormatOptions
 *  @ingroup MainGroup
 *  @brief Options for controlling stream output of ExpressionBase.
 */
class FormatOptions {
    int _width;
    int _precision;
    std::ios_base::fmtflags _flags;
    std::string _delimiter;
    std::string _open;
    std::string _close;
public:
    
    /// @brief Standard constructor.
    explicit FormatOptions(
        int width = 8,
        int precision = 6,
        std::ios_base::fmtflags flags = std::ios_base::fmtflags(0),
        std::string const & delimiter = ", ",
        std::string const & open = "[",
        std::string const & close = "]"
    ) :
        _width(width),
        _precision(precision),
        _flags(flags),
        _delimiter(delimiter),
        _open(open),
        _close(close)
    {}

    /// @brief Format the given expression into the given output stream.
    template <typename Derived>
    void apply(std::ostream & os, ExpressionBase<Derived> const & expr) {
        detail::Formatter<Derived>::apply(*this,os,expr,0);
    }

    template <typename Derived, int N> friend class detail::Formatter;
};

/// @brief Stream output for ExpressionBase using default-constructed FormatOptions.
template <typename Derived>
std::ostream & operator<<(std::ostream & os, ExpressionBase<Derived> const & expr) {
    FormatOptions options;
    options.apply(os,expr);
    return os;
}

namespace detail {

/**
 *  @internal @ingroup ndarrayInternalGroup
 *  @brief Recursive metafunction used in stream output.
 */
template <typename Derived, int N>
class Formatter {
public:
    static void apply(
        FormatOptions const & options,
        std::ostream & os, 
        ExpressionBase<Derived> const & expr,
        int level
    ) {
        os << options._open;
        if (!expr.empty()) {
            typename ExpressionBase<Derived>::Iterator const end = expr.end();
            typename ExpressionBase<Derived>::Iterator iter = expr.begin();
            Formatter<typename ExpressionBase<Derived>::Reference>::apply(options,os,*iter,level+1);
            for (++iter; iter != end; ++iter) {
                os << options._delimiter;
                os << std::endl << std::string(level,' ');
                Formatter<typename ExpressionBase<Derived>::Reference>::apply(options,os,*iter,level+1);
            }
        }
        os << options._close;
    }
};

/**
 *  @internal @ingroup ndarrayInternalGroup
 *  @brief Recursive metafunction used in stream output (1d specialization).
 */
template <typename Derived>
class Formatter<Derived,1> {
public:
    static void apply(
        FormatOptions const & options,
        std::ostream & os, 
        ExpressionBase<Derived> const & expr,
        int level
    ) {
        os << options._open;
        if (!expr.empty()) {
            typename ExpressionBase<Derived>::Iterator const end = expr.end();
            typename ExpressionBase<Derived>::Iterator iter = expr.begin();
            int precision = os.precision(options._precision);
            int width = os.width(options._width);
            std::ios_base::fmtflags flags = os.setf(options._flags,std::ios_base::floatfield);
            os << (*iter);
            for (++iter; iter != end; ++iter) {
                os << options._delimiter << (*iter);
            }
            os.precision(precision);
            os.width(width);
            os.setf(flags);
        }
        os << options._close;
    }
};

} // namespace detail
} // namespace ndarray

#endif // !NDARRAY_formatting_h_INCLUDED
