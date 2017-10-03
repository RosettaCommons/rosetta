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
#ifndef NDARRAY_FFT_fft_fwd_h_INCLUDED
#define NDARRAY_FFT_fft_fwd_h_INCLUDED

/**
 * @file ndarray/fft_fwd.h 
 *
 * @brief Forward declarations and default template parameters for ndarray/fft.
 *
 *  \note This file is not included by the main "ndarray.h" header file.
 */

/**
 *  \defgroup FFTGroup Fourier Transforms
 *
 *  @brief Fast fourier transforms using the FFTW library.
 */

/// @internal \defgroup FFTndarrayInternalGroup Fourier Transform Internals

#include "ndarray_fwd.h"

namespace ndarray {
namespace detail {

template <typename T, bool IsConst=boost::is_const<T>::value> struct FourierTraits;
template <typename T> struct FFTWTraits;

} // namespace detail

template <typename T, int N> class FourierTransform;

} // namespace ndarray

#endif // !NDARRAY_FFT_fft_fwd_h_INCLUDED
