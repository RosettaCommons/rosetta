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
#ifndef NDARRAY_fft_h_INCLUDED
#define NDARRAY_fft_h_INCLUDED

/**
 * @file ndarray/fft.h 
 *
 * @brief Main public header file for ndarray FFT library.
 *
 *  \note This file is not included by the main "ndarray.h" header file.
 */

#include "ndarray.h"
#include "ndarray/fft/FourierTransform.h"
#include "ndarray/fft/FourierOps.h"
#ifndef NDARRAY_FFT_MANUAL_INCLUDE
#include "ndarray/fft/FourierTransform.cc"
#endif

#endif // !NDARRAY_fft_h_INCLUDED
