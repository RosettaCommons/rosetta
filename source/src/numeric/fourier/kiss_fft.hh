// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Frank DiMaio

#ifndef INCLUDED_numeric_fourier_kiss_fft_hh
#define INCLUDED_numeric_fourier_kiss_fft_hh

#include <iostream>
#include <string>
#include <limits>
#include <cmath>
#include <complex>

#include <numeric/fourier/kiss_fft_state.hh>

namespace numeric {
namespace fourier {

inline void kf_cexp(kiss_fft_cpx & x, kiss_fft_scalar phase) {
	x = kiss_fft_cpx( cos(phase) , sin(phase) );
}

////////////////////////
//  kiss_fft_alloc
//
//  Initialize a FFT (or IFFT) algorithm's cfg/state buffer.
//
//  typical usage:      kiss_fft_cfg mycfg=kiss_fft_alloc(1024,0,NULL,NULL);
//
//  The return value from fft_alloc is a cfg buffer used internally
//  by the fft routine or NULL.
//
//  If lenmem is NULL, then kiss_fft_alloc will allocate a cfg buffer using malloc.
//  The returned value should be free()d when done to avoid memory leaks.
//
//  The state can be placed in a user supplied buffer 'mem':
//  If lenmem is not NULL and mem is not NULL and *lenmem is large enough,
//      then the function places the cfg in mem and the size used in *lenmem
//      and returns mem.
//
//  If lenmem is not NULL and ( mem is NULL or *lenmem is not large enough),
//      then the function returns NULL and places the minimum cfg
//      buffer size in *lenmem.
//
//kiss_fft_cfg kiss_fft_alloc(int nfft,int inverse_fft,void * mem,size_t * lenmem);

// kiss_fft(cfg,in_out_buf)
//
// Perform an FFT on a complex input buffer.
// for a forward FFT,
// fin should be  f[0] , f[1] , ... ,f[nfft-1]
// fout will be   F[0] , F[1] , ... ,F[nfft-1]
// Note that each element is complex and can be accessed like
//    f[k].r and f[k].i
void kiss_fft(kiss_fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout);

//
// A more generic version of the above function. It reads its input from every Nth sample.
void kiss_fft_stride(kiss_fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout,int fin_stride);

// Cleans up some memory that gets managed internally. Not necessary to call, but it might clean up
// your compiler output to call this before you exit.
void kiss_fft_cleanup(void);

// Returns the smallest integer k, such that k>=n and k has only "fast" factors (2,3,5)
int kiss_fft_next_fast_size(int n);

// for real ffts, we need an even size
inline int kiss_fftr_next_fast_size_real(int n) {
	return ( kiss_fft_next_fast_size( ((n)+1)>>1)<<1 );
}

///////////////////
// Multidim FFTs
///////////////////
//kiss_fftnd_cfg  kiss_fftnd_alloc(const int *dims,int ndims,int inverse_fft,void*mem,size_t*lenmem);
void kiss_fftnd(kiss_fftnd_cfg  cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout);


///////////////////
// Real FFTs
///////////////////

// nfft must be even
// If you don't care to allocate space, use mem = lenmem = NULL
//kiss_fftr_cfg kiss_fftr_alloc(int nfft,int inverse_fft,void * mem, std::size_t * lenmem);

// input timedata has nfft scalar points
// output freqdata has nfft/2+1 complex points
void kiss_fftr(
           kiss_fftr_cfg cfg,
           const kiss_fft_scalar *timedata,
           kiss_fft_cpx *freqdata);

// input freqdata has  nfft/2+1 complex points
// output timedata has nfft scalar points
void kiss_fftri(
           kiss_fftr_cfg cfg,
           const kiss_fft_cpx *freqdata,
           kiss_fft_scalar *timedata);

///////////////////
// Multidim Real FFTs
///////////////////
// dims[0] must be even
// If you don't care to allocate space, use mem = lenmem = NULL
//kiss_fftndr_cfg  kiss_fftndr_alloc(const int *dims, int ndims, int inverse_fft, void* mem,size_t* lenmem);

// input timedata has dims[0] X dims[1] X ... X  dims[ndims-1] scalar points
// output freqdata has dims[0] X dims[1] X ... X  dims[ndims-1]/2+1 complex points
void kiss_fftndr(
        kiss_fftndr_cfg cfg,
        const kiss_fft_scalar *timedata,
        kiss_fft_cpx *freqdata);

// input and output dimensions are the exact opposite of kiss_fftndr
void kiss_fftndri(
        kiss_fftndr_cfg cfg,
        const kiss_fft_cpx *freqdata,
        kiss_fft_scalar *timedata);

}
}

#endif

