// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Frank DiMaio


//  implementation loosely based on Kiss-FFT's fast fourier transform
//     for licensing info see external/kiss_fft_v1_2_8/COPYING
#if (defined WIN32) || (defined PYROSETTA)
#define _USE_MATH_DEFINES
#include <math.h>
#endif
#include <numeric/fourier/kiss_fft.hh>
#include <cstdlib> //g++ 4.3.2 requires for exit()
// tracer

#include <iostream>

#if (defined WIN32) && (!defined WIN_PYROSETTA)
#include <algorithm>
#endif

namespace numeric {
namespace fourier {

// facbuf is populated by p1,m1,p2,m2, ...
// where
// p[i] * m[i] = m[i-1]
// m0 = n
void kf_factor(int n,int * facbuf) {
	int p=4;
	double floor_sqrt = floor( sqrt((double)n) );

	// factor out powers of 4, powers of 2, then any remaining primes
	do {
		while ( n % p ) {
			switch (p) {
			case 4 : p = 2; break;
			case 2 : p = 3; break;
			default : p += 2; break;
			}
			if ( p > floor_sqrt ) {
				p = n;    // no more factors, skip to end
			}
		}
		n /= p;
		*facbuf++ = p;
		*facbuf++ = n;
	} while (n > 1);
}


/////////////////////////
/// 1D c->c fft
/////////////////////////
kiss_fft_state::kiss_fft_state():
	nfft_ (0),
	inverse_(0),
	factors_() // Zero-initialize C-style array
{}

kiss_fft_state::kiss_fft_state(int n, int inv) {
	// This used to just call resize, but has a conditional
	// that will always query uninitialized values in the ctor case!

	nfft_ = n;
	inverse_ = inv;
	twiddles_.resize(nfft_);

	for ( int i=0; i<nfft_; ++i ) {
		const double pi=M_PI;
		double phase = -2*pi*i / nfft_;
		if ( inverse_ ) {
			phase *= -1;
		}
		kf_cexp(twiddles_[i], phase );
	}
	kf_factor(nfft_,factors_);
}

void kiss_fft_state::resize(int n, int inv) {
	if ( nfft_ == n && inverse_ == inv ) return; // already the correct size

	nfft_ = n;
	inverse_ = inv;
	twiddles_.resize(nfft_);

	for ( int i=0; i<nfft_; ++i ) {
		const double pi=M_PI;
		double phase = -2*pi*i / nfft_;
		if ( inverse_ ) {
			phase *= -1;
		}
		kf_cexp(twiddles_[i], phase );
	}
	kf_factor(nfft_,factors_);
}

/////////////////////////
/// 1D c->c split fft
/////////////////////////
kiss_fftsplit_state::kiss_fftsplit_state() { }

kiss_fftsplit_state::kiss_fftsplit_state(int n, int inv) {
	resize( n, inv );
}

void kiss_fftsplit_state::resize(int n, int inv) {
	if ( substate_.nfft() == n && substate_.inverse() == inv ) return; // already the correct size

	tmpbuf_.resize( n );
	substate_.resize( n , inv );
}

/////////////////////////
/// 1D dft
/////////////////////////
kiss_dct_state::kiss_dct_state() { }

kiss_dct_state::kiss_dct_state(int n, int inv) {
	resize( n, inv );
}

void kiss_dct_state::resize(int n, int inv) {
	if ( substate_.nfft() == n && substate_.inverse() == inv ) return; // already the correct size

	if ( n & 1 ) {
		std::cerr << "Real->Real DCT optimization must be even.\n";
		exit(1);
	}
	tmpbuf_.resize( n );
	super_twiddles_.resize( n );
	substate_.resize( n , inv );

	for ( int i=0; i < n; ++i ) {
		double phase = -M_PI * ((double) i/(2*substate_.nfft()));
		if ( substate_.inverse() ) {
			phase *= -1;
		}
		kf_cexp (super_twiddles_[i],phase);
	}
}

/////////////////////////
/// 1D r->c fft
/////////////////////////
kiss_fftr_state::kiss_fftr_state() { }

kiss_fftr_state::kiss_fftr_state(int n, int inv) {
	resize( n, inv );
}

void kiss_fftr_state::resize(int n, int inv) {
	if ( substate_.nfft() == (n>>1) && substate_.inverse() == inv ) return; // already the correct size

	if ( n & 1 ) {
		std::cerr << "Real FFT optimization must be even.\n";
		exit(1);
	}
	n >>= 1;

	substate_.resize( n , inv );
	tmpbuf_.resize( n );
	super_twiddles_.resize( n/2 );

	for ( int i=0; i < n/2; ++i ) {
		double phase = -M_PI * ((double) (i+1) / substate_.nfft() + .5);
		if ( substate_.inverse() ) {
			phase *= -1;
		}
		kf_cexp (super_twiddles_[i],phase);
	}
}


/////////////////////////
/// ND c->c fft
/////////////////////////
kiss_fftnd_state::kiss_fftnd_state() { }

kiss_fftnd_state::kiss_fftnd_state(std::vector<int> const &n, int inv) {
	resize(n,inv);
}

void kiss_fftnd_state::resize(std::vector<int> const &n, int inv) {
	if ( dims_ == n && inverse_ == inv ) return; // already the correct size

	dims_ = n;
	inverse_ = inv;
	dimprod_ = 1;
	int ndims = dims_.size();

	for ( int i=0; i<ndims; ++i ) {
		dimprod_ *= dims_[i];
	}

	states_.resize( ndims );
	tmpbuf_.resize( dimprod_ );

	for ( int i=0; i<ndims; ++i ) {
		states_[i].resize(dims_[i], inverse_);
	}
}


/////////////////////////
/// ND r->c fft
/////////////////////////
kiss_fftndr_state::kiss_fftndr_state() { }

kiss_fftndr_state::kiss_fftndr_state(std::vector<int> const &n, int inv) {
	inverse_ = inv;
	resize(n,inv);
}

void kiss_fftndr_state::resize(std::vector<int> const &n, int inv) {
	// split vector
	int ndims = n.size();
	std::vector<int> nOther( ndims-1 );
	dimReal_ = n[ ndims-1 ];
	dimOther_ = 1;
	for ( int i=0; i<ndims-1; ++i ) {
		nOther[i] = n[i];
		dimOther_ *= n[i];
	}

	if ( cfg_r_.nfft() == dimReal_ && cfg_nd_.dims_v() == nOther && inverse_ == inv ) return;
	// already the correct size

	cfg_r_.resize( dimReal_, inv );
	cfg_nd_.resize( nOther, inv );
	tmpbuf_.resize(std::max( 2*dimOther_ , dimReal_+2) + dimOther_*(dimReal_+2));
}

}
}

