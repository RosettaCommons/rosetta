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

#ifndef INCLUDED_numeric_fourier_kiss_fft_state_hh
#define INCLUDED_numeric_fourier_kiss_fft_state_hh

#include <iostream>
#include <string>
#include <limits>
#include <cmath>
#include <complex>
#include <vector>

namespace numeric {
namespace fourier {

#define kiss_fft_scalar double
typedef std::complex< kiss_fft_scalar > kiss_fft_cpx;

/////////////////////////////////
/// In addition to state info, the state variables also contain a pointer to a precomputed buffer.
/// To avoid frequent allocation/recalculation of this buffer, it
///    may make sense to cache these buffers
/// If many FFTs of the same size were requested, only a single allocation would be needed
/// However, if many different-sized FFTs were desired, the memory usage of
///    the cache could be quite large.
#define KISSFFT_MAXFACTORS 32 ///


typedef class kiss_fft_state      * kiss_fft_cfg;
typedef class kiss_fftsplit_state * kiss_fftsplit_cfg;
typedef class kiss_dct_state      * kiss_dct_cfg;
typedef class kiss_fftnd_state    * kiss_fftnd_cfg;
typedef class kiss_fftr_state     * kiss_fftr_cfg;
typedef class kiss_fftndr_state   * kiss_fftndr_cfg;

class kiss_fft_state{
public:
	kiss_fft_state();
	kiss_fft_state(int n, int inv);

	void resize(int n, int inv);

	int nfft()    { return nfft_; }
	int inverse() { return inverse_; }

	// we give raw access to the pointers to avoid rewriting significant portions of kiss-fft
	int* factors() { return factors_; };
	kiss_fft_cpx *twiddles() { return (&twiddles_[0]); }

private:
	int nfft_;
	int inverse_;
	int factors_[2*KISSFFT_MAXFACTORS];
	std::vector<kiss_fft_cpx> twiddles_;
};

///////////////////////////
class kiss_dct_state{
public:
	kiss_dct_state();
	kiss_dct_state(int n, int inv);
	void resize(int n, int inv);

	// we give raw access to the pointers to avoid rewriting significant portions of kiss-fft
	int nfft()    { return substate_.nfft(); }
	int inverse() { return substate_.inverse(); }
	kiss_fft_cfg substate() { return &substate_; }
	kiss_fft_cpx *tmpbuf() { return (&tmpbuf_[0]); }
	kiss_fft_cpx *super_twiddles() { return (&super_twiddles_[0]); }

private:
	kiss_fft_state substate_;
	std::vector<kiss_fft_cpx> tmpbuf_;
	std::vector<kiss_fft_cpx> super_twiddles_;
};

///////////////////////////
class kiss_fftsplit_state{
public:
	kiss_fftsplit_state();
	kiss_fftsplit_state(int n, int inv);
	void resize(int n, int inv);

	// we give raw access to the pointers to avoid rewriting significant portions of kiss-fft
	int nfft()    { return substate_.nfft(); }
	int inverse() { return substate_.inverse(); }
	kiss_fft_cfg substate() { return &substate_; }
	kiss_fft_cpx *tmpbuf() { return (&tmpbuf_[0]); }

private:
	kiss_fft_state substate_;
	std::vector<kiss_fft_cpx> tmpbuf_;
};

///////////////////////////
class kiss_fftr_state{
public:
	kiss_fftr_state();
	kiss_fftr_state(int n, int inv);
	void resize(int n, int inv);

	// we give raw access to the pointers to avoid rewriting significant portions of kiss-fft
	int nfft()    { return substate_.nfft(); }
	int inverse() { return substate_.inverse(); }
	kiss_fft_cfg substate() { return &substate_; }
	kiss_fft_cpx *tmpbuf() { return (&tmpbuf_[0]); }
	kiss_fft_cpx *super_twiddles() { return (&super_twiddles_[0]); }

private:
	kiss_fft_state substate_;
	std::vector<kiss_fft_cpx> tmpbuf_;
	std::vector<kiss_fft_cpx> super_twiddles_;
};

////////////////////////////
class kiss_fftnd_state{
public:
	kiss_fftnd_state();
	kiss_fftnd_state(std::vector<int> const &n, int inv);
	void resize(std::vector<int> const &n, int inv);

	// we give raw access to the pointers to avoid rewriting significant portions of kiss-fft
	std::vector<int> dims_v() { return dims_; }
	int* dims() { return &dims_[0]; }
	int ndims() { return dims_.size(); }
	int dimprod() { return dimprod_; }
	int inverse() { return inverse_; }
	kiss_fft_cfg states(int i) { return (&states_[i]); }
	kiss_fft_cpx *tmpbuf() { return (&tmpbuf_[0]); }

private:
	int inverse_;   // for convenience store 'inv' here as well
	std::vector<int> dims_;
	int dimprod_;
	std::vector<kiss_fft_state> states_; // cfg states for each dimension
	std::vector<kiss_fft_cpx> tmpbuf_; // buffer capable of hold the entire input
};


/////////////////////////////
class kiss_fftndr_state {
public:
	kiss_fftndr_state();
	kiss_fftndr_state(std::vector<int> const &n, int inv);
	void resize(std::vector<int> const &n, int inv);

	int dimReal() { return dimReal_; }
	int dimOther() { return dimOther_; }
	int inverse() { return inverse_; }
	kiss_fftr_cfg  cfg_r() { return (&cfg_r_); }
	kiss_fftnd_cfg cfg_nd() { return (&cfg_nd_); }
	kiss_fft_scalar *tmpbuf() { return (&tmpbuf_[0]); }

private:
	int inverse_;   // for convenience store 'inv' here as well
	int dimReal_;
	int dimOther_;
	kiss_fftr_state cfg_r_;
	kiss_fftnd_state cfg_nd_;
	std::vector<kiss_fft_scalar> tmpbuf_;
};

}
}

#endif

