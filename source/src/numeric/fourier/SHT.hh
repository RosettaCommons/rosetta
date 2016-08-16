// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  A few spherical harmonic transform functions
/// @author Frank DiMaio

#ifndef INCLUDED_numeric_fourier_SHT_hh
#define INCLUDED_numeric_fourier_SHT_hh

#include <numeric/fourier/kiss_fft_state.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

#include <utility/vector0.hh>

#include <ObjexxFCL/FArray3D.hh>


#include <complex>

namespace numeric {
namespace fourier {

// NOTE: we use vector0's throughout since FArray's are 0-indexed using the [] operator

// helper function
void
transpose_so3(
	ObjexxFCL::FArray3D< std::complex<double> > &arrayIn,
	ObjexxFCL::FArray3D< std::complex<double> > &arrayOut,
	int m, int n );


// utility class for calculating Wigners, so3 transform
class SO3coeffs {
public:
	/// @brief Constructor
	SO3coeffs() { bw = 0; }

	/// @brief Construct with a given bandwidth, # of radial shells
	void init (int B) { bw = B; }

	// location of f_m^l in the spharm xformed data
	int lm_index(int m, int l);

	// row size in precompute_pml_trans_
	int transposeRowSize(int row, int m);

	// where in the 3D fft'd array Wigner-(m1,m2) transform starts
	int sampLoc( int m1, int m2 );

	// where in the 3D coeff array evaluated inner-products should go
	int coefLoc( int m1, int m2 );

	// Legendre xform
	// requires two 2*bw-length scratch vectors
	void Legendre(
		utility::vector0< double > &data,
		int coeffs_idx,
		int m,
		utility::vector0< double > &result,
		utility::vector0< double > &cos_pml,
		utility::vector0< double > &weights,
		utility::vector0< double > &scratch1,
		utility::vector0< double > &scratch2,
		kiss_dct_cfg dctPlan
	);

	void InvLegendre(
		utility::vector0< double > &coeffs,
		int coeffs_idx,
		int m,
		utility::vector0< double > &result,
		utility::vector0< double > &trans_cos_pml,
		utility::vector0< double > &sin_values,
		utility::vector0< double > &scratch,
		kiss_dct_cfg idctPlan
	);

	// Given orders m1, m2, and a bandwidth bw, generate all Wigner little d functions
	// scratch avoids multiple allocations with repeated calls ... each must be 2*bw in size!
	void
	genWigner_ds( int m1, int m2,
		utility::vector0< double > const &cosEval,
		utility::vector0< double > const &sinEval2,
		utility::vector0< double > const &cosEval2,
		utility::vector0< double > &result, int start_idx,
		utility::vector0< double > &scratch1,
		utility::vector0< double > &scratch2
	);

	// Wigner xform using the precomputed d's
	void wignerSynthesis(
		int m1, int m2,
		ObjexxFCL::FArray3D< std::complex<double> > &coeffs,
		int coeffs_idx,
		utility::vector0< double > &wignersTrans,
		int wigners_idx,
		ObjexxFCL::FArray3D< std::complex<double> > &signal,
		int signal_idx );

	void wignerSynthesisSameSign(
		int m1, int m2,
		ObjexxFCL::FArray3D< std::complex<double> > &coeffs,
		int coeffs_idx,
		utility::vector0< double > &wignersTrans,
		int wigners_idx,
		ObjexxFCL::FArray3D< std::complex<double> > &signal,
		int signal_idx  );

	void wignerSynthesisDiffSign(
		int m1, int m2,
		ObjexxFCL::FArray3D< std::complex<double> > &coeffs,
		int coeffs_idx,
		utility::vector0< double > &wignersTrans,
		int wigners_idx,
		ObjexxFCL::FArray3D< std::complex<double> > &signal,
		int signal_idx,
		ObjexxFCL::FArray3D< std::complex<double> > &scratch );

private:
	int bw;
};


class SHT {
public:
	/// @brief Constructor
	SHT();

	/// @brief Construct with a given bandwidth, # of radial shells
	SHT(int B, int nR);

	/// @brief Destructor
	~SHT();

	/// @brief Initialize with a given bandwidth, # of radial shells
	void init(int B, int nR);


	/// @brief Take spherical harmonic transform of 'sigR', result in 'sigCoefR' and 'sigCoefI'
	void sharm_transform( ObjexxFCL::FArray3D< double > const & sigR,
		ObjexxFCL::FArray3D< double > & sigCoefR,
		ObjexxFCL::FArray3D< double > & sigCoefI);

	/// @brief Take inverse spherical harmonic transform of 'sigCoefR' and 'sigCoefI', result in 'sigR'
	void sharm_invTransform( ObjexxFCL::FArray3D< double > & sigR,
		ObjexxFCL::FArray3D< double > & sigCoefR,  // should be const
		ObjexxFCL::FArray3D< double > & sigCoefI); // should be const

	/// @brief Standardize coefficients 'sigCoefR' and 'sigCoefI'
	void sph_standardize( ObjexxFCL::FArray3D< double > & sigCoefR,
		ObjexxFCL::FArray3D< double > & sigCoefI);

	/// @brief Correlate two signals ('sigCoef' and 'tmpCoef') as a function of rotation of 'tmpCoef', result in 'so3_correlation'
	void so3_correlate( ObjexxFCL::FArray3D< double > & so3_correlation,
		ObjexxFCL::FArray3D< double > & sigCoefR,  // should be const
		ObjexxFCL::FArray3D< double > & sigCoefI,  // should be const
		ObjexxFCL::FArray3D< double > & tmpCoefR,  // should be const
		ObjexxFCL::FArray3D< double > & tmpCoefI  ); // should be const

	/// @brief Convert an index from 'so3_correlation' into a rotation matrix
	void idx_to_rot(int maxloc , numeric::xyzMatrix< double > & thisRot);

	/// @brief Convert an index from 'so3_correlation'  into Euler angles alpha, beta, gamma
	void idx_to_euler(int maxloc , numeric::xyzVector< double > & eulerAngles);

private:

	// helper functions
	// take raw pointers since we are accessing 2D slices of 3D data and do not want to reallocate!

	// @brief in one shell, combine two S^2 spherical coefficients so
	//     the inverse SO(3) transform will result in their rotational correl
	void
	so3CombineCoef(
		utility::vector0< double > &sigCoefR,
		utility::vector0< double > &sigCoefI,
		utility::vector0< double > &tmpCoefR,
		utility::vector0< double > &tmpCoefI,
		ObjexxFCL::FArray3D< std::complex<double> > &so3Coef);

	// @brief do the inverse SO(3) transform
	//  read from so3Coef
	//  write to conv
	void
	inverseSo3( ObjexxFCL::FArray3D< std::complex<double> > &so3Coef, ObjexxFCL::FArray3D< std::complex<double> > &conv );

	// @brief forward transform in one shell
	void
	forwardS2(
		utility::vector0< double > &rdata,
		utility::vector0< double > &rcoeffs,
		utility::vector0< double > &icoeffs );

	// @brief inverse transform in one shell
	void
	inverseS2(
		utility::vector0< double > &rcoeffs,
		utility::vector0< double > &icoeffs,
		utility::vector0< double > &rdata );

	// @brief setup Pml tables
	void
	setup_Pmls( );

	// @brief setup Wigner table
	void
	setup_Wig( );

	// @brief setup weights
	void
	setup_Weights( );


private:
	int bw, nRsteps;
	SO3coeffs so3_;

	// storage space for ffts
	kiss_fft_cfg p1;
	kiss_fftsplit_cfg fftPlan, ifftPlan;
	kiss_dct_cfg dctPlan, idctPlan;

	// precomputed data
	utility::vector0< double > weights_, precompute_wig_trans_;
	utility::vector0< utility::vector0<double> > precompute_pml_;
	utility::vector0< utility::vector0<double> > precompute_pml_trans_;
};

}
}

#endif

