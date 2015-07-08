// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  A few spherical harmonic transform function: very simple, with support for real, even BW only
/// @author Frank DiMaio


#include <numeric/fourier/SHT.hh>
#include <numeric/constants.hh>


#include <string>
#include <cstdlib>

#ifdef WIN32
	#undef min
	#undef max
	#undef PostMessage
#endif

#include <algorithm>

namespace numeric {
namespace fourier {


// Workaround for MSVC bug
#ifdef WIN32
	template <class T> __declspec(noinline) const T& local_max (const T& a, const T& b) { return (a<b)?b:a; }
#else
	template <class T> const T& local_max (const T& a, const T& b) { return (a<b)?b:a; }
#endif


void
transpose_so3(
		ObjexxFCL::FArray3D< std::complex<double> > &arrayIn,
		ObjexxFCL::FArray3D< std::complex<double> > &arrayOut,
		int m, int n ) {
	int mud;

	if( n >= m ) {
		for( int d = n-m+1 ; d < n ; d ++ ) {
			mud = n - d ;
			for( int k = 0 ; k < mud ; k ++ )
				arrayOut[ m*(k+d) + k ] = arrayIn[ n*k + k + d ];
		}
		for( int d = 0 ; d < n-m+1 ; d ++ ) {
			mud = m ;
			for( int k = 0 ; k < mud ; k ++ )
				arrayOut[ m*(k+d) + k ] = arrayIn[ n*k + k + d] ;
		}
		for( int d = -m+1 ; d < 0 ; d ++ ) {
			mud = m + d ;
			for( int k = 0 ; k < mud ; k ++ )
				arrayOut[ m*k + k - d ] = arrayIn[ n*(k-d) + k ] ;
		}
	} else {
		for( int d = 0 ; d < n ; d ++ ) {
			mud = n - d ;
			for( int k = 0 ; k < mud ; k ++ )
				arrayOut[ m*(k+d) + k ] = arrayIn[ n*k + k + d ] ;
		}
		for( int d = n-m ; d < 0 ; d ++ ) {
			mud = n ;
			for( int k = 0 ; k < mud ; k ++ )
				arrayOut[ m*(k) + k - d ] = arrayIn[ n*(k-d) + k ] ;
		}
		for( int d = -m+1 ; d < n-m ; d ++ ) {
			mud = m + d ;
			for( int k = 0 ; k < mud ; k ++ )
				arrayOut[ m*k + k - d ] = arrayIn[ n*(k-d) + k ] ;
		}
	}
}


////////////////////
////////////////////

// lm->index
int
SO3coeffs::lm_index(int m, int l) {
	if( m >= 0 )
		return( m * bw - ( ( m * (m - 1) ) /2 ) + ( l - m ) );
	else
		return( ( ((bw-1)*(bw+2))/2 ) + 1 + ( (bw-1+m)*(bw+m)/2 ) + ( l-std::abs( m ) ) );
}

// Transposed row size returns the number of non-zero coefficients
//   in the transposition of 'precompute_pml_'
int
SO3coeffs::transposeRowSize( int row, int m ) {
	if ( bw % 2 ) {
		if ( m % 2 ) {
			if ( m == 1 )
				return( (bw-row)/2 );
			else if ( row < m - 1 )
				return ( (bw-m+1)/2 );
			else
				return ( transposeRowSize(row, 1) ) ;
		} else {
			if ( m == 0 )
				return( (bw-row)/2 + ((row+1)%2) );
			else if ( row < m )
				return ( (bw-m)/2 + ((row+1)%2) );
			else
				return ( transposeRowSize(row, 0) ) ;
		}
	} else {
		if ( m % 2 ) {
			if ( m == 1 )
				return( (bw-row)/2 );
			else if ( row < m - 1 )
				return ( (bw-m+1)/2 - (row%2) );
			else
				return ( transposeRowSize(row, 1) ) ;
		} else {
			if ( m == 0 )
				return( (bw-row)/2 + (row%2) );
			else if ( row < m )
				return ( (bw-m)/2 );
			else
				return ( transposeRowSize(row, 0) ) ;
		}
	}
}

int
SO3coeffs::sampLoc( int m1, int m2 ) {
	if ( m1 >= 0 ) {
		if ( m2 >= 0 )
			return ( (2*bw)*((2*bw*m1 + m2)) );
		else
			return ( (2*bw)*((2*bw*m1 + (2*bw+m2))) );
	} else {
		if ( m2 >= 0 )
			return ( (2*bw)*((2*bw)*(2*bw+m1) + m2) );
		else
			return ( (2*bw)*((2*bw)*(2*bw+m1) + (2*bw+m2)) );
	}
}


int
SO3coeffs::coefLoc( int m1, int m2 ) {
	int retval=0;

	if ( m1 >= 0 ) {
		if ( m2 >= 0 ) {
			retval = (bw*bw)*m1 - ( (m1-1)*m1*(2*m1-1)/6 ) ;
			for(int k = 0 ; k < m2 ; k ++ ) { retval += bw - std::max(m1,k); }
		} else {
			retval = (bw*bw)*(m1+1) - ( (m1+1)*m1*(2*m1+1)/6 );
			for(int k = m2 ; k < 0 ; k ++ ) { retval -= bw - local_max(m1,-k); }
		}
	} else {
		if ( m2 >= 0 ) {
			retval = (4*bw*bw*bw-bw)/3 - (bw*bw)*(-m1) + ( (-m1+1)*(-m1)*(2*(-m1)+1)/6 ) ;
			for (int k = 0 ; k < m2 ; k ++ ) { retval += bw - std::max(-m1,k); }
		} else {
			retval = (4*bw*bw*bw-bw)/3 - (bw*bw)*(-m1-1) + ( (-m1-1)*(-m1)*(2*(-m1)-1)/6 );
			for (int k = m2 ; k < 0 ; k ++ ) { retval -= bw - std::max(-m1,-k); }
		}
	}

	return retval;
}


//    Given orders m1, m2, and a bandwidth bw,
//      generate all Wigner little d functions
void
SO3coeffs::genWigner_ds( int m1, int m2,
		utility::vector0< double > const &cosEval,
		utility::vector0< double > const &sinEval2,
		utility::vector0< double > const &cosEval2,
		utility::vector0< double > &result, int start_idx,
		utility::vector0< double > &scratch1,
		utility::vector0< double > &scratch2 ) {
	int m = std::max( std::abs( m1 ), std::abs( m2 ) );
	int N = 2*bw ;

	//  initialize arrays for the recurrence
	for (int i=0; i<N; ++i) scratch1[i]=0;

	int cosExp, sinExp;
	int l = std::max( std::abs( m1 ), std::abs( m2 ) ) ;
	int delta = l - std::min( std::abs( m1 ), std::abs( m2 ) ) ;
	double sinSign = 1.0;
	double normFactor = sqrt((2.0*l+1.0)/2.0);
	for ( int i = 0 ; i < delta ; i ++ ) normFactor *= sqrt( (2.0*l - i)/( i + 1.0 ) );

	if ( l == std::abs( m1 ) ) {
		if ( m1 >= 0 ) {
			cosExp = l + m2 ; sinExp = l - m2 ;
			if ( (l - m2) % 2 ) sinSign = -1.0;
		} else {
			cosExp = l - m2 ; sinExp = l + m2 ;
		}
	} else if ( m2 >= 0 ) {
		cosExp = l + m1 ; sinExp = l - m1 ;
	} else {
		cosExp = l - m1 ; sinExp = l + m1 ;
		if ( (l + m1) % 2 ) sinSign = -1.0 ;
	}

	for ( int i = 0 ; i < N ; i++ ) {
		scratch2[i] = normFactor * sinSign * pow( sinEval2[i], (double)sinExp ) * pow( cosEval2[i], (double)cosExp ) ;
		result[ start_idx + (bw-m) * i ] = scratch2[ i ] ;
	}

	for ( int i = 0 ; i < bw-m-1; i++ ) {
		double a1 = sqrt( (2.0*(m+i) + 3.0)/(2.0*(m+i) - 1.0) ) * ((m+i) + 1.0);
		double b1 = sqrt( (2.0*(m+i) + 3.0)/(2.0*(m+i) + 1.0) ) * ((m+i) + 1.0);
		double a2 = sqrt( ((m+i)*(m+i) - m1*m1)*((m+i)*(m+i) - m2*m2) );
		double ab3 = sqrt((((m+i)+1.0)*((m+i)+1.0) - m1*m1) * (((m+i)+1.0)*((m+i)+1.0) - m2*m2));
		double a = (m+i<=0) ? 0.0 : -a1 / ab3 * a2 / (m+i);
		double b = b1*(2.0*(m+i) + 1.0) / ab3;
		double c = (m+i==0) ? 0.0 : -b * m1 * m2 / ( (m+i) * ( (m+i) + 1.0) );

		for ( int j = 0; j<N; ++j) {
			double result_j = c*scratch2[j]+a*scratch1[j]+cosEval[j]*b*scratch2[j];
			result[ start_idx + i+1+(bw-m)*j ] = result_j;
			scratch1[j] = scratch2[j];
			scratch2[j] = result_j;
		}
	}
}


// Computes the Legendre transform
void
SO3coeffs::Legendre(
		utility::vector0< double > &data,
		int data_idx,
		int m,
		utility::vector0< double > &result,
		utility::vector0< double > &cos_pml,
		utility::vector0< double > &weights,
		utility::vector0< double > &scratch1,
		utility::vector0< double > &scratch2,
		kiss_dct_cfg dctPlan ) {
	int n = 2*bw;

	for (int i=0; i<(bw-m); ++i) result[i]=0;

	// apply weights and compute cosine transform
	if ( m % 2 )
		for ( int i = 0; i < n ; ++i )
			scratch1[i] = data[ data_idx+i ] * weights[ 2*bw + i ];
	else
		for ( int i = 0; i < n ; ++i )
			scratch1[i] = data[ data_idx+i ] * weights[ i ];
	kiss_dct( dctPlan, &scratch1[0], &scratch2[0]);

	// normalize
	scratch2[0] *= 1.0/sqrt(2.0) ;
	double cos_scale = 1.0/sqrt(2.0*n);
	for ( int j = 0 ; j < n ; j ++ ) scratch2[j] *= cos_scale;

	//  project
	int toggle = 0 ;
	for ( int i=m; i<bw; i++ ) {
		int pml_index;
		if (m%2==0) {
			pml_index = ((i/2)*((i/2)+1)) - ((m/2)*((m/2)+1));
			if (i%2==1) pml_index += (i/2)+1;
		} else {
			pml_index = (((i-1)/2)*(((i-1)/2)+1)) - (((m-1)/2)*(((m-1)/2)+1));
			if (i%2==0) pml_index += ((i-1)/2)+1;
		}

		for ( int j=0; j < (i/2) ; j ++ )
			result[i-m] += scratch2[(2*j)+toggle] * cos_pml[pml_index+j];

		if ( (((i-m)%2) == 0 ) || ( (m%2) == 0 ))
			result[i-m] += scratch2[(2*(i/2))+toggle] * cos_pml[pml_index+(i/2)];

		toggle = (toggle + 1)%2 ;
	}
}


//  Computes the inverse Legendre transform
void
SO3coeffs::InvLegendre(
		utility::vector0< double > &coeffs,
		int coeffs_idx,
		int m,
		utility::vector0< double > &result,
		utility::vector0< double > &trans_cos_pml,
		utility::vector0< double > &sin_values,
		utility::vector0< double > &scratch,
		kiss_dct_cfg idctPlan ) {

	for (int i=0; i<2*bw; ++i) {
		result[i]=0.0;
	}

	int trans_pml_idx=0, coeffs_offset=0;
	for (int i=0; i<bw; i++) {
		if (i == (bw-1) && m % 2 ) {
			scratch[bw-1] = 0.0;
			break;
		}

		if (i > m)
			coeffs_offset = (i-m) + (m % 2);
		else
			coeffs_offset = (i%2);

		int rowsize = transposeRowSize(i, m);
		for (int j = 0; j < rowsize; ++j)
			scratch[i] += coeffs[coeffs_idx + coeffs_offset + 2*j] * trans_cos_pml[trans_pml_idx+j];

		trans_pml_idx += rowsize;
	}

	double cos_scale = 0.5 / sqrt((double) bw) ;
	for ( int j = 1 ; j < 2*bw ; j ++ ) scratch[j] *= cos_scale ;
	scratch[0] /= sqrt(2.0 * ((double) bw));

	kiss_dct( idctPlan, &scratch[0], &result[0]);

	if ( m % 2 ) {
		for ( int j=0; j<(2*bw); j++)
			result[j] *= sin_values[j];
	}
}

//  synthesize the signal of length n = 2*bw
void
SO3coeffs::wignerSynthesis(
		int m1, int m2,
		ObjexxFCL::FArray3D< std::complex<double> > &coeffs,
		int coeffs_idx,
		utility::vector0< double > &wignersTrans,
		int wigners_idx,
		ObjexxFCL::FArray3D< std::complex<double> > &signal,
		int signal_idx ) {
	int wigners_curr_idx = wigners_idx;

	int m = std::max( std::abs( m1 ) , std::abs( m2 ) ) ;
	int n = 2 * bw ;

	for ( int i = 0 ; i < n ; i ++ ) {
		double realPart = 0.0 ;
		double imagPart = 0.0 ;
		for ( int j = 0 ; j < (bw - m) ; j ++ ) {
			realPart += wignersTrans[wigners_curr_idx] * coeffs[coeffs_idx+j].real();
			imagPart += wignersTrans[wigners_curr_idx] * coeffs[coeffs_idx+j].imag();
			wigners_curr_idx++ ;
		}
		signal[ signal_idx+i ] = std::complex<double>( realPart , imagPart );
	}
}

void
SO3coeffs::wignerSynthesisSameSign(
		int m1, int m2,
		ObjexxFCL::FArray3D< std::complex<double> > &coeffs,
		int coeffs_idx,
		utility::vector0< double > &wignersTrans,
		int wigners_idx,
		ObjexxFCL::FArray3D< std::complex<double> > &signal,
		int signal_idx ) {
	int coeffs_scale=1;

	int m = std::max( std::abs( m1 ) , std::abs( m2 ) ) ;
	int n = 2 * bw ;

	if ( std::abs( m1 - m2 ) % 2 )
		coeffs_scale = -1 ;

	int wigners_curr_idx = wigners_idx;

	for ( int i = 0 ; i < n ; i ++ ) {
		double realPart = 0.0 ;
		double imagPart = 0.0 ;
		for ( int j = 0 ; j < (bw - m) ; j ++ ) {
			realPart += wignersTrans[wigners_curr_idx] * coeffs[coeffs_idx+j].real();
			imagPart += wignersTrans[wigners_curr_idx] * coeffs[coeffs_idx+j].imag();
			wigners_curr_idx++ ;
		}
		signal[ signal_idx+i ] = std::complex<double>( coeffs_scale * realPart , coeffs_scale * imagPart );
	}
}

void
SO3coeffs::wignerSynthesisDiffSign(
		int m1, int m2,
		ObjexxFCL::FArray3D< std::complex<double> > &coeffs,
		int coeffs_idx,
		utility::vector0< double > &wignersTrans,
		int wigners_idx,
		ObjexxFCL::FArray3D< std::complex<double> > &signal,
		int signal_idx,
		ObjexxFCL::FArray3D< std::complex<double> > &scratch ) {
	int coeffs_scale=1;

	int m = std::max( std::abs( m1 ) , std::abs( m2 ) ) ;
	int n = 2 * bw ;

	if ( (m1 < 0 && (m - m2) % 2) || (m1>=0 && (m + m1) % 2 ) )
		coeffs_scale = -1 ;

	// place the ptr at the end of the array
	int wigners_curr_idx = wigners_idx + (bw - m) * n - (bw-m);

	for ( int i=0 ; i < n ; i ++ ) {
		scratch[i] = ((double)coeffs_scale) * coeffs[coeffs_idx+i];
		coeffs_scale *= -1 ;
	}

	for ( int i = 0 ; i < n ; i ++ ) {
		double realPart = 0.0 ;
		double imagPart = 0.0 ;
		for ( int j = 0 ; j < ( bw - m ) ; j ++ ) {
			realPart += wignersTrans[wigners_curr_idx] * scratch[j].real();
			imagPart += wignersTrans[wigners_curr_idx] * scratch[j].imag();
			wigners_curr_idx++ ;
		}
		signal[ signal_idx+i ] = std::complex<double>( realPart , imagPart );
		wigners_curr_idx = wigners_idx + (n - i - 1)*(bw-m) - (bw - m);
	}
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////


// null ctor
SHT::SHT() {
	bw = nRsteps = 0;
	p1 = NULL;
	fftPlan = ifftPlan = NULL;
	dctPlan = idctPlan = NULL;
}

// initialize for a given bandwidth, # of radial shells
SHT::SHT(int B, int nR) {
	bw = nRsteps = 0;
	p1 = NULL;
	fftPlan = ifftPlan = NULL;
	dctPlan = idctPlan = NULL;

	init(B, nR);
}

// destructor
SHT::~SHT() {
	if (p1) delete p1;
	if (fftPlan) delete fftPlan;
	if (ifftPlan) delete ifftPlan;
	if (dctPlan) delete dctPlan;
	if (idctPlan) delete idctPlan;
}

// initialize
void
SHT::init(int B, int nR) {
	if (B == bw && nR == nRsteps) {
		// no reinitialization needed
		return;
	}

	if (B%2 == 1) {
		std::cerr << "Odd bandwidths unsupported!\n" << std::endl;
		exit(1);
	}

	bw = B;
	nRsteps = nR;
	so3_.init(B);

	precompute_pml_.clear();
	precompute_pml_trans_.clear();
	precompute_wig_trans_.clear();

	// fft plans
	if (p1) delete p1;
	if (fftPlan) delete fftPlan;
	if (ifftPlan) delete ifftPlan;
	if (dctPlan) delete dctPlan;
	if (idctPlan) delete idctPlan;

	// precomputed coeffs (keep this order!)
	setup_Weights( );
	setup_Pmls( );
	setup_Wig( );

	// fft scratch space
	p1 = new kiss_fft_state(2*bw,0);
	fftPlan = new kiss_fftsplit_state (2*bw,0);
	ifftPlan = new kiss_fftsplit_state (2*bw,0);
	dctPlan = new kiss_dct_state (2*bw,0);
	idctPlan = new kiss_dct_state (2*bw,1);
}


// spherical harmonic transform
void SHT::sharm_transform(
		ObjexxFCL::FArray3D< double > const & sigR,
		ObjexxFCL::FArray3D< double > & sigCoefR,
		ObjexxFCL::FArray3D< double > & sigCoefI
) {
	int R = sigR.u3(), B2=sigR.u2();
	init (B2/2,R);

	sigCoefR.dimension(bw,bw,nRsteps);
	sigCoefI.dimension(bw,bw,nRsteps);

	// 2D slices.  this allocation may be inefficient
	utility::vector0< double > sigRslice( 2*bw*2*bw );
	utility::vector0< double > sigCoefRslice( bw*bw );
	utility::vector0< double > sigCoefIslice( bw*bw );

	for (int r_idx=0; r_idx<nRsteps; ++r_idx) {
		for (int i=0; i<4*bw*bw; ++i)
			sigRslice[i] = sigR[r_idx*4*bw*bw + i];

		forwardS2( sigRslice, sigCoefRslice, sigCoefIslice);

		for (int i=0; i<bw*bw; ++i) {
			sigCoefR[r_idx*4*bw*bw + i] = sigCoefRslice[i];
			sigCoefI[r_idx*4*bw*bw + i] = sigCoefIslice[i];
		}
	}
}


// spherical harmonic inverse transform
void SHT::sharm_invTransform(
		ObjexxFCL::FArray3D< double > & sigR,
		ObjexxFCL::FArray3D< double > & sigCoefR,
		ObjexxFCL::FArray3D< double > & sigCoefI
) {
	int R = sigCoefR.u3(), B=sigCoefR.u2();
	init (B,R);

	sigR.dimension(2*bw,2*bw,nRsteps);

	// 2D slices.  this allocation may be inefficient
	utility::vector0< double > sigRslice( 2*bw*2*bw );
	utility::vector0< double > sigCoefRslice( bw*bw );
	utility::vector0< double > sigCoefIslice( bw*bw );

	for (int r_idx=0; r_idx<nRsteps; ++r_idx) {
		for (int i=0; i<bw*bw; ++i) {
			sigCoefRslice[i] = sigCoefR[r_idx*4*bw*bw + i];
			sigCoefIslice[i] = sigCoefI[r_idx*4*bw*bw + i];
		}

		inverseS2 ( sigCoefRslice, sigCoefIslice, sigRslice );

		for (int i=0; i<4*bw*bw; ++i)
			sigR[r_idx*4*bw*bw + i] = sigRslice[i];
	}
}

// standardize the signal given spharm coefficients
// wangyr: report s and s2 in order to debug
void SHT::sph_standardize(
	ObjexxFCL::FArray3D< double > & sigCoefR,
	ObjexxFCL::FArray3D< double > & sigCoefI,
	double & s,
	double & s2
)
{
	// don't care if init() has been called
	int R = sigCoefR.u3(), B=sigCoefR.u2();

	double b_scaleFactor = sqrt(1.0/sqrt(0.015*B)); // this normalizes map s.t. correlation with itself ~ 1

	double wt=0;
	s = s2 = 0;

	for (int r_idx=1; r_idx<=R; ++r_idx) {
		// standardize
		double thisRWt = numeric::constants::d::pi * ( 4.0 * square((double)r_idx) + 1.0/3.0 );

		double sumCoef2 = 0.0;
		for (int i=1; i<=B; ++i)
		for (int j=1; j<=B; ++j)
			if (i != 1 || j!= 1)
				sumCoef2 += square(sigCoefR(j,i,r_idx)) + square(sigCoefI(j,i,r_idx));

		wt += thisRWt;
		s  += thisRWt * sigCoefR(1,1,r_idx);
		s2 += thisRWt * sumCoef2;
	}

	s = s/wt;
	s2 = sqrt(2 * b_scaleFactor * wt / s2);

	//std::cout << "Standardizing s: " << s << " , s2: " << s2 << std::endl;
	for (int r_idx=1; r_idx<=R; ++r_idx) {
		for (int i=1; i<=B; ++i)
		for (int j=1; j<=B; ++j) {
			if (i == 1 && j== 1) {
				// adjust mean
				sigCoefR(1,1,r_idx) = -s + sigCoefR(1,1,r_idx);
				sigCoefI(1,1,r_idx) = 0.0 ;
			} else {
				// adjust stdev
				sigCoefR(j,i,r_idx) = s2 * sigCoefR(j,i,r_idx) ;
				sigCoefI(j,i,r_idx) = s2 * sigCoefI(j,i,r_idx) ;
			}
		}
	}
}

// correlate two signals as a function of rotation
void SHT::so3_correlate(
		ObjexxFCL::FArray3D< double > & so3_correlation,
		ObjexxFCL::FArray3D< double > & sigCoefR,  // should be const
		ObjexxFCL::FArray3D< double > & sigCoefI,  // should be const
		ObjexxFCL::FArray3D< double > & tmpCoefR,  // should be const
		ObjexxFCL::FArray3D< double > & tmpCoefI
) {
	int R = sigCoefR.u3(), B=sigCoefR.u2();
	init (B,R);

	double sumRWt = 0.0;
	for (int r_idx=1; r_idx<=nRsteps; ++r_idx)
		sumRWt += numeric::constants::d::pi * ( 4.0 * square((double)r_idx) + 1.0/3.0 );

	ObjexxFCL::FArray3D< std::complex<double> > so3Sig, so3Coef;
	so3Sig.dimension(2*bw,2*bw,2*bw);
	so3Coef.dimension(2*bw,2*bw,2*bw);
	so3_correlation.dimension(2*bw,2*bw,2*bw);

	for (int i=0; i<8*bw*bw*bw; ++i) so3_correlation[i] = 0.0;

	// 2D slices.  we could save some time by not reallocating ...
	utility::vector0< double > sigCoefR_slice(bw*bw);
	utility::vector0< double > sigCoefI_slice(bw*bw);
	utility::vector0< double > tmpCoefR_slice(bw*bw);
	utility::vector0< double > tmpCoefI_slice(bw*bw);

	for (int r_idx=0; r_idx<nRsteps; ++r_idx) {
		for (int i=0; i<bw*bw; ++i) {
			sigCoefR_slice[i] = sigCoefR[r_idx*bw*bw+i];
			sigCoefI_slice[i] = sigCoefI[r_idx*bw*bw+i];
			tmpCoefR_slice[i] = tmpCoefR[r_idx*bw*bw+i];
			tmpCoefI_slice[i] = tmpCoefI[r_idx*bw*bw+i];
		}

		so3CombineCoef( sigCoefR_slice, sigCoefI_slice, tmpCoefR_slice, tmpCoefI_slice, so3Coef );

		// inverse so(3)
		inverseSo3( so3Coef, so3Sig );

		double thisRWt = numeric::constants::d::pi * ( 4.0 * square((r_idx+1.0)) + 1.0/3.0 ) / sumRWt;

		for (int i=0; i<8*bw*bw*bw; ++i) {
			so3_correlation[i] += thisRWt * so3Sig[i].real();
		}
	}
}

// convert an index from 'so3Sig' into a rotation matrix
void SHT::idx_to_rot(int maxloc , numeric::xyzMatrix< double > & thisRot) {
	// uses current BW
	if (bw == 0) {
		thisRot = numeric::xyzMatrix< double >::identity();
		return;
	}

	int ii = maxloc / (4*bw*bw);
	int tmp = maxloc - (ii*4*bw*bw);
	int jj = tmp / (2*bw);
	tmp = maxloc - (ii*4*bw*bw) - jj*(2*bw);
	int kk = tmp ;

	double a = numeric::constants::d::pi*kk/((double) bw);
	double b = numeric::constants::d::pi*(2*ii+1)/(4.0*bw);
	double g = numeric::constants::d::pi*jj/((double) bw);

	thisRot(1,1) = (- sin(a)*cos(b)*sin(g) + cos(a)*cos(g));
	thisRot(1,2) = (cos(a)*cos(b)*sin(g) + sin(a)*cos(g));
	thisRot(1,3) = (sin(b)*sin(g));
	thisRot(2,1) = (-sin(a)*cos(b)*cos(g) - cos(a)*sin(g));
	thisRot(2,2) = (cos(a)*cos(b)*cos(g) - sin(a)*sin(g));
	thisRot(2,3) = (sin(b)*cos(g));
	thisRot(3,1) = (sin(a)*sin(b));
	thisRot(3,2) = (-cos(a)*sin(b));
	thisRot(3,3) = (cos(b));
}

// convert an index from 'so3Sig' into Euler angles alpha, beta, gamma
void SHT::idx_to_euler(int maxloc , numeric::xyzVector< double > & eulerAngles) {
	if (bw == 0) {
		eulerAngles = numeric::xyzVector< double >(0,0,0);
		return;
	}

	int ii = maxloc / (4*bw*bw);
	int tmp = maxloc - (ii*4*bw*bw);
	int jj = tmp / (2*bw);
	tmp = maxloc - (ii*4*bw*bw) - jj*(2*bw);
	int kk = tmp ;

	eulerAngles[0] = numeric::constants::d::pi*kk/((double) bw);
	eulerAngles[1] = numeric::constants::d::pi*(2*ii+1)/(4.0*bw);
	eulerAngles[2] = numeric::constants::d::pi*jj/((double) bw);
}


//   so3CombineCoef: combine the S^2 spherical coefficients
//     of two signals such that inverse SO(3) transform will
//     result in the correlation of the two
void
SHT::so3CombineCoef(
		utility::vector0< double > &sigCoefR,
		utility::vector0< double > &sigCoefI,
		utility::vector0< double > &tmpCoefR,
		utility::vector0< double > &tmpCoefI,
		ObjexxFCL::FArray3D< std::complex<double> > &so3Coef ) {

	int currsign, index ;
	double tmpSigCoefR, tmpSigCoefI ;
	double tmpPatCoefR, tmpPatCoefI ;

	// set all so3coefs to 0
	//int ncoeffs = (4*bw*bw*bw-bw)/3;
	so3Coef.dimension(2*bw,2*bw,2*bw);
	for (int i=0; i < 8*bw*bw*bw; ++i) so3Coef[i] = 0.0;

	for( int l = 0 ; l < bw ; l ++ ) {
		double wigNorm = 2.*numeric::constants::d::pi*sqrt(2./(2.*l+1.)) ;
		for( int m1 = -l ; m1 <= l ; m1 ++ ) {
			// grab signal coefficient, at this degree and order
			index = so3_.lm_index(-m1, l) ;
			tmpSigCoefR = sigCoefR[ index ];
			tmpSigCoefI = sigCoefI[ index ];

			currsign = ((m1 + l) % 2) ? -1 : 1 ;

			for( int m2 = -l ; m2 <= l ; m2 ++ ) {
				index = so3_.lm_index( -m2, l );

				tmpPatCoefR = currsign *   tmpCoefR[ index ];
				tmpPatCoefI = currsign * (-tmpCoefI[ index ]);

				// now multiply the signal coef by the pattern coef,
				// and save it in the so3 coefficient array

				// get the index
				index = ( l - std::max(std::abs(m1),std::abs(m2)) );
				if ( m1 >= 0 && m2 >= 0 ) {
					index += m1 * ( 6*bw*bw + 3*m1 - 2*m1*m1 - 1) / 6;
					for( int k = 0 ; k < m2 ; k ++ ) index += bw - std::max( m1,k ) ;
				} else if (m1 >= 0 && m2 < 0) { //  mp < 0
					index += -(m1*(1+m1)*(1+2*m1))/6 + bw*bw*(1+m1);
					for(  int k = m2; k < 0 ; k ++ ) index -= bw - std::max( m1,-k ) ;
				} else if (m1 < 0 && m2 >= 0) {  // m < 0
					index += (4*bw*bw*bw-bw)/3 - (m1)*(1-m1)*(1-2*m1) / 6 + bw*bw*m1 ;
					for( int  k = 0 ; k < m2 ; k ++ ) index += bw - std::max( -m1,k ) ;
				} else {
					index += (4*bw*bw*bw-bw)/3 + (m1+1)*(6*bw*bw-m1-2*m1*m1) / 6;
					for( int k = m2; k < 0 ; k ++ ) index -= bw - std::max( -m1,-k );
				}

				so3Coef[ index ] = std::complex<double> (
						wigNorm * ( tmpSigCoefR*tmpPatCoefR - tmpSigCoefI*tmpPatCoefI ),
						wigNorm * ( tmpSigCoefR*tmpPatCoefI + tmpSigCoefI*tmpPatCoefR ) ) ;

				// multiply currsign factor by -1
				currsign *= -1 ;
			}
		}
	}
}


// inverse SO3 transform
//   used to get ccs after calling above function
void
SHT::inverseSo3(
		ObjexxFCL::FArray3D< std::complex<double> > &so3Coef,
		ObjexxFCL::FArray3D< std::complex<double> > &so3Sig ) {

	int sampHere, coefHere, sampHere2 ;
	int n = 2 * bw ;

	// scratch space
	ObjexxFCL::FArray3D< std::complex<double> > scratch;
	scratch.dimension( 2*bw, 2*bw, 2*bw );

	// keep track of where in precompute_wig_trans_ we are pointing
	int wignerIdx = 0;

	// f_(0,0)
	sampHere = so3_.sampLoc( 0, 0 );
	coefHere = so3_.coefLoc( 0, 0 );
	so3_.wignerSynthesis( 0, 0, so3Coef, coefHere, precompute_wig_trans_, wignerIdx, so3Sig, sampHere ) ;

	wignerIdx += n * bw ;

	for ( int m1 = 1 ; m1 < bw ; m1 ++ ) {
		// f_(m1,m1)
		sampHere = so3_.sampLoc( m1, m1 ) ;
		coefHere = so3_.coefLoc( m1, m1 ) ;
		so3_.wignerSynthesis( m1, m1, so3Coef, coefHere, precompute_wig_trans_, wignerIdx, so3Sig, sampHere ) ;

		// f_(-m1,-m1)
		sampHere = so3_.sampLoc( m1, m1 );
		sampHere2 = so3_.sampLoc( -m1, -m1 );
		for ( int j = 0 ; j < 2*bw ; j ++ ) {
			so3Sig[sampHere2+j] = std::complex<double>( so3Sig[sampHere+j].real(), -so3Sig[sampHere+j].imag() );
		}

		// f_(-m1,m1)
		sampHere = so3_.sampLoc( -m1, m1 ) ;
		coefHere = so3_.coefLoc( -m1, m1 ) ;
		so3_.wignerSynthesisDiffSign( -m1, m1, so3Coef, coefHere, precompute_wig_trans_, wignerIdx, so3Sig, sampHere, scratch ) ;

		// f_(m1,-m1)
		sampHere = so3_.sampLoc( -m1, m1 );
		sampHere2 = so3_.sampLoc( m1, -m1 );
		for ( int j = 0 ; j < 2*bw ; j ++ ) {
			so3Sig[sampHere2+j] = std::complex<double>( so3Sig[sampHere+j].real(), -so3Sig[sampHere+j].imag() );
		}

		wignerIdx += n*(bw-m1) ;
	}

	for ( int m1 = 1 ; m1 < bw ; m1 ++ ) {
		// f_(m1,0)
		sampHere = so3_.sampLoc( m1, 0 ) ;
		coefHere = so3_.coefLoc( m1, 0 ) ;
		so3_.wignerSynthesis( m1, 0, so3Coef, coefHere, precompute_wig_trans_, wignerIdx, so3Sig, sampHere ) ;

		// f_(-m1,0)
		sampHere = so3_.sampLoc( m1, 0 );
		sampHere2 = so3_.sampLoc( -m1, 0 );
		for ( int j = 0 ; j < 2*bw ; j ++ ) {
			so3Sig[sampHere2+j] = std::complex<double>( so3Sig[sampHere+j].real(), -so3Sig[sampHere+j].imag() );
		}

		// f_(0,m1)
		sampHere = so3_.sampLoc( 0, m1 ) ;
		coefHere = so3_.coefLoc( 0, m1 ) ;
		so3_.wignerSynthesisSameSign( 0, m1, so3Coef, coefHere, precompute_wig_trans_, wignerIdx, so3Sig, sampHere ) ;

		// f_(0,-m1)
		sampHere = so3_.sampLoc( 0, m1 );
		sampHere2 = so3_.sampLoc( 0, -m1 );

		for ( int j = 0 ; j < 2*bw ; j ++ ) {
			so3Sig[sampHere2+j] = std::complex<double>( so3Sig[sampHere+j].real(), -so3Sig[sampHere+j].imag() );
		}

		//  advance Wigners ptr
		wignerIdx += n*(bw-m1) ;
	}

	//  1 <= m1 <= bw-1
	//  m1+1 <= m2 <= bw-1
	for ( int m1 = 1 ; m1 < bw ; m1 ++ ) {
		for ( int m2 = m1 + 1 ; m2 < bw ; m2 ++ ) {

			// f_(m1,m2)
			sampHere = so3_.sampLoc( m1, m2 ) ;
			coefHere = so3_.coefLoc( m1, m2 ) ;
			so3_.wignerSynthesis( m1, m2, so3Coef, coefHere, precompute_wig_trans_, wignerIdx, so3Sig, sampHere ) ;

			// f_(-m1,-m2)
			sampHere = so3_.sampLoc( m1, m2 );
			sampHere2 = so3_.sampLoc( -m1, -m2 );
			for ( int j = 0 ; j < 2*bw ; j ++ ) {
				so3Sig[sampHere2+j] = std::complex<double>( so3Sig[sampHere+j].real(), -so3Sig[sampHere+j].imag() );
			}

			// f_(m1,-m2)
			sampHere = so3_.sampLoc( m1, -m2 ) ;
			coefHere = so3_.coefLoc( m1, -m2 ) ;
			so3_.wignerSynthesisDiffSign( m1, -m2, so3Coef, coefHere, precompute_wig_trans_, wignerIdx, so3Sig, sampHere, scratch ) ;

			// f_(-m1,m2)
			sampHere = so3_.sampLoc( m1, -m2 );
			sampHere2 = so3_.sampLoc( -m1, m2 );

			for ( int j = 0 ; j < 2*bw ; j ++ ) {
				so3Sig[sampHere2+j] = std::complex<double>( so3Sig[sampHere+j].real(), -so3Sig[sampHere+j].imag() );
			}

			// f_(m2,m2)
			sampHere = so3_.sampLoc( m2, m1 ) ;
			coefHere = so3_.coefLoc( m2, m1 ) ;
			so3_.wignerSynthesisSameSign( m2, m1, so3Coef, coefHere, precompute_wig_trans_, wignerIdx, so3Sig, sampHere ) ;

			// f_(-m1,m2)
			sampHere = so3_.sampLoc( m2, m1 );
			sampHere2 = so3_.sampLoc( -m2, -m1 );

			for ( int j = 0 ; j < 2*bw ; j ++ ) {
				so3Sig[sampHere2+j] = std::complex<double>( so3Sig[sampHere+j].real(), -so3Sig[sampHere+j].imag() );
			}


			// f_(m2,-m1)
			sampHere = so3_.sampLoc( m2, -m1 ) ;
			coefHere = so3_.coefLoc( m2, -m1 ) ;
			so3_.wignerSynthesisDiffSign( m1, -m2, so3Coef, coefHere, precompute_wig_trans_, wignerIdx, so3Sig, sampHere, scratch ) ;

			// f_(-m2,m1)
			sampHere = so3_.sampLoc( m2, -m1 );
			sampHere2 = so3_.sampLoc( -m2, m1 );

			for ( int j = 0 ; j < 2*bw ; j ++ ) {
				so3Sig[sampHere2+j] = std::complex<double>( so3Sig[sampHere+j].real(), -so3Sig[sampHere+j].imag() );
			}

			//  advance Wigners ptr
			wignerIdx += n*(bw-m2) ;
		}
	}

	int dataIndex = (n)*(bw);

	for ( int m1 = 0 ; m1 < bw  ; m1 ++ ) {
		for (int i=0; i<n; ++i) so3Sig[dataIndex+i]=0.0;
		dataIndex += (2*n)*(bw) ;
	}

	dataIndex = bw*n*n;
	for (int i=0; i<n*n; ++i) so3Sig[dataIndex+i]=0.0;
	dataIndex += n*n + n*bw;

	for ( int m1 = 1 ; m1 < bw  ; m1 ++ ) {
		for (int i=0; i<n; ++i) so3Sig[dataIndex+i]=0.0;
		dataIndex += (2*n)*(bw) ;
	}

	// transpose
	transpose_so3( so3Sig, scratch, n, n*n ) ;

	// FFT each row: take n^2-many FFTs, each of length n.
	for (int i=0; i<n*n; ++i)
		kiss_fft( p1, &scratch[i*n], &so3Sig[i*n] );

	// transpose
	transpose_so3( so3Sig, scratch, n, n*n ) ;

	// FFT again
	for (int i=0; i<n*n; ++i) {
		kiss_fft( p1, &scratch[i*n], &so3Sig[i*n] );
	}

	//  normalize
	double dn = ( bw / numeric::constants::d::pi )/( n );
	for ( int j = 0 ; j < n*n*n; j++ ) so3Sig[j] *= dn;
}


// generate the cosine series for Pml at a specified value of m
void
SHT::setup_Pmls( ) {
	int k;

	precompute_pml_.resize( bw );
	precompute_pml_trans_.resize( bw );

	// scratch space
	utility::vector0< double > prev(bw), prevprev(bw);
	utility::vector0< double > arccosres(bw), x_i(bw), eval_args(bw), cosres(bw);

	for (int m=0; m<bw; m++) {
		int size_m = (bw/2-m%2)*(bw/2-m%2+1) + (bw/2)*(m%2) - (m/2)*((m/2)+1);
		precompute_pml_[m].resize( size_m, 0 );

		int tableptr = 0;

		kiss_dct_state p(bw,0);

		for (int i=0; i<bw; i++) {
			x_i[i] = cos((( 2.0*((double)i)+1.0 ) * numeric::constants::d::pi) / (2.0*bw));
			eval_args[i] = (( 2.0*((double)i)+1.0 ) * numeric::constants::d::pi) / (2.0*bw);
		}

		//  set initial values of first two Pmls
		for (int i=0; i<bw; i++) prevprev[i] = 0.0;

		if (m == 0) {
			for (int i=0; i<bw; i++)
				prev[i] = 1.0/sqrt(2.0);
		} else {
			//so3_.Pmm_L2(m, &eval_args[0], &prev[0]);
			double mcons = sqrt(m+0.5);
			for (int i=0; i<m; i++) mcons *= sqrt((m-(i/2.0))/(m-i));
			if ( m != 0 ) mcons *= pow(2.0,-m/2.0);

			if ((m%2) != 0) {
				for (int i=0; i<bw; i++) prev[i] = -mcons * pow(sin(eval_args[i]),(double) m-1);
			} else {
				for (int i=0; i<bw; i++) prev[i] = mcons * pow(sin(eval_args[i]),(double) m);
			}
		}

		if ((m % 2) == 0)
			k = m;
		else
			k = m-1;

		for (int i=0; i<bw; ++i) arccosres[i] = prev[i];

		kiss_dct( &p, &arccosres[0], &cosres[0]);
		cosres[0] *= 1.0/sqrt(2.0);
		double cos_scale = 1.0 / sqrt((double) bw);
		for ( int i = 0 ; i < bw ; i ++ ) cosres[i] *= cos_scale ;

		for (int i=0; i<=k; i+=2)
			precompute_pml_[m][i/2] = cosres[i];
		tableptr += k/2+1; // jump

		for (int i=0; i<bw-m-1; i++) {
			double A = sqrt( ((2.0*(m+i)+3.0)/(2.0*(m+i)+1)) * ((i+1.0)/(2.0*m+i+1.0)) ) * ((2*(m+i)+1.0)/(i+1.0));
			double C =  (m+i==0) ? 0.0 : -1.0 *
				sqrt( ((2.0*(m+i)+3.0)/(2.0*(m+i)-1)) * ((i+1.0)/(2.0*m+i+1.0)) * (i/(2.0*m+i)) ) *
				((2.0*m+i)/(i+1.0)) ;

			for (int j=0; j<bw; ++j) {
				arccosres[j] = C*prevprev[j]+A*prev[j]*x_i[j];
			}

			kiss_dct( &p, &arccosres[0], &cosres[0]);
			cosres[0] *= 1.0/sqrt(2.0);
			for ( int j = 0 ; j < bw ; j ++ ) cosres[j] *= cos_scale ;

			k++;

			if ( i % 2 )
				for (int j=0; j<=k; j+=2)
					precompute_pml_[m][tableptr+j/2] = cosres[j];
			else
				for (int j=1; j<=k; j+=2)
					precompute_pml_[m][tableptr+j/2] = cosres[j];

			//  update tableptr
			tableptr += k/2+1;

			for (int j=0; j<bw; ++j) {
				prevprev[j] = prev[j];
				prev[j] = arccosres[j];
			}
		}
	}


	// transpose
	for (int m=0; m<bw; m++) {
		precompute_pml_trans_[m].resize( precompute_pml_[m].size() );

		if ( m == bw - 1 ) {
			int ncoeffs = precompute_pml_[m].size();
			for (int i=0; i<ncoeffs; ++i) {
				precompute_pml_trans_[m][i]=precompute_pml_[m][i];
			}
		} else {
			
			int trans_tableptr=0, tableptr=0;
			
			for (int row = 0; row < bw; row++) {
				//  if m odd, no need to do last row - all zeroes
				if (row == (bw-1)) {
					if ( m % 2 )
						return;
				}

				//  compute the starting point for values in cos_pml_table
				if (row <= m) {
					if ((row % 2) == 0)
						tableptr = row/2;
					else
						tableptr = (m/2) + 1 + (row/2);
				} else {
					tableptr = 0;
					int rowlimit = (m%2==0)?row:row+1;
					for (int i=m; i<=rowlimit; i++) {
						int rowsize = (i<m) ? 0 : ( (m%2==0) ? ((i/2)+1) : (((i-1)/2)+1) );
						tableptr += rowsize;
					}
					tableptr--; // back 1
				}

				int stride = row + 2;
				if (row <= m) stride = m + 2 - (m % 2) + (row % 2);

				int rowsize = so3_.transposeRowSize(row, m);
				int costable_offset = 0;
				for (int i=0; i < rowsize; i++) {
					precompute_pml_trans_[m][trans_tableptr+i] = precompute_pml_[m][tableptr+costable_offset];
					costable_offset += stride;
					stride += 2;
				} //  closes i loop

				trans_tableptr += rowsize;
			} //  closes row loop
		}
	}
}


// precompute wigner d's up to the input bandwidth
void
SHT::setup_Wig( ) {
	precompute_wig_trans_.resize( ((bw*bw)*(2+3*bw+bw*bw))/3 );
	int wigner_idx = 0;

	int N = 2*bw ;

	// scratch space
	utility::vector0< double > sinPts(N);
	utility::vector0< double > cosPts(N);
	utility::vector0< double > sinPts2(N);
	utility::vector0< double > cosPts2(N);
	utility::vector0< double > scratch1(N);
	utility::vector0< double > scratch2(N);

	//  precompute sines and cosines
	for (int i=0; i<N; i++) {
		sinPts[i] = sin( (( 2.0*((double)i)+1.0 ) * numeric::constants::d::pi) / (2.0*N));
		cosPts[i] = cos( (( 2.0*((double)i)+1.0 ) * numeric::constants::d::pi) / (2.0*N));
		sinPts2[i] = sin( 0.5 * (( 2.0*((double)i)+1.0 ) * numeric::constants::d::pi) / (2.0*N));
		cosPts2[i] = cos( 0.5 * (( 2.0*((double)i)+1.0 ) * numeric::constants::d::pi) / (2.0*N));
	}

	//  precompute Wigner d's for m1 = m2 = 0
	so3_.genWigner_ds( 0, 0, cosPts, sinPts2, cosPts2, precompute_wig_trans_, wigner_idx, scratch1, scratch2 ) ;
	wigner_idx += bw * N;

	// precompute Wigner d's for abs(m1)=abs(m2)
	for ( int m1 = 1 ; m1 < bw ; m1 ++ ) {
		so3_.genWigner_ds( m1, m1, cosPts, sinPts2, cosPts2, precompute_wig_trans_, wigner_idx, scratch1, scratch2 ) ;
		wigner_idx += N * ( bw - m1 );
	}

	// precompute Wigner d's for one order being 0, and the other not
	for ( int m1 = 1 ; m1 < bw ; m1 ++ ) {
		so3_.genWigner_ds( m1, 0, cosPts, sinPts2, cosPts2, precompute_wig_trans_, wigner_idx, scratch1, scratch2 ) ;
		wigner_idx += N * ( bw - m1 );
	}

	// precompute Wigner d's for general m1, m2
	for ( int m1 = 1 ; m1 < bw ; m1 ++ ) {
		for ( int m2 = m1 + 1 ; m2 < bw ; m2 ++ ) {
			so3_.genWigner_ds( m1, m2, cosPts, sinPts2, cosPts2, precompute_wig_trans_, wigner_idx, scratch1, scratch2 ) ;
			wigner_idx += N * ( bw - m2 );
		}
	}
}

// make weights for Legendre transforms
void
SHT::setup_Weights(  ) {
	double scale = numeric::constants::d::pi/((double)(4*bw)) ;

	weights_.resize( 4*bw );

	for ( int j = 0 ; j < 2*bw ; j ++ ) {
		double tmpsum = 0.0 ;
		for ( int k = 0 ; k < bw ; k ++ )
			tmpsum += 1.0/((double)(2*k+1)) * sin((double)((2*j+1)*(2*k+1))*scale);
		tmpsum *= sin((double)(2*j+1)*scale) * 2./((double) bw) ;

		weights_[j] = tmpsum ;
		weights_[j+2*bw] = tmpsum * sin((double)(2*j+1)*scale);
	}
}


// forward SH transform in one 2D shell
void
SHT::forwardS2(
		utility::vector0< double > &rdata,
		utility::vector0< double > &rcoeffs,
		utility::vector0< double > &icoeffs ) {
	int N = 2*bw ;

	// scratch space
	utility::vector0< double > rres(N*N);
	utility::vector0< double > ires(N*N);
	utility::vector0< double > fltres(bw);
	utility::vector0< double > eval_pts(N);
	utility::vector0< double > scratch1(N);
	utility::vector0< double > scratch2(N);
	utility::vector0< double > idata(N,0.0); //dummy

	//  do the FFTs along phi
	for (int i=0; i<2*bw; ++i) {
		kiss_fft_split( fftPlan, &rdata[i*N], &idata[0], &rres[i], &ires[i], 1, 2*bw );
	}

	// normalize
	double scale = 1./ (2.0*bw) * sqrt( 2.0*numeric::constants::d::pi );
	for( int i = 0 ; i < N*N ; i ++ ) {
		rres[i] *= scale;
		ires[i] *= scale;
	}

	//  point to start of output data buffers
	int coeffs_ptr = 0;

	for (int m=0; m<bw; m++) {
		so3_.Legendre(rres, m*N, m, fltres, precompute_pml_[m], weights_, scratch1, scratch2, dctPlan);
		for (int i=0; i<bw-m; ++i) rcoeffs[coeffs_ptr+i] = fltres[i];
		so3_.Legendre(ires, m*N, m, fltres, precompute_pml_[m], weights_, scratch1, scratch2, dctPlan);
		for (int i=0; i<bw-m; ++i) icoeffs[coeffs_ptr+i] = fltres[i];

		coeffs_ptr += bw-m;
	}

	//  upper coefficients (since input is real, use symmetry)
	double pow_one = 1.0;
	for(int i = 1; i < bw; i++){
		pow_one *= -1.0;
		for(int j = i; j < bw; j++){
			int source = so3_.lm_index(i, j);
			int target = so3_.lm_index(-i, j);

			rcoeffs[target] = pow_one * rcoeffs[source];
			icoeffs[target] = -pow_one * icoeffs[source];
		}
	}
}


// inverse SH transform in one 2D shell
void
SHT::inverseS2(
		utility::vector0< double > &rcoeffs,
		utility::vector0< double > &icoeffs,
		utility::vector0< double > &rdata ) {
	int N = 2*bw ;

	// scratch
	utility::vector0< double > rfourdata(N*N);
	utility::vector0< double > ifourdata(N*N);
	utility::vector0< double > rinvfltres(N);
	utility::vector0< double > iminvfltres(N);
	utility::vector0< double > sin_values(N);
	utility::vector0< double > scratch(N);

	// filler (real signal only)
	utility::vector0< double > idata(N*N, 0.0);

	for (int i=0; i<N; i++) {
		sin_values[i] = sin( (( 2.0*((double)i)+1.0 ) * numeric::constants::d::pi) / (2.0*N) );
	}

	//  Now do all of the inverse Legendre transforms
	int coeffs_idx = 0;

	for (int m=0; m<bw; m++) {
		so3_.InvLegendre(rcoeffs, coeffs_idx, m, rinvfltres, precompute_pml_trans_[m], sin_values, scratch, idctPlan );
		so3_.InvLegendre(icoeffs, coeffs_idx, m, iminvfltres, precompute_pml_trans_[m], sin_values, scratch, idctPlan );

		//  will store normal, then tranpose before doing inverse fft
		for (int i=0; i<N; ++i) {
			rfourdata[m*N+i] = rinvfltres[i];
			ifourdata[m*N+i] = iminvfltres[i];
		}

		coeffs_idx += bw-m;
	}

	//  now fill in zero values where m = bw (from problem definition)
	for (int i=0; i<N; ++i) {
		rfourdata[bw*N+i] = 0.0;
		ifourdata[bw*N+i] = 0.0;
	}

	// the data is real, so invf(l,-m) = conj(invf(l,m)),
	for(int m = bw + 1; m < N; m++){
		for (int i=0; i<N; ++i) {
			rfourdata[m*N+i] = rfourdata[(N-m)*N+i];
			ifourdata[m*N+i] = ifourdata[(N-m)*N+i];
		}
		for(int i = 0; i < N; i++)
			ifourdata[(m*N)+i] *= -1.0;
	}

	//  normalize
	double sin_scale = 1./(sqrt(2.*numeric::constants::d::pi) );
	for(int i=0;i<4*bw*bw;i++) {
		rfourdata[i] *= sin_scale ;
		ifourdata[i] *= sin_scale ;
	}

	for (int i=0; i<2*bw; ++i) {
		kiss_fft_split( ifftPlan, &ifourdata[i], &rfourdata[i], &idata[0], &rdata[2*bw*i], 1, 2*bw ); // idata ignored
	}
}


}
}
