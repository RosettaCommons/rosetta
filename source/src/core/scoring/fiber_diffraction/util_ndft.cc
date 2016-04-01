// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffraction/util_ndft.cc
/// @brief  FiberDiffractionDens scoring
/// @author Ingemar Andre


#include <core/scoring/fiber_diffraction/util_ndft.hh>
#include <basic/Tracer.hh>
#include <utility/basic_sys_util.hh>

#ifdef WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif

namespace core {
namespace scoring {
namespace fiber_diffraction {

static basic::Tracer TR("core.scoring.fiber_diffraction.util_ndft");

void
FBnl_R(
	ObjexxFCL::FArray3D< float > & fourier_in,
	utility::vector0< utility::vector0 < int > >::iterator & nvals,
	core::Size & lmax,
	core::Size & rindex,
	//core::Real & c,
	ObjexxFCL::FArray3D< std::complex<float> > & fourier_out
) {
	int NZ=fourier_in.size3();      // z
	int ML=lmax+1;                  // l
	int NPHI=fourier_in.size2();    // phi
	ObjexxFCL::FArray2D< std::complex<float> > f1( ML, NPHI ); // ( l, z )
	for ( int nphi=1; nphi<=NPHI; ++nphi ) {                        // loop over each phi value

		ObjexxFCL::FArray1D< std::complex <float> > f0(NZ);
		ObjexxFCL::FArray1D< float > x0(ML);
		ObjexxFCL::FArray1D< std::complex <float> > f2(ML);

		for ( int nz=1; nz<=NZ; ++nz ) {                // loop over each z value
			std::complex< float> fij(fourier_in(rindex,nphi,nz)/(NZ), 0.0);
			f0(nz) = fij;
		}

		// layer lines start from 0...
		for ( int i=1; i<=ML; ++i ) {
			x0(i)=-float(i-1)/float(NZ);
		}

		ndft_1d(f0,x0,f2);

		// copy the result to f1
		for ( int ml=1; ml<=ML; ++ml ) {
			f1(ml,nphi)=f2(ml);
		}
	}

	for ( int ml=1; ml<=ML; ++ml ) {
		// number of bessel orders for layer line
		int MN (nvals[ml-1].size() );

		ObjexxFCL::FArray1D< int > nu;
		nu.dimension( nvals[ml-1].size() );
		for ( Size i=1; i<= nvals[ml-1].size(); ++i ) {
			nu(i) = nvals[ml-1][i-1];
		}

		ObjexxFCL::FArray1D< float > phi;
		phi.dimension(NPHI);
		for ( int i=0; i< NPHI; ++i ) {
			phi(i+1) = 2*M_PI*i/NPHI;
		}
		ObjexxFCL::FArray1D< std::complex <float> > fin;
		fin.dimension(NPHI);
		for ( int nphi=1; nphi<=NPHI; ++nphi ) {
			fin(nphi) = f1(ml,nphi);
		}
		ObjexxFCL::FArray1D< std::complex <float> > fout;
		fout.dimension(MN);

		ndft(fin,phi,nu,fout);

		for ( int mn=1; mn<=MN; ++mn ) {
			fourier_out(rindex,ml,mn)=fout(mn);
		}
	}
}


void
ft_nfft(
	ObjexxFCL::FArray3D< float > & fourier_in,
	utility::vector0< utility::vector0 < int > >::iterator & nvals,
	core::Size & lmax,
	core::Size & total_rvals,
	ObjexxFCL::FArray3D< std::complex<float> > & fourier_out
) {
	//Main FFT loop
	for ( Size ri=1; ri<= total_rvals; ++ri ) {
		FBnl_R( fourier_in, nvals, lmax, ri, fourier_out );
	}
}

void
ndft(
	ObjexxFCL::FArray1D< std::complex <float> > & fin,
	ObjexxFCL::FArray1D< float > & phi,
	ObjexxFCL::FArray1D< int > & n,
	ObjexxFCL::FArray1D< std::complex <float> > & fout
) {
	Size Nt=fin.size();
	Size Nphi=phi.size();
	Size Nn=n.size();

	runtime_assert(Nt==Nphi);

	for ( Size i=1; i<=Nn; ++i ) {
		fout(i)=0.0;
		for ( Size j=1; j<=Nphi; ++j ) {
			double re = cos(-phi(j)*n(i))/Nphi;
			double im = sin(-phi(j)*n(i))/Nphi;
			std::complex< float> wij( re, im );
			fout(i)+= wij*fin(j);
		}
	}
}

void
ndft_1d(
	ObjexxFCL::FArray1D < std::complex< float> >  & f0,
	ObjexxFCL::FArray1D< float > & x0,
	ObjexxFCL::FArray1D< std::complex< float> > & f2
) {
	Size NZ=f0.size();
	Size ML=f2.size();
	int m;
	for ( Size l=1; l<=ML; ++l ) {
		f2(l)=0.0;
		for ( Size k=1; k<=NZ; ++k ) {
			m = k - (NZ/2+1);
			float arg1(1.0);
			float arg2(-2*M_PI*m*x0(l));
			f2(l) += f0(k)* std::polar<float>( arg1, arg2 );
		}
	}
}


} // fiber_diffraction
} // scoring
} // core
