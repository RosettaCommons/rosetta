// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  Wrapper class performing ffts on 1,2, and 3 dimensional farrays
/// @author Frank DiMaio

// Package headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>


// Project headers
#include <numeric/fourier/kiss_fft.hh>
#include <numeric/fourier/FFT.hh>


namespace numeric {
namespace fourier {

/// @brief 1D fft c->c double
void fft(ObjexxFCL::FArray1D< std::complex<double> > &X , ObjexxFCL::FArray1D< std::complex<double> > &fX) {
	kiss_fft_state fft_params;
	fft_params.resize( X.I1().size(), 0 );
	fX.dimension(X.I1().size());  // resizes "on demand"
	kiss_fft(&fft_params, &X[0], &fX[0] );
}

/// @brief 1D inverse fft c->c double
void ifft(ObjexxFCL::FArray1D< std::complex<double> > &fX , ObjexxFCL::FArray1D< std::complex<double> > &X) {
	kiss_fft_state ifft_params;
	ifft_params.resize( fX.I1().size(), 1 );
	X.dimension(fX.I1().size());
	kiss_fft(&ifft_params, &fX[0], &X[0] );

	// rescale
	int dimsProd = X.I1().size();
	for ( int i=0; i<dimsProd; ++i ) X[i] /= (double)dimsProd;
}

/// @brief 1D fft r->c double
void fft(ObjexxFCL::FArray1D< double > &X , ObjexxFCL::FArray1D< std::complex<double> > &fX) {
	ObjexxFCL::FArray1D< std::complex<double> > Xcpx;
	Xcpx.dimension (X.I1().size());
	int dimsProd = X.I1().size();
	for ( int i=0; i<dimsProd; ++i ) Xcpx[i] = X[i];
	fft(Xcpx,fX);
}

/// @brief 1D inverse ifft c->r double
void ifft(ObjexxFCL::FArray1D< std::complex<double> > &fX , ObjexxFCL::FArray1D< double > &X) {
	ObjexxFCL::FArray1D< std::complex<double> > Xcpx;
	ifft( fX,Xcpx );
	X.dimension(fX.I1().size());
	int dimsProd = X.I1().size();
	for ( int i=0; i<dimsProd; ++i ) X[i] = (float)Xcpx[i].real();
}

/// @brief 1D fft r->c float .. wraps double version
void fft(ObjexxFCL::FArray1D< float > &X , ObjexxFCL::FArray1D< std::complex<double> > &fX) {
	ObjexxFCL::FArray1D< std::complex<double> > Xcpx;
	Xcpx.dimension (X.I1().size());
	int dimsProd = X.I1().size();
	for ( int i=0; i<dimsProd; ++i ) Xcpx[i] = X[i];
	fft(Xcpx,fX);
}

/// @brief 1D inverse ifft c->r float ... wraps double version
void ifft(ObjexxFCL::FArray1D< std::complex<double> > &fX , ObjexxFCL::FArray1D< float > &X) {
	ObjexxFCL::FArray1D< std::complex<double> > Xcpx;
	ifft( fX,Xcpx );
	X.dimension(fX.I1().size());
	int dimsProd = X.I1().size();
	for ( int i=0; i<dimsProd; ++i ) X[i] = (float)Xcpx[i].real();
}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

/// @brief 2D fft c->c double
void fft2(ObjexxFCL::FArray2D< std::complex<double> > &X , ObjexxFCL::FArray2D< std::complex<double> > &fX) {
	kiss_fftnd_state fft_params;
	std::vector< int > dims(2);
	dims[0]=X.I2().size(); dims[1]=X.I1().size();  /// FArray dimensions in reverse order
	fft_params.resize( dims, 0 );
	fX.dimension(X.I1().size(),X.I2().size());
	kiss_fftnd(&fft_params, &X[0], &fX[0] );
}

/// @brief 2D inverse fft c->c double
void ifft2(ObjexxFCL::FArray2D< std::complex<double> > &fX , ObjexxFCL::FArray2D< std::complex<double> > &X) {
	kiss_fftnd_state ifft_params;
	std::vector< int > dims(2);
	dims[0]=fX.I2().size(); dims[1]=fX.I1().size();  /// FArray dimensions in reverse order
	ifft_params.resize( dims, 1 );
	X.dimension(fX.I1().size(),fX.I2().size());
	kiss_fftnd(&ifft_params, &fX[0], &X[0] );

	// rescale
	int dimsProd = X.I1().size()*X.I2().size();
	for ( int i=0; i<dimsProd; ++i ) X[i] /= (double)dimsProd;
}

/// @brief 2D fft r->c double
void fft2(ObjexxFCL::FArray2D< double > &X , ObjexxFCL::FArray2D< std::complex<double> > &fX) {
	ObjexxFCL::FArray2D< std::complex<double> > Xcpx;
	Xcpx.dimension (X.I1().size(),X.I2().size());
	int dimsProd = X.I1().size()*X.I2().size();
	for ( int i=0; i<dimsProd; ++i ) Xcpx[i] = X[i];
	fft2(Xcpx,fX);
}

/// @brief 2D inverse ifft c->r double
void ifft2(ObjexxFCL::FArray2D< std::complex<double> > &fX , ObjexxFCL::FArray2D< double > &X) {
	ObjexxFCL::FArray2D< std::complex<double> > Xcpx;
	ifft2( fX,Xcpx );
	X.dimension(fX.I1().size(),fX.I2().size());
	int dimsProd = fX.I1().size()*fX.I2().size();
	for ( int i=0; i<dimsProd; ++i ) X[i] = (double)Xcpx[i].real();
}

/// @brief 2D fft r->c float
void fft2(ObjexxFCL::FArray2D< float > &X , ObjexxFCL::FArray2D< std::complex<double> > &fX) {
	ObjexxFCL::FArray2D< std::complex<double> > Xcpx;
	Xcpx.dimension (X.I1().size(),X.I2().size());
	int dimsProd = X.I1().size()*X.I2().size();
	for ( int i=0; i<dimsProd; ++i ) Xcpx[i] = X[i];
	fft2(Xcpx,fX);
}

/// @brief 2D inverse ifft c->r float
void ifft2(ObjexxFCL::FArray2D< std::complex<double> > &fX , ObjexxFCL::FArray2D< float > &X) {
	ObjexxFCL::FArray2D< std::complex<double> > Xcpx;
	ifft2( fX,Xcpx );
	X.dimension(fX.I1().size(),fX.I2().size());
	int dimsProd = fX.I1().size()*fX.I2().size();
	for ( int i=0; i<dimsProd; ++i ) X[i] = (float)Xcpx[i].real();
}


//////////////////////////////////////////////////
//////////////////////////////////////////////////

/// @brief 3D fft c->c double
void fft3(ObjexxFCL::FArray3D< std::complex<double> > const &X , ObjexxFCL::FArray3D< std::complex<double> > &fX) {
	kiss_fftnd_state fft_params;
	std::vector< int > dims(3);
	dims[0]=X.I3().size(); dims[1]=X.I2().size(); dims[2]=X.I1().size();  /// FArray dimensions in reverse order
	fft_params.resize( dims, 0 );
	fX.dimension(X.I1().size(),X.I2().size(),X.I3().size());
	kiss_fftnd(&fft_params, &X[0], &fX[0] );
}

/// @brief 3D inverse fft c->c double
void ifft3(ObjexxFCL::FArray3D< std::complex<double> > const &fX , ObjexxFCL::FArray3D< std::complex<double> > &X) {
	kiss_fftnd_state ifft_params;
	std::vector< int > dims(3);
	dims[0]=fX.I3().size(); dims[1]=fX.I2().size(); dims[2]=fX.I1().size();  /// FArray dimensions in reverse order
	ifft_params.resize( dims, 1 );
	X.dimension(fX.I1().size(),fX.I2().size(),fX.I3().size());
	kiss_fftnd(&ifft_params, &fX[0], &X[0] );

	// rescale
	int dimsProd = X.I1().size()*X.I2().size()*X.I3().size();
	for ( int i=0; i<dimsProd; ++i ) X[i] /= (double)dimsProd;
}

/////////////////////////////////////
//3D FFT and inverse-FFT with dynamic variable allocation.
//Avoiding static variables occpupy unnecessary space when the data is large.
/// @brief 3D fft c->c double with no static
void fft3_dynamic(ObjexxFCL::FArray3D< std::complex<double> > &X , ObjexxFCL::FArray3D< std::complex<double> > &fX) {
	kiss_fftnd_state fft_params;
	std::vector< int > dims(3);
	dims[0]=X.I3().size(); dims[1]=X.I2().size(); dims[2]=X.I1().size();  /// FArray dimensions in reverse order
	fft_params.resize( dims, 0 );
	fX.dimension(X.I1().size(),X.I2().size(),X.I3().size());
	kiss_fftnd(&fft_params, &X[0], &fX[0] );
}

/// @brief 3D inverse fft c->c doublewith no static
void ifft3_dynamic(ObjexxFCL::FArray3D< std::complex<double> > &fX , ObjexxFCL::FArray3D< std::complex<double> > &X) {
	kiss_fftnd_state ifft_params;
	std::vector< int > dims(3);
	dims[0]=fX.I3().size(); dims[1]=fX.I2().size(); dims[2]=fX.I1().size();  /// FArray dimensions in reverse order
	ifft_params.resize( dims, 1 );
	X.dimension(fX.I1().size(),fX.I2().size(),fX.I3().size());
	kiss_fftnd(&ifft_params, &fX[0], &X[0] );

	// rescale
	int dimsProd = X.I1().size()*X.I2().size()*X.I3().size();
	for ( int i=0; i<dimsProd; ++i ) X[i] /= (double)dimsProd;
}
//////////////////////////////////

/// @brief 3D fft r->c double
void fft3(ObjexxFCL::FArray3D< double > const &X , ObjexxFCL::FArray3D< std::complex<double> > &fX) {
	ObjexxFCL::FArray3D< std::complex<double> > Xcpx;
	Xcpx.dimension (X.I1().size(),X.I2().size(),X.I3().size());
	int dimsProd = X.I1().size()*X.I2().size()*X.I3().size();
	for ( int i=0; i<dimsProd; ++i ) Xcpx[i] = X[i];
	fft3(Xcpx,fX);
}

/// @brief 3D inverse ifft c->r double
void ifft3(ObjexxFCL::FArray3D< std::complex<double> > const &fX , ObjexxFCL::FArray3D< double > &X) {
	ObjexxFCL::FArray3D< std::complex<double> > Xcpx;
	ifft3( fX,Xcpx );
	X.dimension(fX.I1().size(),fX.I2().size(),fX.I3().size());
	int dimsProd = fX.I1().size()*fX.I2().size()*fX.I3().size();
	for ( int i=0; i<dimsProd; ++i ) X[i] = (double) Xcpx[i].real();
}

/// @brief 3D fft r->c float
void fft3(ObjexxFCL::FArray3D< float > const &X , ObjexxFCL::FArray3D< std::complex<double> > &fX) {
	ObjexxFCL::FArray3D< std::complex<double> > Xcpx;
	Xcpx.dimension (X.I1().size(),X.I2().size(),X.I3().size());
	int dimsProd = X.I1().size()*X.I2().size()*X.I3().size();
	for ( int i=0; i<dimsProd; ++i ) Xcpx[i] = X[i];
	fft3(Xcpx,fX);
}

/// @brief 3D inverse ifft c->r float
void ifft3(ObjexxFCL::FArray3D< std::complex<double> > const &fX , ObjexxFCL::FArray3D< float > &X) {
	ObjexxFCL::FArray3D< std::complex<double> > Xcpx;
	ifft3( fX,Xcpx );
	X.dimension(fX.I1().size(),fX.I2().size(),fX.I3().size());
	int dimsProd = fX.I1().size()*fX.I2().size()*fX.I3().size();
	for ( int i=0; i<dimsProd; ++i ) X[i] = (float) Xcpx[i].real();
}


}
}

