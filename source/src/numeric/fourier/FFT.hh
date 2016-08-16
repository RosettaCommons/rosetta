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

#ifndef INCLUDED_numeric_fourier_FFT_hh
#define INCLUDED_numeric_fourier_FFT_hh

// Package headers

// Project headers

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>

// C++ headers
#include <complex>


namespace numeric {
namespace fourier {

/// @brief 1D fft c->c double
void fft(ObjexxFCL::FArray1D< std::complex<double> > &X , ObjexxFCL::FArray1D< std::complex<double> > &fX);

/// @brief 1D inverse fft c->c double
void ifft(ObjexxFCL::FArray1D< std::complex<double> > &fX , ObjexxFCL::FArray1D< std::complex<double> > &X);

/// @brief 1D fft r->c float
void fft(ObjexxFCL::FArray1D< float >  &X , ObjexxFCL::FArray1D< std::complex<double> > &fX);

/// @brief 1D fft r->c double
void fft(ObjexxFCL::FArray1D< double >  &X , ObjexxFCL::FArray1D< std::complex<double> > &fX);

/// @brief 1D inverse ifft c->r float
void ifft(ObjexxFCL::FArray1D< std::complex<double> >  &fX , ObjexxFCL::FArray1D< float > &X);

/// @brief 1D inverse ifft c->r double
void ifft(ObjexxFCL::FArray1D< std::complex<double> >  &fX , ObjexxFCL::FArray1D< double > &X);


/// @brief 2D fft c->c double
void fft2(ObjexxFCL::FArray2D< std::complex<double> > &X , ObjexxFCL::FArray2D< std::complex<double> > &fX);

/// @brief 2D inverse fft c->c double
void ifft2(ObjexxFCL::FArray2D< std::complex<double> > &fX , ObjexxFCL::FArray2D< std::complex<double> > &X);

/// @brief 2D fft r->c float
void fft2(ObjexxFCL::FArray2D< float >  &X , ObjexxFCL::FArray2D< std::complex<double> > &fX);

/// @brief 2D fft r->c double
void fft2(ObjexxFCL::FArray2D< double >  &X , ObjexxFCL::FArray2D< std::complex<double> > &fX);

/// @brief 2D inverse ifft c->r float
void ifft2(ObjexxFCL::FArray2D< std::complex<double> >  &fX , ObjexxFCL::FArray2D< float > &X);

/// @brief 2D inverse ifft c->r double
void ifft2(ObjexxFCL::FArray2D< std::complex<double> >  &fX , ObjexxFCL::FArray2D< double > &X);


/// @brief 3D fft c->c double
void fft3(ObjexxFCL::FArray3D< std::complex<double> > const &X , ObjexxFCL::FArray3D< std::complex<double> > &fX);

/// @brief 3D inverse fft c->c double
void ifft3(ObjexxFCL::FArray3D< std::complex<double> > const &fX , ObjexxFCL::FArray3D< std::complex<double> > &X);

/// @brief 3D fft r->c float
void fft3(ObjexxFCL::FArray3D< float >  const &X , ObjexxFCL::FArray3D< std::complex<double> > &fX);

/// @brief 3D fft r->c double
void fft3(ObjexxFCL::FArray3D< double >  const &X , ObjexxFCL::FArray3D< std::complex<double> > &fX);

/// @brief 3D inverse ifft c->r float
void ifft3(ObjexxFCL::FArray3D< std::complex<double> > const &fX , ObjexxFCL::FArray3D< float > &X);

/// @brief 3D inverse ifft c->r double
void ifft3(ObjexxFCL::FArray3D< std::complex<double> > const &fX , ObjexxFCL::FArray3D< double > &X);

/// @brief 3D fft c->c double with no static
void fft3_dynamic(ObjexxFCL::FArray3D< std::complex<double> > &X , ObjexxFCL::FArray3D< std::complex<double> > &fX);

/// @brief 3D inverse fft c->c double with no static
void ifft3_dynamic(ObjexxFCL::FArray3D< std::complex<double> > &fX , ObjexxFCL::FArray3D< std::complex<double> > &X);
}
}

#endif
