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

#ifndef INCLUDED_protocols_cryst_util_hh
#define INCLUDED_protocols_cryst_util_hh

#include <core/types.hh>

#include <fstream>
#include <iostream>
#include <math.h>

#include <sstream>
#include <string>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/fourier/FFT.hh>
#include <numeric/constants.hh>


#define DEG2RAD 0.0174532925199433
#define RAD2DEG 57.295779513082323

namespace protocols {
namespace cryst {


inline core::Real absfpart( core::Real x) {
	return ( std::fabs( x-std::floor( x+0.5) ) );
}

inline core::Size denom( core::Real x ) {
	if ( absfpart(x) <= 1e-4 ) { return 1; }
	for ( core::Size i=2; i<=6; ++i ) {
		if ( absfpart(i*x) <= 1e-4 ) { return i;}
	}

	utility_exit_with_message( "error in denom()");
	return 0;
}

inline int pos_mod(int x,int y) {
	int r=x%y; if ( r<0 ) r+=y;
	return r;
}
inline core::Real pos_mod(core::Real x,core::Real y) {
	core::Real r=std::fmod(x,y); if ( r<0 ) r+=y;
	return r;
}
inline int min_mod(int x,int y) {
	int r=x%y; if ( r<-y/2 ) r+=y;
	if ( r>=y/2 ) r-=y;
	return r;
}
inline core::Real min_mod(core::Real x,core::Real y) {
	core::Real r=std::fmod(x,y); if ( r<-0.5*y ) r+=y;
	if ( r>=0.5*y ) r-=y;
	return r;
}


}
}

#endif
