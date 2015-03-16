// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/func/CircularSplineFunc.cc
/// @brief  Similar to spline func but periodic on [0,360)
/// @author fpd

// Unit Headers
#include <core/scoring/func/CircularSplineFunc.hh>

// Package Headers
#include <core/scoring/func/Func.fwd.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/database/open.hh>

// Project Headers
#include <basic/Tracer.hh>

// Utility and Numeric Headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <numeric/interpolation/spline/Interpolator.hh>
#include <numeric/interpolation/util.hh>
#include <numeric/MathVector.hh>

// C++ Headers
#include <iostream>
#include <sstream>
#include <string>

static thread_local basic::Tracer TR( "core.scoring.constraints.CircularSplineFunc" );

namespace core {
namespace scoring {
namespace func {

CircularSplineFunc::CircularSplineFunc( core::Real weight_in, utility::vector1< core::Real> energies_in ) {
	weight_ = weight_in;
	train( energies_in );
}

void
CircularSplineFunc::train( utility::vector1< core::Real> energies_in ) {
	using namespace numeric;
	using namespace numeric::interpolation::spline;

	runtime_assert( energies_in.size() == 36 );

	numeric::MathVector< core::Real > energies( 36 );
	for ( Size jj = 0; jj < 36; ++jj ) energies( jj ) = energies_in[jj+1];

	BorderFlag periodic_boundary = e_Periodic;
	Real start_val = 5.0; // grid is shifted by five degrees.
	Real delta = 10.0; // grid is 10 degrees wide
	std::pair< Real, Real > unused = std::make_pair( 0.0, 0.0 ); // unused

	spline_.train( periodic_boundary, start_val, delta, energies, unused );
}

// Read in data (e.g., experimental distance), weight, and histogram filename.  Bind filename to stream.
void CircularSplineFunc::read_data( std::istream &in) {
	using namespace numeric;
	using namespace numeric::interpolation::spline;

	numeric::MathVector< core::Real > energies( 36 );
	in >> weight_;
	for ( Size jj = 0; jj < 36; ++jj ) {
		in >> energies( jj );
	}

	BorderFlag periodic_boundary = e_Periodic;
	Real start_val = 5.0; // grid is shifted by five degrees.
	Real delta = 10.0; // grid is 10 degrees wide
	std::pair< Real, Real > unused = std::make_pair( 0.0, 0.0 ); // unused

	spline_.train( periodic_boundary, start_val, delta, energies, unused );
}

core::Real CircularSplineFunc::func( core::Real const x) const {
	core::Real fX = spline_.F(x);
	return weight_*fX;
}

core::Real CircularSplineFunc::dfunc( core::Real const x) const {
	core::Real dfX = spline_.dF(x);
	return weight_*dfX;
}

/// @brief show the definition of this CircularSplineFunc to the specified output stream.
void CircularSplineFunc::show_definition( std::ostream &out ) const {
	out << "CircularSplineFunc" << std::endl;
}

/// @brief show some sort of stringified representation of the violations for this constraint.
core::Size CircularSplineFunc::show_violations( std::ostream &out, core::Real x, core::Size verbose_level, core::Real threshold) const {
	return Func::show_violations( out, x, verbose_level, threshold);
} // show_violations()

} // constraints
} // scoring
} // core
