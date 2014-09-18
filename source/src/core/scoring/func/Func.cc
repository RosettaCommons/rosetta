// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/Func.cc
/// @brief Definition for functions used in definition of constraints.
/// @author Andrew Leaver-Fay
/// @author James Thompson
/// @author Oliver Lange

#include <core/scoring/func/Func.hh>
#include <core/scoring/func/Func.fwd.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

// AUTO-REMOVED #include <numeric/angle.functions.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

// C++ Headers
//#include <cstdlib>
//#include <iostream>
//#include <map>
//#include <utility>

namespace core {
namespace scoring {
namespace func {

/// @details Auto-generated virtual destructor
Func::~Func() {}

void
Func::read_data( std::istream& ) {
	basic::Tracer tr( "core.scoring.constraints.Func" );
	tr.Warning << " Base clase Func::read_data stubbed out ---  virtual function not overloaded " << std::endl;
}

Real
Func::estimate_dfunc( Real const r ) const {
	Real h = 1e-6;
	return estimate_dfunc( r, h );
}

Real Func::estimate_dfunc( Real const r, Real const h ) const {
	return ( (func(r+h) - func(r-h)) / (2*h) );
}

void Func::show( std::ostream& out ) const {
	using namespace ObjexxFCL::format;

	Real start = 2;
	Real end   = 20;
	Real res   = 0.5;
	int width  = 10;
	out << A( width, "r" )
			<< A( width, "func" )
			<< A( width, "dfunc")
			<< A( width, "dfunc_est" )
			<< std::endl;
	for ( Real r = start; r <= end; r += res ) {
		out << I( width, r )
				<< F( width, 3, func(r)  )
				<< F( width, 3, dfunc(r) )
				<< F( width, 3, estimate_dfunc(r) )
				<< std::endl;
	}
} // virtual void show( std::ostream& out )

void Func::show_definition( std::ostream &out ) const {
	out << "Func::show_def() stubbed out" << std::endl;
}

Size Func::show_violations( std::ostream& out, Real r, Size verbose_level,  Real threshold  ) const {
	Real f = func(r);
	if ( verbose_level > 100 ) out << f << std::endl;
	return f > threshold;
}

std::ostream& operator<<( std::ostream& out, const Func& f ) {
	f.show( out );
	return out;
}

} // namespace constraints
} // namespace scoring
} // namespace core
