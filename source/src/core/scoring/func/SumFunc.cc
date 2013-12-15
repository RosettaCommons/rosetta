// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/SumFunc.cc
/// @brief Weighted constraint function that reweights other constraints
/// by a constant scalar value.
/// @author James Thompson, Greg Taylor


#include <core/scoring/func/SumFunc.hh>
#include <core/scoring/func/FuncFactory.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

// AUTO-REMOVED #include <numeric/angle.functions.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>

// C++ Headers
// AUTO-REMOVED #include <string>

#include <sstream>


namespace core {
namespace scoring {
namespace func {

	Real
	SumFunc::func( Real const x ) const {
		Real f( 0.0 );
		for( utility::vector1< FuncOP >::const_iterator it = funcs_.begin(), end = funcs_.end();
					it != end; ++it
		) {
			f += (*it)->func( x );
		}

		return f;
	}

	Real
	SumFunc::dfunc( Real const x ) const {
		Real df( 0.0 );
		for( utility::vector1< FuncOP >::const_iterator it = funcs_.begin(), end = funcs_.end();
					it != end; ++it
		) {
			df += (*it)->dfunc( x );
		}

		return df;
	}

	void
	SumFunc::read_data( std::istream& in ) {
		FuncFactory func_factory;
		std::string func_type;

		Size n_funcs( 0 );
		in >> n_funcs;
		if ( in.fail() ) return;

		assert( n_funcs >= 1 );

		for ( Size i = 1; i <= n_funcs; ++i ) {
			in >> func_type;
			FuncOP current_func;

			current_func = func_factory.func_types_[ func_type ]->clone();
    	current_func->read_data( in );
			add_func( current_func );
		}
	} // read_data

} // namespace constraints
} // namespace scoring
} // namespace core
