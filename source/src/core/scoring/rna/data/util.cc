// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/rna/data/util.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/rna/data/util.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "core.scoring.rna.data.util" );

/////////////////////////////////////////////////////////////////////////////
// This is super-silly -- should replace with a class with
//  useful grid lookup and fill-in -- check other grid &
//  density map classes in Rosetta. DO THIS WHEN EXPANDING
//  chemical mapping scores beyond DMS!
//
// -- rhiju, 2014
/////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace rna {
namespace data {

	core::Size
	lookup_idx( Real const value, utility::vector1< Real > & values ) {
		if ( values.has_value( value ) ) return values.index( value );

		// values must be read in in increasing order -- sanity check.
		for ( Size n = 1; n <= values.size(); n++ ) runtime_assert( value > values[n] );

		values.push_back( value );
		return values.size();
	}

	//////////////////////////////////////////////////////////////////////////////////
	Size
	get_bool_idx( bool const value, utility::vector1< bool > const & values ){
		for ( Size n = 1; n <= values.size(); n++ ){
			if ( value == values[ n ] ) return n;
		}
		return 0;
	}

	//////////////////////////////////////////////////////////////////////////////////
	// Later upgrade to return how far into bin the value is, then carry out bilinear
	//  interpolation & derivs.
	Size
	get_idx( Real const value, utility::vector1< Real > const & values ){
		for ( Size n = 1; n < values.size(); n++ ){
			if ( value >= values[ n ] && value < values[ n+1 ] ) return n;
		}
		return values.size();
	}


} //data
} //rna
} //scoring
} //core
