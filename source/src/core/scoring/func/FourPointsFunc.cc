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


#include <core/scoring/func/FourPointsFunc.hh>

#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace func {

FourPointsFunc::FourPointsFunc() : points_( 4, Vector(0.0) ) {}

/// @brief set the coordinate for one of the four atoms
void FourPointsFunc::xyz( Size atomno, Vector const & coord )
{
	if ( atomno < 1 || atomno > 4 ) {
		utility_exit_with_message( "Error: FourPointsFunc can only store coordinates for four atoms; id " + utility::to_string( atomno ) + " is invalid." );
	}
	points_[ atomno ] = coord;
}

FourPointsFunc::~FourPointsFunc() {}

Vector const &
FourPointsFunc::operator()( AtomID const & id ) const
{
	if ( id.rsd() != 1 ) {
		utility_exit_with_message( "Error: invalid AtomID for FourPointsFunc.  Must request residu 1.  Requested residue " + utility::to_string( id.rsd() ) + " instead");
	}
	if ( id.atomno() < 1 || id.atomno() > 4 ) {
		utility_exit_with_message( "Error: FourPointsFunc can only store coordinates for four atoms; id " + utility::to_string( id.atomno() ) + " is invalid." );
	}
	return points_[ id.atomno() ];
}

conformation::Residue const &
FourPointsFunc::residue( Size ) const
{
	utility_exit_with_message( "FourPointsFunc does not implement a residue() method" );

	// unreachable -- appease compiler
	chemical::ResidueType rt( NULL, NULL, NULL, NULL );
	static conformation::Residue r( rt, true );
	return r;
}

} // constraints
} // scoring
} // core
