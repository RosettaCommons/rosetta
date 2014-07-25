// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/conformation/MembranePlanes.cc
///
/// @brief 		Specification for Membrane Planes Residue Definitions
///	@details	When using the membrane code, users can optionally view the membrane planes
///				defined by the center/normal positions. This object will store the positions of anchoring residues
///				which will be used to draw CGO planes in PyMOl, representing the membrane planes.
///				This data should not be used in dynamic simulations, it's only purpose is for visualization.
///				Last Modified: 7/23/14
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <core/conformation/membrane/MembranePlanes.hh>

// Project Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <iostream>

namespace core {
namespace conformation {
namespace membrane {

using namespace core;

////////////////////
/// Constructors ///
////////////////////

/// @brief Custom Constructor
/// @details Create a MembranePlanes object detailing the positions of
/// the upper and lower membrane planes
MembranePlanes::MembranePlanes(
	utility::vector1< Size > top_points,
	utility::vector1< Size > bottom_points
	) :
	utility::pointer::ReferenceCount(),
	top_points_( top_points ),
	bottom_points_( bottom_points )
{
	// TODO: Add error checking!
}

/// @brief Copy Constructor
/// @details Make a deep copy of this membrane planes object
MembranePlanes::MembranePlanes( MembranePlanes const & src ) :
	utility::pointer::ReferenceCount(),
	top_points_( src.top_points_ ),
	bottom_points_( src.bottom_points_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this obejct overriding "="
MembranePlanes &
MembranePlanes::operator=( MembranePlanes const & src ) {
	
	// Abort self-assignment.
	if (this == &src) {
		return *this;
	}
	
	// Otherwise, create a new object
	return *( new MembranePlanes( *this ) );
	
}

/// @brief Destructor
MembranePlanes::~MembranePlanes() {}

/// @brief Show Membrane Planes Info Data
void
MembranePlanes::show( std::ostream & output ) const {
	
	output << "Membrane Plane Points" << std::endl;
	
	// Print top point residue positions
	output << "Upper Plane defined at: ";
	for ( core::Size i = 1; i <= top_points_.size(); ++i ) {
		output << utility::to_string( top_points_[i] ) << " ";
	}
	output << " " << std::endl;
	
	// Pring bottom point residue positions
	output << "Lower Plane defined at: ";
	for ( core::Size i = 1; i <= bottom_points_.size(); ++i ) {
		output << utility::to_string( bottom_points_[i] ) << " ";
	}
	output << " " << std::endl;
}

///////////////////////////
/// Data Access Methods ///
///////////////////////////

/// @brief Access point residues defining top membrane plane
utility::vector1< Size >
MembranePlanes::top_points() {
	return top_points_;
}

/// @brief Access point residues defining bottom membrane planes
utility::vector1< Size >
MembranePlanes::bottom_points() {
	return bottom_points_;
}

/// @brief Default Constructor
MembranePlanes::MembranePlanes() :
	utility::pointer::ReferenceCount()
{}


} // membrane
} // conformation
} // core
