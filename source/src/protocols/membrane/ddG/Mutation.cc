// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		protocols/membrane/ddG/Mutation.hh
///
/// @brief      Informaiton about a Point Mutation
/// @details    Object for scoring information about a point mutation - includes
///				the position to mutate to, and the specified amino acid.
///				Last Modified: 7/18/14
///
/// @remarks	This is pretty general and can probably live somewhere in core pack?
///				Not sure if anyone needs anything like this
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <protocols/membrane/ddG/Mutation.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/AA.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

using namespace core;
using namespace core::chemical;

namespace protocols {
namespace membrane {
namespace ddG {
	
/// @brief Default Constructor
/// @details Private default constructor - mutation to alanine at posiiton 1. Don't use
Mutation::Mutation() :
	utility::pointer::ReferenceCount(),
	position_( 1 ),
	aa_( core::chemical::aa_ala )
{}

/// @brief Custom Constructor
/// @details Construct a new mutation object from a position and amino acid
Mutation::Mutation( Size position, AA aa ) :
	utility::pointer::ReferenceCount(),
	position_( position ),
	aa_( aa )
{}

/// @brief Copy Constructor
/// @details Create a deep copy of this object
Mutation::Mutation( Mutation const & src ) :
	utility::pointer::ReferenceCount( src ),
	position_( src.position_ ),
	aa_( src.aa_ )
{}

/// @brief Destructor
Mutation::~Mutation() {}

/// @brief Get the posiiton at which to make this mutation
Size
Mutation::position() const {
	return position_;
}

/// @brief Get the amino acid to change this position to
AA
Mutation::aa() const {
	return aa_;
}

} // ddG
} // membrane
} // protocols

