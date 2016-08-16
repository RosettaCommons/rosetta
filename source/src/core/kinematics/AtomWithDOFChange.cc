// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/AtomWithDOFChange.hh
/// @brief  Data structure for output-sensitie refold data class declaration
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/kinematics/AtomWithDOFChange.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

/// @brief Automatically generated serialization method
template< class Archive >
void
core::kinematics::AtomWithDOFChange::save( Archive & arc ) const {
	arc( CEREAL_NVP( atomid_ ) ); // id::AtomID
	arc( CEREAL_NVP( reached_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::kinematics::AtomWithDOFChange::load( Archive & arc ) {
	arc( atomid_ ); // id::AtomID
	arc( reached_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::kinematics::AtomWithDOFChange );
#endif // SERIALIZATION


