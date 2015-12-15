// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/TorsionID_Range.cc
/// @brief  Conformation torsion identifier class
/// @author Colin Smith


// Unit headers
#include <core/id/TorsionID_Range.hh>

// C++ headers


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#endif // SERIALIZATION


namespace core {
namespace id {


/// @brief stream << TorsionID_Range
std::ostream &
operator <<( std::ostream & os, TorsionID_Range const & a )
{
	os << a.torsion_id() << " min= " << a.min() << " max= " << a.max() << ' ';
	return os;
}


} // namespace id
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::id::TorsionID_Range::save( Archive & arc ) const {
	arc( CEREAL_NVP( torsion_id_ ) ); // class core::id::TorsionID
	arc( CEREAL_NVP( min_ ) ); // core::Real
	arc( CEREAL_NVP( max_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::id::TorsionID_Range::load( Archive & arc ) {
	arc( torsion_id_ ); // class core::id::TorsionID
	arc( min_ ); // core::Real
	arc( max_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::id::TorsionID_Range );
#endif // SERIALIZATION

