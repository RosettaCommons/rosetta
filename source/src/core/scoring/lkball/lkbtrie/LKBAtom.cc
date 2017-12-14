// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/hbtrie/LKBAtom.hh
/// @brief
/// @author

// Unit Headers
#include <core/scoring/lkball/lkbtrie/LKBAtom.hh>

// Project Headers
#include <core/types.hh>

// STL Headers
#include <iostream>

// Numceric Headers
#include <numeric/xyzVector.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/fixedsizearray1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace lkball {
namespace lkbtrie {


LKBAtom::LKBAtom() : n_attached_waters_( 0 ) {}

LKBAtom::~LKBAtom() = default;

/// @brief send a description of the atom to standard out
void
LKBAtom::print() const { print( std::cout ); }

/// @brief send a description of the atom to an output stream
void
LKBAtom::print( std::ostream & os ) const
{
	os << "LKBAtom" <<  " ";
	os << "(" << base_.xyz().x();
	os << ", " << base_.xyz().y();
	os << ", " << base_.xyz().z() << ")";
}

std::ostream & operator << ( std::ostream & os, LKBAtom const & atom )
{
	atom.print( os );
	return os;
}

} // namespace lkbtrie
} // namespace lkball
} // namespace scoring
} // namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::lkball::lkbtrie::LKBAtom::save( Archive & arc ) const {
	arc( CEREAL_NVP( base_ ) ); // conformation::Atom
	arc( CEREAL_NVP( n_attached_waters_ ) ); // Size
	arc( CEREAL_NVP( waters_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( atom_weights_ ) ); // utility::vector1<Real>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::lkball::lkbtrie::LKBAtom::load( Archive & arc ) {
	arc( base_ ); // conformation::Atom
	arc( n_attached_waters_ ); // Size
	arc( waters_ ); // utility::vector1<Vector>
	arc( atom_weights_ ); // utility::vector1<Real>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::lkball::lkbtrie::LKBAtom );
#endif // SERIALIZATION

