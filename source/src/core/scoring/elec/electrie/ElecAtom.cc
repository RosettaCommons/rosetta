// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/hbtrie/ElecAtom.hh
/// @brief
/// @author

// Unit Headers
#include <core/scoring/elec/electrie/ElecAtom.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>

// STL Headers
#include <iostream>

// Numceric Headers
#include <numeric/xyzVector.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace elec {
namespace electrie {



ElecAtom::ElecAtom() : frac_(0.0), isbb_( false ), is_hydrogen_( false ), is_wat_( false ), charge_( 0.0 ) {}

ElecAtom::ElecAtom( conformation::Residue const & res, Size atom_index )
:
	base_( res.atom( atom_index ) ),
	frac_(0.0),
	isbb_( res.atom_is_backbone( atom_index ) ),
	is_hydrogen_( false ),
	is_wat_( false ), // hydrate/SPaDES protocol
	charge_( res.atomic_charge( atom_index ) )
{}

ElecAtom::~ElecAtom() = default;

/// @brief send a description of the atom to standard out
void
ElecAtom::print() const { print( std::cout ); }

/// @brief send a description of the atom to an output stream
void
ElecAtom::print( std::ostream & os ) const
{
	os << "ElecAtom" <<  " ";
	os << "(" << base_.xyz().x();
	os << ", " << base_.xyz().y();
	os << ", " << base_.xyz().z() << ")";
}

std::ostream & operator << ( std::ostream & os, ElecAtom const & atom )
{
	atom.print( os );
	return os;
}

} // namespace electrie
} // namespace elec
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::elec::electrie::ElecAtom::save( Archive & arc ) const {
	arc( CEREAL_NVP( base_ )  );
	arc( CEREAL_NVP( offsite_ ) ); // vector1
	arc( CEREAL_NVP( frac_ ) ); // Real
	arc( CEREAL_NVP( isbb_ ) ); // _Bool
	arc( CEREAL_NVP( is_hydrogen_ ) ); // _Bool
	arc( CEREAL_NVP( is_wat_ ) ); // _Bool
	arc( CEREAL_NVP( charge_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::elec::electrie::ElecAtom::load( Archive & arc ) {
	arc( base_ );
	arc( offsite_ ); // vector1
	arc( frac_ ); // Real
	arc( isbb_ ); // _Bool
	arc( is_hydrogen_ ); // _Bool
	arc( is_wat_ ); // _Bool
	arc( charge_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::elec::electrie::ElecAtom );

#endif // SERIALIZATION
