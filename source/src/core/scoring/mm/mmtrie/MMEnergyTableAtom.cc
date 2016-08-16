// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/mm/mmtrie/MMEnergyTableAtom.cc
/// @brief  Implimentation for the MMEnergyTableAtom. Heavily coppied from the EtableAtom.cc
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit Headers
#include <core/scoring/mm/mmtrie/MMEnergyTableAtom.hh>

// Project Headers
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>

// STL Headers
#include <iostream>

// Numceric Headers
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace mm {
namespace mmtrie {


MMEnergyTableAtom::MMEnergyTableAtom() : parent(), is_hydrogen_( false ) {}

MMEnergyTableAtom::MMEnergyTableAtom( conformation::Residue const & res, Size atom_index )
:
	parent( res.atom( atom_index ) ),
	is_hydrogen_( false )
{}

MMEnergyTableAtom::~MMEnergyTableAtom() {}

/// @brief send a description of the atom to standard out
void
MMEnergyTableAtom::print() const { print( std::cout ); }

/// @brief send a description of the atom to an output stream
void
MMEnergyTableAtom::print( std::ostream & os ) const
{
	os << "mm atom type" << mm_type() << " ";
	os << "(" << xyz().x();
	os << ", " << xyz().y();
	os << ", " << xyz().z() << ")" << std::endl;
}

std::ostream & operator << ( std::ostream & os, MMEnergyTableAtom const & atom )
{
	atom.print( os );
	return os;
}

} // namespace mmtrie
} // namespace mm
} // namespace scoring
} // namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::mm::mmtrie::MMEnergyTableAtom::save( Archive & arc ) const {
	arc( cereal::base_class< core::conformation::Atom >( this ) );
	arc( CEREAL_NVP( is_hydrogen_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::mm::mmtrie::MMEnergyTableAtom::load( Archive & arc ) {
	arc( cereal::base_class< core::conformation::Atom >( this ) );
	arc( is_hydrogen_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::mm::mmtrie::MMEnergyTableAtom );
CEREAL_REGISTER_TYPE( core::scoring::mm::mmtrie::MMEnergyTableAtom )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_mm_mmtrie_MMEnergyTableAtom )
#endif // SERIALIZATION
