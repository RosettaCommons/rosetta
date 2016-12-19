// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A class for defining atom parameters, known as atom_types
///
/// @details
/// This class contains the "chemical" information for atoms. This does not contain the actual
/// xyz coordinates of the class (xyz found in core/conformation/Orbital.hh. The atom_type properties
/// are assigned by the class OrbitalSet which is initiated from the ChemicalManager. Orbital type properties
/// are currently are read in from the file located chemical/atom_type_sets/fa_standard/atom_properties.txt.
/// These properties contain the the properties of LJ_RADIUS, LJ_WDEPTH, LK_DGRFREE, LK_LAMBDA, LK_VOLUME.
/// These properties are used in the scoring function fa_atr, fa_rep, fa_sol, which is located in the Etable (core/scoring/etable/Etable.hh)
/// Additional parameters are acceptor/donor, hybridization, and orbital parameters.
///
///
///
/// @author
/// Phil Bradley
/// Steven Combs - comments
///
///
/////////////////////////////////////////////////////////////////////////

// Rosetta headers
#include <core/chemical/Orbital.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

void
Orbital::remap_atom_vds( std::map< VD, VD > const & old_to_new ) {
	icoor_.remap_atom_vds( old_to_new );
	new_icoor_.remap_atom_vds( old_to_new );
}

void
Orbital::print(
	std::ostream & out
) const {
	out << "Name: " << name() << std::endl;
	out << "orbital_type_index: " << orbital_type_index() << std::endl;
	//out << "xyz: " << xyz() << std::endl;
	//out << "icoor: " << icoor() << std::endl;
	//out << "new_icoor: " << new_icoor() << std::endl;
	out << std::endl;
}

std::ostream &
operator<< (std::ostream & out, const Orbital & orbital ){
	orbital.print( out );
	return out;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


} // pose
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::Orbital::save( Archive & arc ) const {
	arc( CEREAL_NVP( name_ ) ); // std::string
	arc( CEREAL_NVP( orbital_type_index_ ) ); // Size
	arc( CEREAL_NVP( ideal_xyz_ ) ); // Vector
	arc( CEREAL_NVP( icoor_ ) ); // orbitals::ICoorOrbitalData
	arc( CEREAL_NVP( new_icoor_ ) ); // orbitals::ICoorOrbitalData
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::Orbital::load( Archive & arc ) {
	arc( name_ ); // std::string
	arc( orbital_type_index_ ); // Size
	arc( ideal_xyz_ ); // Vector
	arc( icoor_ ); // orbitals::ICoorOrbitalData
	arc( new_icoor_ ); // orbitals::ICoorOrbitalData
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::Orbital );
#endif // SERIALIZATION
