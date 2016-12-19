// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/MMAtomTypeSet.hh
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


#ifndef INCLUDED_core_chemical_MMAtomTypeSet_hh
#define INCLUDED_core_chemical_MMAtomTypeSet_hh


// Unit headers
#include <core/chemical/MMAtomTypeSet.fwd.hh>

// Project headers

// Utility headers
#include <utility/exit.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>

#include <core/chemical/MMAtomType.fwd.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
#include <utility/serialization/serialization.hh>
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {


/// @brief A set of MMAtomTypes
///
/// @details This class contains a vector of pointers each of which points to an
/// MMAtomType and the vector index is looked up by an atom_name string
/// in a map.
///
class MMAtomTypeSet : public utility::pointer::ReferenceCount {

public:
	MMAtomTypeSet( std::string const & name = "" );
	~MMAtomTypeSet() override;

	/// @brief What the ChemicalManager knows this as, if relevant
	std::string const &
	name() const { return name_; }

	/// @brief Number of MM atom types in the set
	Size
	n_atomtypes() const
	{
		return atoms_.size();
	}


	/// @brief Check if there is an atom_type associated with an atom_type_name string
	bool
	contains_atom_type( std::string const & atom_type_name ) const
	{
		auto
			iter( atom_type_index_.find( atom_type_name ) );
		return iter != atom_type_index_.end();
	}


	/// @brief Lookup the atom_type by the atom_type_name string
	int
	atom_type_index( std::string const & atom_type_name ) const
	{
		auto
			iter( atom_type_index_.find( atom_type_name ) );
		if ( iter == atom_type_index_.end() ) {
			utility_exit_with_message("unrecognized mm_atom_type_name "+atom_type_name);
		}
		return iter->second;
	}


	/// @brief Lookup an MMAtomType by 1-based indexing
	MMAtomType const &
	operator[] ( Size const index ) const
	{
		return *( atoms_[ index ] );
	}


	/// @brief Load the MMAtomTypeSet from a file
	void
	read_file( std::string const & filename );

	/// @brief Print all of the names of all of the MMAtomTypes in the set. Usefull for debuging.
	void
	print_all_types();


	// data
private:

	/// @brief What the ChemicalManager knows this as, if relevant
	std::string name_;

	/// @brief atom_type_index_ lookup map
	///
	/// @details atom_type_index_ allows lookup of the atom type index by a string
	std::map< std::string, int > atom_type_index_;

	/// @brief a collection of MMAtomTypes,
	///
	/// @details MMAtomType has data of atom properties, and it can be
	/// looked up by atom_type_index.
	utility::vector1< MMAtomTypeOP > atoms_;

};

#ifdef    SERIALIZATION
/// @brief Serialize an ElementSet - assumes that they're all stored in the ChemicalManager
template < class Archive >
void serialize_mm_atom_type_set( Archive & arc, MMAtomTypeSetCOP ptr );

/// @brief Deserialize an AtomTypeSet - assumes that they're all stored in the ChemicalManager
template < class Archive >
void deserialize_mm_atom_type_set( Archive & arc, MMAtomTypeSetCOP & ptr );
#endif // SERIALIZATION

} // chemical
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_MMAtomTypeSet )

SPECIAL_COP_SERIALIZATION_HANDLING( core::chemical::MMAtomTypeSet, core::chemical::serialize_mm_atom_type_set, core::chemical::deserialize_mm_atom_type_set )
#endif // SERIALIZATION

#endif // INCLUDED_core_chemical_MMAtomTypeSet_HH
