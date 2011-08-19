// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/CSDAtomTypeSet.hh
/// @author Ian W. Davis


#ifndef INCLUDED_core_chemical_CSDAtomTypeSet_hh
#define INCLUDED_core_chemical_CSDAtomTypeSet_hh


// Unit headers
#include <core/chemical/CSDAtomTypeSet.fwd.hh>

// Project headers
#include <core/chemical/CSDAtomType.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>


namespace core {
namespace chemical {


/// @brief A set of CSDAtomTypes
///
/// @details This class contains a vector of pointers each of
/// which points to an CSDAtomType and the vector index is looked
/// up by an atom_name string in a map.
class CSDAtomTypeSet : public utility::pointer::ReferenceCount {

public:
	CSDAtomTypeSet();
	~CSDAtomTypeSet();

	/// @brief Number of atom types in the set
	Size
	n_atomtypes() const;

	/// @brief Check if there is an atom_type associated with an
	/// atom_type_name string
	bool
	contains_atom_type( std::string const & atom_type_name ) const;

	/// @brief Lookup the atom_type by the atom_type_name string
	int
	atom_type_index( std::string const & atom_type_name ) const;

	/// @brief Lookup an CSDAtomType by 1-based indexing
	CSDAtomType const &
	operator[] ( Size const index ) const;

	/// @brief Load the CSDAtomTypeSet from a file
	void
	read_file( std::string const & filename );

	/// @brief Print all of the names of all of the CSDAtomTypes in
	/// the set. Useful for debugging.
	void
	print_all_types();

	// data
private:

	/// @brief atom_type_index_ lookup map
	///
	/// @details atom_type_index_ allows lookup of the atom type
	/// index by a string
	std::map< std::string, int > atom_type_index_;

	/// @brief a collection of CSDAtomTypes,
	///
	/// @details CSDAtomType has data of atom properties, and it can be
	/// looked up by atom_type_index.
	utility::vector1< CSDAtomType* > atoms_;

};

} // chemical
} // core

#endif // INCLUDED_core_chemical_CSDAtomTypeSet_HH
