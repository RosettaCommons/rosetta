// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/IdealBondLengthSet.hh
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


#ifndef INCLUDED_core_chemical_IdealBondLengthSet_hh
#define INCLUDED_core_chemical_IdealBondLengthSet_hh


// Unit headers
#include <core/chemical/IdealBondLengthSet.fwd.hh>

// Project headers

// Utility headers

#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>

#include <utility/vector1_bool.hh>

#include <core/types.hh>
#include <core/chemical/Element.fwd.hh>


namespace core {
namespace chemical {

typedef int AtomTypeIndex;
typedef Real BondLength;

/// @brief A set of Elements
///
/// @details This class contains a vector of pointers each of which points to an
/// Element and the vector index is looked up by an element_name string
/// in a map.
///
class IdealBondLengthSet : public utility::pointer::ReferenceCount {

public:
	IdealBondLengthSet();
	virtual ~IdealBondLengthSet();

	/// @brief Check if an ideal bond length is known for this pair of atom types...
	bool contains_bond_length( std::string const & atom_type_name1, std::string const & atom_type_name2) const;
	bool contains_bond_length( AtomTypeIndex atom_type_index1, AtomTypeIndex atom_type_index2) const;

	/// @brief Lookup the element index by the element_symbol string
	BondLength get_bond_length( std::string const & atom_type_name1, std::string const & atom_type_name2) const;
	BondLength get_bond_length( AtomTypeIndex const atom_type_index1, AtomTypeIndex const atom_type_index2) const;

	/// @brief Load the IdealBondLengthSet from a file
	void
	read_file( std::string const & filename );

	/// @brief Print all of the symbols of all of the Elements in the set. Usefull for debuging.
	void
	print_all_bond_lengths();


	// data
private:
	void
	add_bond_length(
		std::string const & atom_type_name1,
		std::string const & atom_type_name2,
		BondLength const length
	);

	void
	add_bond_length(
		AtomTypeIndex const atom_type_index1,
		AtomTypeIndex const atom_type_index2,
		Real const length
	);

	/// @details pair of atom type indices and the ideal length between these atom types
	std::map< std::pair<AtomTypeIndex, AtomTypeIndex>, BondLength > bond_lengths_;

};

} // chemical
} // core

#endif // INCLUDED_core_chemical_IdealBondLengthSet_HH
