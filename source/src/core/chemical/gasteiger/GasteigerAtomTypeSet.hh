// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/gasteiger/GasteigerAtomTypeSet.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_gasteiger_GasteigerAtomTypeSet_hh
#define INCLUDED_core_chemical_gasteiger_GasteigerAtomTypeSet_hh


// Unit headers
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.fwd.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.fwd.hh>
#include <core/chemical/ElementSet.hh>

// Project headers

// Utility headers
#include <utility/exit.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers

#include <map>

#include <utility/vector1.hh>

#include <core/types.hh>

namespace core {
namespace chemical {
namespace gasteiger {

/// @brief A set of Bcl Atom types
///
/// @details This class contains a vector of pointers each of which points to an
/// Atom and the vector index is looked up by an element_name string
/// in a map.
///
class GasteigerAtomTypeSet : public utility::pointer::ReferenceCount {

private:

	GasteigerAtomTypeSet(); // private default constructor - Require element type set when instantiated.

public:

	GasteigerAtomTypeSet( ElementSetCAP element_set );
	GasteigerAtomTypeSet( GasteigerAtomTypeSet const & other );
	virtual ~GasteigerAtomTypeSet();

	/// @brief Number of atom types in the set
	Size
	n_types() const
	{
		return atom_types_.size();
	}

	/// @brief Return the associated element type set.
	ElementSetCAP element_set() const;

	/// @brief Check if there is an element_type associated with an element_symbol string
	bool
	contains_atom_type( std::string const & atom_type_name ) const;

	/// @brief Lookup the element index by the element_symbol string
	Size
	atom_type_index( std::string const & atom_type_name ) const;

	/// @brief Lookup the element index by the element_symbol string
	GasteigerAtomTypeDataCOP
	atom_type( std::string const & atom_type_name ) const;

	/// @brief Lookup an Atom by 1-based indexing
	GasteigerAtomTypeDataCOP
	operator[] ( Size const index ) const;


	/// @brief Load the AtomSet from a file
	void
	read_file( std::string const & filename );

	// data
private:

	/// @brief The associated element type set.
	ElementSetCAP element_set_;

	/// @brief allows lookup of the atom type by its name
	std::map< std::string, core::Size > atom_type_index_;

	/// @brief the set of atom types,
	utility::vector1< GasteigerAtomTypeDataOP > atom_types_;

};

} // gasteiger
} // chemical
} // core

#endif
