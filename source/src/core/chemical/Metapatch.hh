// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/Metapatch.hh
/// @author Andy watkins

#ifndef INCLUDED_core_chemical_Metapatch_hh
#define INCLUDED_core_chemical_Metapatch_hh

// Unit headers
#include <core/chemical/Metapatch.fwd.hh>
#include <core/chemical/Patch.hh>

// Package headers
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/chemical/AtomProperty.hh>

// Utility headers
#include <utility/vector1.hh>


namespace core {
namespace chemical {

////////////////////////////////////////////////////////////////////////////////////////
/// @brief A class patching basic ResidueType to create variant types, containing multiple PatchCase
class Metapatch : public utility::pointer::ReferenceCount {
public:

	Metapatch();

	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Metapatch();
	/// @brief constructor from file
	void
	read_file( std::string const & filename );

	/// can I operate on this residue type?
	virtual
	bool
	applies_to( ResidueType const & rsd ) const
	{
		return selector_[ rsd ];
	}

	PatchCOP
	get_one_patch( /*ResidueType const & rsd_type, */std::string const & atom_name ) const;

	inline
	bool
	meets_requirements( ResidueType const & r, Size i ) const {
		if ( pertinent_property_ == NO_ATOM_PROPERTY ) return true;
		return ( r.atom(i).has_property( pertinent_property_ ) );
	}

	/// @brief unique name of this patch, eg Nter-simple, Cter-full, Phospho, ... ?
	virtual
	std::string const &
	name() const
	{
		return name_;
	}

	/// @brief the variant types created by applying this patch
	virtual
	utility::vector1< std::string > const &
	types() const
	{
		return types_;
	}

	utility::vector1< std::string >
	atoms( ResidueType const & rsd_type ) const;

	/// private data
private:
	/// name of the patch
	std::string name_;

	/// @brief variant types created by the created patches
	utility::vector1< std::string > types_;

	/// @brief different cases to which the derived patches patch will be applied slightly differently, e.g., N-terminus patch to PRO and GLY
	utility::vector1< PatchCaseOP > cases_;

	chemical::AtomProperty pertinent_property_;

	utility::vector1< std::string > case_lines_;

	/// @brief criteria to select ResidueTypes to which the patches are applied
	ResidueTypeSelector selector_;

};

} // chemical
} // core

#endif
