// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/Patch.hh
/// @author Phil Bradley

// See Patch.cc to understand what's going on.

#ifndef INCLUDED_core_chemical_Patch_hh
#define INCLUDED_core_chemical_Patch_hh

// Unit headers
#include <core/chemical/Patch.fwd.hh>

// Package headers
#include <core/chemical/PatchOperation.hh>
#include <core/chemical/ResidueTypeSelector.hh>

// Utility headers
#include <utility/vector1.hh>


namespace core {
namespace chemical {

/// @brief the string used to create new residue names after patching
extern std::string const PATCH_LINKER;

/// @brief helper function, returns the base residue name prior to any patching
std::string
residue_type_base_name( ResidueType const & rsd_type );

/// @brief helper function, returns the name of all added patches
std::string
residue_type_all_patches_name( ResidueType const & rsd_type );

///  @brief  A single case of a patch, eg proline Nterminus is a case of NtermProteinFull
class PatchCase : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~PatchCase();

	/// @brief whether the PatchCase is applicable to this ResidueType?
	bool
	applies_to( ResidueType const & rsd ) const
	{
		return selector_[ rsd ];
	}

	/// @brief returns patched residue, 0 if patch failed
	virtual
	ResidueTypeOP
	apply( ResidueType const & rsd_in, bool const instantiate = true ) const;

	/// @brief add one more operation in this PatchCase
	void
	add_operation( PatchOperationOP operation )
	{
		operations_.push_back( operation );
	}

	/// @brief to which ResidueTypes this PatchCase applies to?
	ResidueTypeSelector &
	selector()
	{
		return selector_;
	}

	/// @brief returns list of added atom names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	adds_atoms() const;

	/// @brief returns list of deleted atom names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	deletes_atoms() const;

	/// @brief returns list of added property names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	adds_properties() const;

	/// @brief returns list of deleted property names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	deletes_properties() const;

	/// @brief returns list of deleted variant names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	deletes_variants() const;

	// data:
private:
	/// @brief to which ResidueTypes this PatchCase applies to?
	ResidueTypeSelector selector_;
	/// @brief operations to done in this PatchCase
	utility::vector1< PatchOperationOP > operations_;
};


////////////////////////////////////////////////////////////////////////////////////////
/// @brief A class patching basic ResidueType to create variant types, containing multiple PatchCase
class Patch : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Patch();
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

	/// do I replace this residue type?
	virtual
	bool
	replaces( ResidueType const & rsd ) const
	{
		return 	applies_to( rsd ) && replaces_residue_type_;
	}

	/// @brief returns patched residue, 0 if patch failed
	virtual
	ResidueTypeOP
	apply( ResidueType const & rsd_type, bool const instantiate = true ) const;


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

	/// @brief returns list of added atom names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	adds_atoms( ResidueType const & rsd_in ) const;

	/// @brief returns list of deleted atom names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	deletes_atoms( ResidueType const & rsd_in ) const;

	/// @brief returns list of added property names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	adds_properties( ResidueType const & rsd_in ) const;

	/// @brief returns list of deleted property names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	deletes_properties( ResidueType const & rsd_in ) const;

	/// @brief returns list of deleted variant names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	deletes_variants( ResidueType const & rsd_in ) const;

	/// private data
private:
	/// name of the patch
	std::string name_;

	/// @brief variant types created by the patch
	utility::vector1< std::string > types_;

	/// @brief criteria to select ResidueTypes to which the patch is applied
	ResidueTypeSelector selector_;

	/// @brief different cases to which the patch is applied slightly differently, e.g., N-terminus patch to PRO and GLY
	utility::vector1< PatchCaseOP > cases_;

	/// @brief if set this patch will not change the name of the ResidueType and returns true for replaces()
	bool replaces_residue_type_;

};

} // chemical
} // core

#endif
