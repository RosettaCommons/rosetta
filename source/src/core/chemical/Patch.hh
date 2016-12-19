// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/Patch.hh
/// @author Phil Bradley

// See Patch.cc to understand what's going on.

#ifndef INCLUDED_core_chemical_Patch_hh
#define INCLUDED_core_chemical_Patch_hh

// Unit headers
#include <core/chemical/Patch.fwd.hh>

// Package headers
#include <core/chemical/PatchOperation.hh> // MSVC can't use fwd header
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

// Utility headers
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

/// @brief the string used to create new residue names after patching
std::string const PATCH_LINKER =  ":";

/// @brief handy function, return the first word from a line
std::string tag_from_line( std::string const & line );

/// @brief helper function, returns the base residue name prior to any patching
std::string
residue_type_base_name( ResidueType const & rsd_type );

/// @brief helper function, returns the name of all added patches
std::string
residue_type_all_patches_name( ResidueType const & rsd_type );

utility::vector1< std::string > get_patch_names( ResidueType const & rsd_type );

///  @brief  A single case of a patch, eg proline Nterminus is a case of NtermProteinFull
class PatchCase : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~PatchCase() override;

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

	void set_selector( ResidueTypeSelector const & selector ) { selector_ = selector; }

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

	/// @brief returns new name3, if changed
	std::string
	generates_new_name3() const;

	/// @brief returns interchangeability group, if set.
	std::string
	generates_interchangeability_group() const;

	bool
	may_change_aa() const;

	// data:
private:
	/// @brief to which ResidueTypes this PatchCase applies to?
	ResidueTypeSelector selector_;
	/// @brief operations to done in this PatchCase
	utility::vector1< PatchOperationOP > operations_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/// @brief create a PatchCase from input lines
/// @details add selector_ from lines enclosed by "BEGIN_SELECTOR" and "END_SELECTOR".\n
/// add operations_ from each input line containing a single operation
PatchCaseOP
case_from_lines(
	utility::vector1< std::string > const & lines,
	TypeSetMode res_type_set_mode = INVALID_t,
	std::string const & patch_name = ""
);

////////////////////////////////////////////////////////////////////////////////////////
/// @brief A class patching basic ResidueType to create variant types, containing multiple PatchCase
class Patch : public utility::pointer::ReferenceCount {
public:
	Patch() {}
	Patch( TypeSetMode res_type_set_mode ) : res_type_set_mode_(res_type_set_mode) {}

	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~Patch() override;
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
		return  applies_to( rsd ) && replaces_residue_type_;
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

	void set_selector( ResidueTypeSelector const & selector ) { selector_ = selector; }

	void set_name( std::string name ) { name_ = name; }

	/// @brief the variant types created by applying this patch
	virtual
	utility::vector1< std::string > const &
	types() const
	{
		return types_;
	}

	inline void types( utility::vector1< std::string > const & types ) {
		types_ = types;
	}

	inline
	void
	replaces_residue_type( bool replaces ) { replaces_residue_type_ = replaces; }

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

	inline void
	add_case( PatchCaseOP pcase ) { cases_.push_back( pcase ); }

	/// @brief returns new name3, if changed. Only one new name3 allowed.
	std::string
	generates_new_name3( ResidueType const & rsd_in ) const;

	/// @brief returns new interchangeability_group, if changed. Only one new interchangeability_group allowed.
	std::string
	generates_interchangeability_group( ResidueType const & rsd_in ) const;

	/// @brief returns new AA, if changed.
	chemical::AA
	generates_aa( ResidueType const & rsd_in ) const;

	/// private data
private:
	/// mode of the residuetypeset to which this patch belongs
	TypeSetMode res_type_set_mode_;

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

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // chemical
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_Patch )
#endif // SERIALIZATION


#endif
