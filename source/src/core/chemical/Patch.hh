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
#include <core/chemical/VariantType.hh>
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

	/// @brief Returns list of added property enums.  Useful for identifying patches that go with PDB residues.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	utility::vector1< ResidueProperty >
	adds_properties_enums() const;

	/// @brief returns list of deleted property names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	deletes_properties() const;

	/// @brief Returns list of deleted property enums.  Useful for identifying patches that go with PDB residues.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	utility::vector1< ResidueProperty >
	deletes_properties_enums() const;

	/// @brief returns list of deleted variant names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	deletes_variants() const;

	/// @brief Returns list of deleted VariantTypes.  Doesn't support on-the-fly VariantTypes.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	utility::vector1< core::chemical::VariantType >
	deletes_variants_by_enum() const;

	/// @brief returns new name3, if changed
	std::string
	generates_new_name3() const;

	/// @brief returns interchangeability group, if set.
	std::string
	generates_interchangeability_group() const;

	bool
	may_change_aa() const;

	/// @brief Can the patch case change the connections for atom on the ResidueType?
	bool
	changes_connections_on( ResidueType const & rsd_in, std::string const & atom ) const;

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

	/// @brief Constructor
	///
	Patch();

	/// @brief Constructor with ResidueTypeSet mode.
	///
	Patch( TypeSetMode res_type_set_mode );

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

	/// @brief Returns the name of the residueType after applying the patch
	/// (Theoretical - doesn't actually check if the patch would be applied.)
	std::string
	patched_name( ResidueType const & rsd ) const;

	void set_selector( ResidueTypeSelector const & selector ) { selector_ = selector; }

	void set_name( std::string const & name ) {
		name_ = name;
		// Ugh.  String-parsing.  Figure this out ONCE, dammit!  VKM.
		// This *is* the once. If you're allowing the name to be re-set, you need to
		// re-determine if it's a metapatch. Note that this function is only called
		// during construction so in any particular patch it's only called once.
		is_metapatch_ = (name_.substr( 0, 3 ) == "MP-");
	}

	/// @brief The variant types created by applying this patch
	/// @brief Use custom_types() for a vector of strings of custom types.
	virtual
	utility::vector1< core::chemical::VariantType > const &
	types() const
	{
		return types_;
	}

	/// @brief The custom variant types added by this patch.
	/// @details This must be a vector of strings, since the custom types are generated at runtime and can't be
	/// enumerated at compile time.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	virtual
	utility::vector1 < std::string > const &
	custom_types() const {
		return custom_types_;
	}

	/// @brief Set the variant types created by applying this patch.
	/// @details For custom variant types, use the add_custom_type() function, which takes a string.
	inline void types( utility::vector1< core::chemical::VariantType > const & types ) {
		types_ = types;
	}

	/// @brief Add a VariantType to the list that this patch applies.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	virtual
	void add_type( core::chemical::VariantType const type );

	/// @brief Add a custom VariantType to the list that this patch applies.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	virtual
	void add_custom_type( std::string const & custom_type );

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

	/// @brief Returns list of added properties, by enum.  Useful for identifying patches that go with PDB residues.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	utility::vector1< ResidueProperty >
	adds_properties_enums( ResidueType const & rsd_in ) const;

	/// @brief returns list of deleted property names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	deletes_properties( ResidueType const & rsd_in ) const;

	/// @brief returns list of deleted property enums, useful for identifying patches that go with PDB residues.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	utility::vector1< ResidueProperty >
	deletes_properties_enums( ResidueType const & rsd_in ) const;

	/// @brief returns list of deleted variant names, useful for identifying patches that go with PDB residues
	utility::vector1< std::string >
	deletes_variants( ResidueType const & rsd_in ) const;

	/// @brief Returns list of deleted VariantTypes.  Doesn't support on-the-fly VariantTypes.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	utility::vector1< core::chemical::VariantType > deletes_variants_by_enum( ResidueType const & rsd_type ) const;

	/// @brief Does the patch potentially change the connections for the given atom on the ResidueType
	bool
	changes_connections_on( ResidueType const & rsd_in, std::string const & atom ) const;

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

	/// @brief Is this a metapatch?
	///
	inline bool is_metapatch() const { return is_metapatch_; }

	/// private data
private:
	/// mode of the residuetypeset to which this patch belongs
	TypeSetMode res_type_set_mode_ = core::chemical::FULL_ATOM_t;

	/// name of the patch
	std::string name_ = "";

	/// Is this a metapatch?
	bool is_metapatch_ = false;

	/// @brief Variant types created by the patch.
	utility::vector1< core::chemical::VariantType > types_;

	/// @brief Custom variant types created by this patch.
	utility::vector1 < std::string > custom_types_;

	/// @brief criteria to select ResidueTypes to which the patch is applied
	ResidueTypeSelector selector_;

	/// @brief different cases to which the patch is applied slightly differently, e.g., N-terminus patch to PRO and GLY
	utility::vector1< PatchCaseOP > cases_;

	/// @brief if set this patch will not change the name of the ResidueType and returns true for replaces()
	bool replaces_residue_type_ = false;

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
