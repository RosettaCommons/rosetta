// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/ResidueTypeSet.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_ResidueTypeSet_hh
#define INCLUDED_core_chemical_ResidueTypeSet_hh


// Unit headers
#include <core/chemical/ResidueTypeSet.fwd.hh>

// Package headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/AA.hh>
//#include <core/chemical/ResidueTypeSelector.fwd.hh>
#include <core/chemical/ResidueTypeSetCache.fwd.hh>
#include <core/chemical/MergeBehaviorManager.fwd.hh>

// STL headers
#include <list>

#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Metapatch.fwd.hh>
#include <core/chemical/Patch.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>
#include <set>

//Auto Headers
//#include <core/chemical/Adduct.fwd.hh>
//
//#ifdef WIN32
//#include <core/chemical/Adduct.hh>
//#endif

namespace core {
namespace chemical {

///  A collection of ResidueType defined
/**
One thing that is not nailed down is whether a single ResidueSet can have ResidueType's with
different AtomTypeSets. I've left open this possibility currently although there isnt any code
for this yet (PB-07/07)
**/


class ResidueTypeSet : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< ResidueTypeSet >
{
public:
	typedef std::list< AA >::const_iterator AAsIter;
	typedef std::map< std::string, ResidueTypeCOP >::const_iterator const_residue_iterator;

public:


	/// @brief default c-tor
	ResidueTypeSet( TypeSetCategory category = INVALID_t );

	/// @brief constructor from directory
	ResidueTypeSet(
		std::string const & name,
		std::string const & directory
	);

	~ResidueTypeSet() override;

	void init(
		std::vector< std::string > const & extra_res_param_files = std::vector< std::string >(),
		std::vector< std::string > const & extra_patch_files = std::vector< std::string >()
	);

	/// self pointers
	inline ResidueTypeSetCOP get_self_ptr() const { return shared_from_this(); }
	inline ResidueTypeSetOP  get_self_ptr() { return shared_from_this(); }
	inline ResidueTypeSetCAP get_self_weak_ptr() const { return ResidueTypeSetCAP( shared_from_this() ); }
	//inline ResidueTypeSetAP  get_self_weak_ptr() { return ResidueTypeSetAP( shared_from_this() ); }

	/// @brief name of the residue type set
	/// @details The difference between a ResidueTypeSet *name* and a ResidueTypeSet *category* is that a
	/// a ResidueTypeSet *name* should uniquely identify a ResidueTypeSet (at lease those within the ChemicalManger)
	/// but more than one ResidueTypeSet may have the same *category*.
	/// The type specifies what compatibility class (full atom, centroid) the ResidueTypeSet has.
	/// Generally speaking, the *name* should only be used when interacting with the user.
	std::string const &
	name() const {
		return name_;
	}

	/// @brief The type of the ResidueTypeSet
	/// @details The difference between a ResidueTypeSet *name* and a ResidueTypeSet *category* is that a
	/// a ResidueTypeSet *name* should uniquely identify a ResidueTypeSet (at lease those within the ChemicalManger)
	/// but more than one ResidueTypeSet may have the same *category*.
	/// The type specifies what compatibility class (full atom, centroid) the ResidueTypeSet has.
	TypeSetCategory
	category() const {
		return category_;
	}

	AtomTypeSetCOP atom_type_set() const { return atom_types_; }
	ElementSetCOP element_set() const { return elements_; }
	MMAtomTypeSetCOP mm_atom_type_set() const { return mm_atom_types_; }
	orbitals::OrbitalTypeSetCOP orbital_type_set() const { return orbital_types_; }

	void atom_type_set(AtomTypeSetCOP atom_types) {
		runtime_assert( !atom_types_ ); // Don't change a set default.
		atom_types_ = atom_types;
	}
	void element_set(ElementSetCOP elements) {
		runtime_assert( !elements_ ); // Don't change a set default.
		elements_ = elements;
	}
	void mm_atom_type_set(MMAtomTypeSetCOP mm_atom_types) {
		runtime_assert( !mm_atom_types_ ); // Don't change a set default.
		mm_atom_types_ = mm_atom_types;
	}
	void orbital_type_set(orbitals::OrbitalTypeSetCOP orbital_types) {
		runtime_assert( !orbital_types_ ); // Don't change a set default.
		orbital_types_ = orbital_types;
	}

	/// @brief adds a new residue type to the set, one that CANNOT be generated from a base_residue_type and patches
	void
	add_custom_residue_type( ResidueTypeOP new_type );

	/// @brief adds a new residue type to the set, one that CANNOT be generated from a base_residue_type and patches
	void
	add_custom_residue_type( std::string const &  filename );

	/// @brief delete a custom residue type from the set (Use with care)
	void
	apply_patches(
		utility::vector1< std::string > const & patch_filenames,
		utility::vector1< std::string > const & metapatch_filenames
	);

	/// @brief delete a custom residue type from the set (Use with care)
	void
	remove_custom_residue_type( std::string const & name );

	/// @brief adds new residue types, ones that CANNOT be generated from a base_residue_type and patches
	void
	read_files_for_custom_residue_types(
		utility::vector1< std::string > const & filenames
	);

	/// @brief delete a custom residue type from the set (Use with *extreme* care)
	void
	remove_base_residue_type_DO_NOT_USE( std::string const & name );

	/// @brief query ResidueType by its unique residue id.
	///
	/// @details since residue id is unique, it only returns
	/// one residue type or exit without match.
	virtual
	ResidueType const &
	name_map( std::string const & name ) const;

	/// @brief Get ResidueType by exact name, returning COP
	/// Will return null pointer for no matches
	ResidueTypeCOP
	name_mapOP( std::string const & name ) const;

	/// @brief Get the base ResidueType with the given aa type and variants
	/// @details Returns 0 if one does not exist.
	/// The returned type will have at least all the variants given, but may have more
	/// if a minimal variant type isn't availible.
	ResidueTypeCOP
	get_representative_type_aa( AA aa,
		utility::vector1< std::string > const & variants ) const;

	ResidueTypeCOP
	get_representative_type_aa( AA aa ) const;

	/// @brief Get the base ResidueType with the given name1 and variants
	/// @details Returns 0 if one does not exist.
	/// The returned type will have at least all the variants given, but may have more
	/// if a minimal variant type isn't availible.
	ResidueTypeCOP
	get_representative_type_name1( char name1,
		utility::vector1< std::string > const & variants ) const;

	ResidueTypeCOP
	get_representative_type_name1( char name1 ) const;

	/// @brief Get the base ResidueType with the given name3 and variants
	/// @details Returns 0 if one does not exist.
	/// The returned type will have at least all the variants given, but may have more
	/// if a minimal variant type isn't availible.
	ResidueTypeCOP
	get_representative_type_name3( std::string const &  name3,
		utility::vector1< std::string > const & variants ) const;

	ResidueTypeCOP
	get_representative_type_name3( std::string const &  name3 ) const;

	ResidueTypeCOP
	get_representative_type_base_name( std::string const & base_name ) const;
	
	/// @brief Gets all non-patched types with the given aa type
	ResidueTypeCOPs
	get_base_types_aa( AA aa ) const;

	/// @brief Get all non-patched ResidueTypes with the given name1
	ResidueTypeCOPs
	get_base_types_name1( char name1 ) const;

	/// @brief Get all non-patched ResidueTypes with the given name3
	ResidueTypeCOPs
	get_base_types_name3( std::string const &  name3 ) const;

	/// @brief Given a D-residue, get its L-equivalent.
	/// @details Returns NULL if there is no equivalent, true otherwise.  Throws an error if this is not a D-residue.
	/// Preserves variant types.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	ResidueTypeCOP get_d_equivalent( ResidueTypeCOP l_rsd ) const;

	/// @brief Given an L-residue, get its D-equivalent.
	/// @details Returns NULL if there is no equivalent, true otherwise.  Throws an error if this is not an L-residue.
	/// Preserves variant types.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	ResidueTypeCOP get_l_equivalent( ResidueTypeCOP d_rsd ) const;

	/// @brief Given a residue, get its mirror-image type.
	/// @details Returns the same residue if this is an ACHIRAL type (e.g. gly), the D-equivalent for an L-residue, the L-equivalent of a D-residue,
	/// or NULL if this is an L-residue with no D-equivalent (or a D- with no L-equivalent).  Preserves variant types.
	ResidueTypeCOP get_mirrored_type( ResidueTypeCOP original_rsd ) const;

	/// @brief Check if a base type (like "SER") generates any types with another name3 (like "SEP")
	bool
	generates_patched_residue_type_with_name3( std::string const & base_residue_name, std::string const & name3 ) const;

	/// @brief Check if a base type (like "CYS") generates any types with a new interchangeability group (like "SCY" (via cys_acetylated))
	bool
	generates_patched_residue_type_with_interchangeability_group( std::string const & base_residue_name, std::string const & interchangeability_group ) const;

	/// @brief Gets all types with the given aa type and variants
	/// @details The number of variants must match exactly. Variants can be custom variants.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	ResidueTypeCOPs
	get_all_types_with_variants_aa( AA aa, utility::vector1< std::string > const & variants ) const;

	/// @brief Gets all types with the given aa type and variants, making exceptions for some variants.
	/// @details The number of variants must match exactly. Variants can be custom variants, but exceptions must
	///           be standard types, listed in VariantType.hh.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	ResidueTypeCOPs
	get_all_types_with_variants_aa(
		AA aa,
		utility::vector1< std::string > const & variants,
		utility::vector1< VariantType > const & exceptions
	) const;

	/// @brief Get all non-patched ResidueTypes with the given name1
	/// @details The number of variants must match exactly.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	ResidueTypeCOPs
	get_all_types_with_variants_name1( char name1, utility::vector1< std::string > const & variants ) const;

	/// @brief Get all non-patched ResidueTypes with the given name3
	/// @details The number of variants must match exactly.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	ResidueTypeCOPs
	get_all_types_with_variants_name3(
		std::string const & name3,
		utility::vector1< std::string > const & variants
	) const;

	/// @brief query if a ResidueType of the unique residue id (name) is present.
	bool has_name( std::string const & name ) const;

	/// @brief query if any ResidueTypes in the set have a "name3" tat matches the input name3
	bool has_name3( std::string const & name3 ) const;

	/// @brief Does this ResidueTypeSet have ResidueTypes with the given interchangeability group?
	bool
	has_interchangeability_group( std::string const & name ) const;

	/// @brief Query a variant ResidueType by its base ResidueType and VariantType
	ResidueType const &
	get_residue_type_with_variant_added(
		ResidueType const & init_rsd,
		VariantType const new_type
	) const;

	/// @brief return the residuetype we get from variant rsd type after removing the desired variant type
	ResidueType const &
	get_residue_type_with_variant_removed(
		ResidueType const & init_rsd,
		VariantType const old_type
	) const;

	/// @brief accessor for database_directory
	std::string const &
	database_directory() const;

	/// @brief accessor for merge behavior manager
	MergeBehaviorManager const &
	merge_behavior_manager() const;

	/// @brief the residues with no patches applied
	ResidueTypeCOPs base_residue_types() const { return base_residue_types_; }

	/// @brief the residues with no patches applied
	ResidueTypeCOPs custom_residue_types() const { return custom_residue_types_; }

	/// @brief the patches
	utility::vector1< PatchCOP > const & patches() const { return patches_; }

	/// @brief the metapatches
	utility::vector1< MetapatchCOP > const & metapatches() const { return metapatches_; }

	/// @brief the patches, index by name.
	std::map< std::string, utility::vector1< PatchCOP > > const & patch_map() const { return patch_map_; }

	MetapatchCOP
	metapatch( std::string name ) const { return metapatch_map_.find( name )->second; }

	//////////////////
	// private methods
private:

	/// @brief read a list of residue types
	void
	read_list_of_residues(
		std::string const & list_filename
	);

	void
	read_files(
		utility::vector1< std::string > const & filenames
	);

	void
	update_info_on_name3_and_interchangeability_group( ResidueTypeCOPs base_residue_types );

	void
	generate_all_residue_types();

	bool
	generate_residue_type( std::string const & rsd_name ) const;

	void
	figure_out_last_patch_from_name( std::string const & rsd_name,
		std::string & rsd_name_base,
		std::string & patch_name ) const;

	/// @brief helper function used during replacing residue types after, e.g., orbitals.
	bool
	update_base_residue_types_if_replaced( ResidueTypeCOP rsd_type, ResidueTypeCOP rsd_type_new );

private:

	/// @brief From a file, read which IDs shouldn't be loaded from the components.
	void
	load_shadowed_ids( std::string const & directory, std::string const & file = "shadow_list.txt" );

	/// @brief Attempt to lazily load the given residue type from data.
	bool
	lazy_load_base_type( std::string const & rsd_base_name ) const;

	/// @brief Load a residue type from the components dictionary.
	ResidueTypeOP
	load_pdb_component( std::string const & pdb_id ) const;

	//////////////////
	// data
private:

	/// @brief What does the ChemicalManager call this ResidueTypeSet?
	std::string name_;

	/// @brief What sort of TypeSet is this?
	TypeSetCategory category_;

	// The default subsidiary typesets, typically specified in the database summary file.
	// You can add a residue type with a different subsidiary typeset, but you'll have to
	// construct it yourself.
	AtomTypeSetCOP atom_types_;
	ElementSetCOP elements_;
	MMAtomTypeSetCOP mm_atom_types_;
	orbitals::OrbitalTypeSetCOP orbital_types_;

	/// @brief residue types with no patches applied, read in from database.
	ResidueTypeCOPs base_residue_types_;

	/// @brief information on residue types whose name3's can be changed by patches.
	std::map< std::string, std::set< std::string > > name3_generated_by_base_residue_name_;
	/// @brief interchangeability groups that appear upon patch application.
	std::map< std::string, std::set< std::string > > interchangeability_group_generated_by_base_residue_name_;

	/// @brief new residue types added at runtime with, e.g., custom variants.
	ResidueTypeCOPs custom_residue_types_;

	/// @brief the database directory of the generating files ---> allows to use cached dunbrack libs
	const std::string database_directory_;

	MergeBehaviorManagerCOP merge_behavior_manager_;

	/// @brief the patches
	utility::vector1< PatchCOP > patches_;
	utility::vector1< MetapatchCOP > metapatches_;

	/// @brief patches indexed by name
	std::map< std::string, utility::vector1< PatchCOP > > patch_map_;
	std::map< std::string, MetapatchCOP > metapatch_map_;

	/// @brief Which components shouldn't be loaded from the components file.
	std::set< std::string > shadowed_ids_;

	/// @brief data for lazy loading of PDB components
	std::string pdb_components_filename_;

	/// @brief all cached residue_type information including generated residue_types, name3_map, etc.
	/// By making the following an OP (instead of COP) the cache effectively becomes mutable even when in a
	/// const ResidueTypeSet.
	ResidueTypeSetCacheOP cache_;

	/// @brief A list of L-chirality base types with an equivalent D-chirality base type.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	std::map < ResidueTypeCOP /*L-type*/, ResidueTypeCOP /*D-type*/> l_to_d_mapping_;

	/// @brief A list of D-chirality base types with an equivalent L-chirality base type.
	/// @details For reverse searches.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	std::map < ResidueTypeCOP /*D-type*/, ResidueTypeCOP /*L-type*/> d_to_l_mapping_;

private:
	// uncopyable
	ResidueTypeSet( ResidueTypeSet const & );
	ResidueTypeSet const & operator = ( ResidueTypeSet const & );
};

} // chemical
} // core


#endif
