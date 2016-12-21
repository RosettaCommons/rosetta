// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/ResidueTypeSet.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)
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

#include <basic/datacache/CacheableData.hh>

#include <map>
#include <set>

#ifdef    SERIALIZATION
#include <utility/serialization/serialization.hh>
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

/**
One thing that is not nailed down is whether a single ResidueTypeSet can have ResidueTypes with different AtomTypeSets. (PB-07/07)

The answer is NO -- or at least a ResidueTypeSet shouldn't have ResidueTypes with different mode() values,
which in most cases are determined by the AtomTypeSet being used. -
Though two different AtomTypeSets with the same mode may be permitted. (RM-11/16)
**/

/// @brief An abstract interface to a set of ResidueTypes
class ResidueTypeSet : public utility::pointer::enable_shared_from_this< ResidueTypeSet >
{

public:

	/// @brief default c-tor
	ResidueTypeSet( TypeSetMode mode = INVALID_t );

	~ResidueTypeSet();

	// Const pointers only - subclasses can implement the modifiable version.
	inline ResidueTypeSetCOP get_self_ptr() const { return shared_from_this(); }
	inline ResidueTypeSetCAP get_self_weak_ptr() const { return ResidueTypeSetCAP( shared_from_this() ); }

	AtomTypeSetCOP atom_type_set() const { return atom_types_; }
	ElementSetCOP element_set() const { return elements_; }
	MMAtomTypeSetCOP mm_atom_type_set() const { return mm_atom_types_; }
	orbitals::OrbitalTypeSetCOP orbital_type_set() const { return orbital_types_; }

protected: // We want resetting to go through the child classes

	void atom_type_set(AtomTypeSetCOP atom_types); // In cc file for error reporting.

	void element_set(ElementSetCOP elements) {
		elements_ = elements;
	}
	void mm_atom_type_set(MMAtomTypeSetCOP mm_atom_types) {
		mm_atom_types_ = mm_atom_types;
	}
	void orbital_type_set(orbitals::OrbitalTypeSetCOP orbital_types) {
		orbital_types_ = orbital_types;
	}

public:

	/// @brief The type of the ResidueTypeSet
	/// @details The difference between a ResidueTypeSet *name* and a ResidueTypeSet *mode* is that a
	/// a ResidueTypeSet *name* should uniquely identify a ResidueTypeSet (at lease those within the ChemicalManger)
	/// but more than one ResidueTypeSet may have the same *mode*.
	/// The mode specifies what compatibility class (full atom, centroid) the ResidueTypeSet has.
	TypeSetMode
	mode() const {
		return mode_;
	}

	/// @brief query ResidueType by its unique residue id.
	///
	/// @details since residue id is unique, it only returns
	/// one residue type or exit without match.
	virtual
	ResidueType const &
	name_map( std::string const & name ) const;

	/// @brief Get ResidueType by exact name, returning COP
	/// Will return null pointer for no matches
	virtual
	ResidueTypeCOP
	name_mapOP( std::string const & name ) const;

	/// @brief query if a ResidueType of the unique residue id (name) is present.
	virtual
	bool
	has_name( std::string const & name ) const = 0;

	/// @brief query if any ResidueTypes in the set have a "name3" tat matches the input name3
	bool has_name3( std::string const & name3 ) const;

	/// @brief Does this ResidueTypeSet have ResidueTypes with the given interchangeability group?
	bool
	has_interchangeability_group( std::string const & name ) const;

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
	virtual
	bool
	generates_patched_residue_type_with_name3( std::string const & base_residue_name, std::string const & name3 ) const;

	/// @brief Check if a base type (like "CYS") generates any types with a new interchangeability group (like "SCY" (via cys_acetylated))
	virtual
	bool
	generates_patched_residue_type_with_interchangeability_group( std::string const & base_residue_name, std::string const & interchangeability_group ) const;

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

	/// @brief accessor for merge behavior manager
	MergeBehaviorManager const &
	merge_behavior_manager() const;

	/// @brief The list of ResidueTypes that don't have any patches, but can be patched.
	virtual
	ResidueTypeCOPs base_residue_types() const { return base_residue_types_; }

	/// @brief The list of ResidueTypes which shouldn't get patches applied to them
	virtual
	ResidueTypeCOPs unpatchable_residue_types() const { return unpatchable_residue_types_; }

	/// @brief the patches
	virtual
	utility::vector1< PatchCOP > patches() const { return patches_; }

	/// @brief the metapatches
	virtual
	utility::vector1< MetapatchCOP > metapatches() const { return metapatches_; }

	/// @brief the patches, index by name.
	virtual
	std::map< std::string, utility::vector1< PatchCOP > > patch_map() const { return patch_map_; }

	/// @brief Do we have this metapatch?
	virtual
	bool
	has_metapatch( std::string const & name ) const {
		return metapatch_map_.find( name ) != metapatch_map_.end();
	}

	virtual
	MetapatchCOP
	metapatch( std::string const & name ) const {
		if ( ! has_metapatch( name ) ) {
			utility_exit_with_message(  "Metapatch " + name + " not in the metapatch map!" );
		}
		return metapatch_map_.find( name )->second;
	}

	/// @brief Gets all types with the given aa type and variants
	/// @details The number of variants must match exactly. Variants can be custom variants.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	virtual
	ResidueTypeCOPs
	get_all_types_with_variants_aa( AA aa, utility::vector1< std::string > const & variants ) const;

	/// @brief Gets all types with the given aa type and variants, making exceptions for some variants.
	/// @details The number of variants must match exactly. Variants can be custom variants, but exceptions must
	///           be standard types, listed in VariantType.hh.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	virtual
	ResidueTypeCOPs
	get_all_types_with_variants_aa( AA aa,
		utility::vector1< std::string > const & variants,
		utility::vector1< VariantType > const & exceptions ) const;


	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////

	//////////////////
	// protected methods

	//////////////////////////////////////////////////////////////////////////
	// These modification functions shouldn't be publically exposed in the
	// base class, as GlobalResidueTypeSets should have *no* public methods
	// which modify it, save for the initialization.
	//
	// Other subclasses might make these availible, though.
	/////////////////////////////////////////////////////////////////////////

protected:

	/// @brief Centralize the steps for preparing the ResidueType for addition to the RTS
	/// (e.g. if there's any additional modifications that need to get done.)
	void
	prep_restype( ResidueTypeOP new_type );

	/// @brief adds a new base residue type to the set, one that isn't patched, but can be.
	virtual
	void
	add_base_residue_type( ResidueTypeOP new_type );

	/// @brief adds a new residue type to the set, one that can be patched
	virtual
	void
	add_base_residue_type( std::string const &  filename );

	/// @brief adds new residue types, ones that can be patched
	virtual
	void
	read_files_for_base_residue_types(
		utility::vector1< std::string > const & filenames
	);

	/// @brief adds a new residue type to the set, one that CANNOT be generated from a base_residue_type and patches, and shouldn't have patches applied
	virtual
	void
	add_unpatchable_residue_type( ResidueTypeOP new_type );

	/// @brief adds a new residue type to the set, one that CANNOT be generated from a base_residue_type and patches, and shouldn't have patches applied
	virtual
	void
	add_unpatchable_residue_type( std::string const &  filename );

	/// @brief adds new residue types, ones that CANNOT be generated from a base_residue_type and patches, and shouldn't have patches applied
	virtual
	void
	read_files_for_unpatchable_residue_types(
		utility::vector1< std::string > const & filenames
	);

	/// @brief delete an base residue type from the set (Use with care)
	/// Currently will not remove any patched types
	virtual
	void
	remove_base_residue_type( std::string const & name );

	/// @brief delete an unpatchable residue type from the set (Use with care)
	virtual
	void
	remove_unpatchable_residue_type( std::string const & name );

protected:

	/// @brief helper function used during replacing residue types after, e.g., orbitals.
	bool
	update_base_residue_types_if_replaced( ResidueTypeCOP rsd_type, ResidueTypeCOP rsd_type_new );

	void
	mode( TypeSetMode setting ) { mode_ = setting; }

	// The alterable cache object
	ResidueTypeSetCacheOP
	cache_object() const { return cache_; }

	virtual
	bool
	generate_residue_type( std::string const & rsd_name ) const;

	/// @brief Attempt to lazily load the given residue type from data.
	virtual
	bool
	lazy_load_base_type( std::string const & rsd_base_name ) const = 0;

	void
	figure_out_last_patch_from_name( std::string const & rsd_name,
		std::string & rsd_name_base,
		std::string & patch_name ) const;

	virtual
	void
	add_patches(
		utility::vector1< std::string > const & patch_filenames,
		utility::vector1< std::string > const & metapatch_filenames
	);

	virtual
	void
	set_merge_behavior_manager( MergeBehaviorManagerCOP mbm);

	/// @brief A list of L-chirality base types with an equivalent D-chirality base type.
	std::map < ResidueTypeCOP /*L-type*/, ResidueTypeCOP /*D-type*/> &
	l_to_d_mapping() { return l_to_d_mapping_; }

	/// @brief A list of D-chirality base types with an equivalent L-chirality base type.
	std::map < ResidueTypeCOP /*D-type*/, ResidueTypeCOP /*L-type*/> &
	d_to_l_mapping() { return d_to_l_mapping_; }

protected:
	// External code shouldn't be able to copy a ResidueTypeSet, although
	// derived types may expose it.
	ResidueTypeSet( ResidueTypeSet const & );
	ResidueTypeSet & operator = ( ResidueTypeSet const & ) = delete;

#ifdef    SERIALIZATION
public:
	/// @brief Does this ResidueTypeSet *directly* have the ResidueType?
	bool
	has( ResidueTypeCOP restype ) const;
#endif // SERIALIZATION

	//////////////////
	// data
private:

	// The default subsidiary typesets, typically specified in the database summary file.
	// You can add a residue type with a different subsidiary typeset, but you'll have to
	// construct it yourself.
	AtomTypeSetCOP atom_types_;
	ElementSetCOP elements_;
	MMAtomTypeSetCOP mm_atom_types_;
	orbitals::OrbitalTypeSetCOP orbital_types_;

	/// @brief What sort of TypeSet is this?
	TypeSetMode mode_;

	MergeBehaviorManagerCOP merge_behavior_manager_;

	/// @brief all cached residue_type information including generated residue_types, name3_map, etc.
	/// By making the following an OP (instead of COP) the cache effectively becomes mutable even when in a
	/// const ResidueTypeSet.
	ResidueTypeSetCacheOP cache_;

	/// @brief ResidueTypes that don't have any patches, but can be patched.
	/// @details The ResidueTypes represented here are also present in the ResidueTypeSetCache
	ResidueTypeCOPs base_residue_types_;

	/// @brief ResidueTypes which shouldn't get patches applied to them
	/// @details The ResidueTypes represented here are also present in the ResidueTypeSetCache
	ResidueTypeCOPs unpatchable_residue_types_;

	/// @brief the patches
	utility::vector1< PatchCOP > patches_;
	utility::vector1< MetapatchCOP > metapatches_;

	/// @brief patches indexed by name
	std::map< std::string, utility::vector1< PatchCOP > > patch_map_;
	std::map< std::string, MetapatchCOP > metapatch_map_;

	/// @brief A list of L-chirality base types with an equivalent D-chirality base type.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	std::map < ResidueTypeCOP /*L-type*/, ResidueTypeCOP /*D-type*/> l_to_d_mapping_;

	/// @brief A list of D-chirality base types with an equivalent L-chirality base type.
	/// @details For reverse searches.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	std::map < ResidueTypeCOP /*D-type*/, ResidueTypeCOP /*L-type*/> d_to_l_mapping_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // chemical
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_ResidueTypeSet )

#include <core/chemical/ResidueTypeSet.srlz.hh>
SPECIAL_COP_SERIALIZATION_HANDLING( core::chemical::ResidueTypeSet, core::chemical::serialize_residue_type_set, core::chemical::deserialize_residue_type_set )

#endif // SERIALIZATION

#endif
