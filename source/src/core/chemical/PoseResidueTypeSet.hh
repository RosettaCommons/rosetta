// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core/chemical/PoseResidueTypeSet.hh
/// @brief A ResidueTypeSet which can be cached in the Pose
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_PoseResidueTypeSet_hh
#define INCLUDED_core_chemical_PoseResidueTypeSet_hh

#include <core/chemical/PoseResidueTypeSet.fwd.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ header
#include <string>
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

///@brief A ResidueTypeSet which can be cached in the Pose
class PoseResidueTypeSet : public core::chemical::ResidueTypeSet {

public:

	PoseResidueTypeSet();
	PoseResidueTypeSet(PoseResidueTypeSet const & src);
	PoseResidueTypeSet(ResidueTypeSetCOP deflt_rts);

	// Note: this is *not* from the base class.
	virtual
	PoseResidueTypeSetOP
	clone() const;

	virtual ~PoseResidueTypeSet();

	/// @brief Set the default fall-back ResidueTypeSet
	void
	default_rts(core::chemical::ResidueTypeSetCOP setting);

	/// @brief What is the default fall-back ResidueTypesSet
	core::chemical::ResidueTypeSetCOP
	default_rts() const;

	ResidueTypeCOP get_d_equivalent( ResidueTypeCOP l_rsd ) const override;

	ResidueTypeCOP get_l_equivalent( ResidueTypeCOP d_rsd ) const override;

	ResidueTypeCOP get_mirrored_type( ResidueTypeCOP original_rsd ) const override;

	/// @brief Check if a base type (like "SER") generates any types with another name3 (like "SEP")
	bool
	generates_patched_residue_type_with_name3( std::string const & base_residue_name, std::string const & name3 ) const override;

	/// @brief Check if a base type (like "CYS") generates any types with a new interchangeability group (like "SCY" (via cys_acetylated))
	bool
	generates_patched_residue_type_with_interchangeability_group( std::string const & base_residue_name, std::string const & interchangeability_group ) const override;

	/// @brief the residues with no patches applied
	ResidueTypeCOPs
	base_residue_types() const override;

	/// @brief the residues with no patches applied
	ResidueTypeCOPs
	unpatchable_residue_types() const override;

	/// @brief the patches
	utility::vector1< PatchCOP > patches() const override;

	/// @brief the metapatches
	utility::vector1< MetapatchCOP > metapatches() const override;

	/// @brief the patches, index by name.
	std::map< std::string, utility::vector1< PatchCOP > > const & patch_map() const override { return all_patch_map_; }

	/// @brief the metapatches, index by name.
	std::map< std::string, MetapatchCOP > const & metapatch_map() const override { return all_metapatch_map_; }

	/// @brief Do we have this metapatch?
	bool
	has_metapatch( std::string const & name ) const override;

	MetapatchCOP
	metapatch( std::string const & name ) const override;

	/// @brief Gets all types with the given aa type and variants
	/// @details The number of variants must match exactly. Variants can be custom variants.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	ResidueTypeCOPs
	get_all_types_with_variants_aa( AA aa, utility::vector1< std::string > const & variants ) const override;

	/// @brief Gets all types with the given aa type and variants, making exceptions for some variants.
	/// @details The number of variants must match exactly. Variants can be custom variants, but exceptions must
	///           be standard types, listed in VariantType.hh.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	ResidueTypeCOPs
	get_all_types_with_variants_aa( AA aa,
		utility::vector1< std::string > const & variants,
		utility::vector1< VariantType > const & exceptions ) const override;

	////////////////////////////////////////////////////////////////////////////////
	// Protected methods in the parent which we're hoisting to public for this class
public:

	// These could be hoisted by using statements, but PyRosetta currently has issues with that.
	void add_base_residue_type( ResidueTypeOP new_type ) override;
	void add_base_residue_type( std::string const &  filename ) override;
	void read_files_for_base_residue_types( utility::vector1< std::string > const & filenames ) override;
	void add_unpatchable_residue_type( ResidueTypeOP new_type ) override;
	void add_unpatchable_residue_type( std::string const &  filename ) override;
	void read_files_for_unpatchable_residue_types( utility::vector1< std::string > const & filenames ) override;
	void remove_base_residue_type( std::string const & name ) override;
	void remove_unpatchable_residue_type( std::string const & name ) override;
	void add_patches(
		utility::vector1< std::string > const & patch_filenames,
		utility::vector1< std::string > const & metapatch_filenames
	) override;
	void set_merge_behavior_manager( MergeBehaviorManagerCOP mbm) override;

	// Hoist the getters
	AtomTypeSetCOP atom_type_set() const { return ResidueTypeSet::atom_type_set(); }
	ElementSetCOP element_set() const { return ResidueTypeSet::element_set(); }
	MMAtomTypeSetCOP mm_atom_type_set() const { return ResidueTypeSet::mm_atom_type_set(); }
	orbitals::OrbitalTypeSetCOP orbital_type_set() const { return ResidueTypeSet::orbital_type_set(); }

	void atom_type_set(AtomTypeSetCOP atom_types);
	void element_set(ElementSetCOP elements);
	void mm_atom_type_set(MMAtomTypeSetCOP mm_atom_types);
	void orbital_type_set(orbitals::OrbitalTypeSetCOP orbital_types);

protected:

	/// @brief Add a patch object to the RTS
	void
	add_patch(PatchCOP p) override;

	/// @brief Add a metapatch object to the RTS
	void
	add_metapatch(MetapatchCOP p) override;

	/// @brief Generate the ResidueType from either/both of the local RTS and the default RTS.
	/// @details Assumes that we already have a write lock to the cache object for this RTS
	/// (but not for the default RTS).
	ResidueTypeCOP
	generate_residue_type_write_locked( std::string const & name_in ) const override;

	/// @brief Attempt to lazily load the given residue type from data.
	bool
	lazy_load_base_type( std::string const & rsd_base_name ) const override;


private:

	core::chemical::ResidueTypeSetCOP default_rts_;

	/// @brief patches indexed by name
	/// Includes both the patches from this map and the default_rts_ patches.
	std::map< std::string, utility::vector1< PatchCOP > > all_patch_map_;
	std::map< std::string, MetapatchCOP > all_metapatch_map_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //core
} //chemical


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_PoseResidueTypeSet )
#endif // SERIALIZATION


#endif //INCLUDED_core_chemical_PoseResidueTypeSet_hh





