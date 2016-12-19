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

	/// @brief Get ResidueType by exact name, returning COP
	/// Will return null pointer for no matches
	virtual
	ResidueTypeCOP
	name_mapOP( std::string const & name ) const override;

	/// @brief query if a ResidueType of the unique residue id (name) is present.
	virtual
	bool
	has_name( std::string const & name ) const override;

	/// @brief Check if a base type (like "SER") generates any types with another name3 (like "SEP")
	virtual
	bool
	generates_patched_residue_type_with_name3( std::string const & base_residue_name, std::string const & name3 ) const override;

	/// @brief Check if a base type (like "CYS") generates any types with a new interchangeability group (like "SCY" (via cys_acetylated))
	virtual
	bool
	generates_patched_residue_type_with_interchangeability_group( std::string const & base_residue_name, std::string const & interchangeability_group ) const override;

	/// @brief the residues with no patches applied
	virtual
	ResidueTypeCOPs
	base_residue_types() const override;

	/// @brief the residues with no patches applied
	virtual
	ResidueTypeCOPs
	unpatchable_residue_types() const override;

	/// @brief the patches
	virtual
	utility::vector1< PatchCOP > patches() const override;

	/// @brief the metapatches
	virtual
	utility::vector1< MetapatchCOP > metapatches() const override;

	/// @brief the patches, index by name.
	virtual
	std::map< std::string, utility::vector1< PatchCOP > > patch_map() const override;

	/// @brief Do we have this metapatch?
	virtual
	bool
	has_metapatch( std::string const & name ) const override;

	virtual
	MetapatchCOP
	metapatch( std::string const & name ) const override;

	/// @brief Gets all types with the given aa type and variants
	/// @details The number of variants must match exactly. Variants can be custom variants.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	virtual
	ResidueTypeCOPs
	get_all_types_with_variants_aa( AA aa, utility::vector1< std::string > const & variants ) const override;

	/// @brief Gets all types with the given aa type and variants, making exceptions for some variants.
	/// @details The number of variants must match exactly. Variants can be custom variants, but exceptions must
	///           be standard types, listed in VariantType.hh.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	virtual
	ResidueTypeCOPs
	get_all_types_with_variants_aa( AA aa,
		utility::vector1< std::string > const & variants,
		utility::vector1< VariantType > const & exceptions ) const override;

	////////////////////////////////////////////////////////////////////////////////
	// Protected methods in the parent which we're hoisting to public for this class
public:

	// These could be hoisted by using statements, but PyRosetta currently has issues with that.
	virtual void add_base_residue_type( ResidueTypeOP new_type ) override;
	virtual void add_base_residue_type( std::string const &  filename ) override;
	virtual void read_files_for_base_residue_types( utility::vector1< std::string > const & filenames ) override;
	virtual void add_unpatchable_residue_type( ResidueTypeOP new_type ) override;
	virtual void add_unpatchable_residue_type( std::string const &  filename ) override;
	virtual void read_files_for_unpatchable_residue_types( utility::vector1< std::string > const & filenames ) override;
	virtual void remove_base_residue_type( std::string const & name ) override;
	virtual void remove_unpatchable_residue_type( std::string const & name ) override;
	virtual void add_patches(
		utility::vector1< std::string > const & patch_filenames,
		utility::vector1< std::string > const & metapatch_filenames
	) override;
	virtual void set_merge_behavior_manager( MergeBehaviorManagerCOP mbm) override;

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

	/// @brief Attempt to lazily load the given residue type from data.
	bool
	lazy_load_base_type( std::string const & rsd_base_name ) const override;


private:

	core::chemical::ResidueTypeSetCOP default_rts_;

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





