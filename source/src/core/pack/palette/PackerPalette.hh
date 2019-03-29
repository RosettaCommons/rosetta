// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/PackerPalette.hh
/// @brief  PackerPalette: a class for storing the set of ResidueTypes
/// and VariantTypes that the packer uses by default, in the absence of any
/// TaskOperations that limit the set actually used.
/// @details The PackerPalette says, "Here are the types that you're
/// allowed to use, and which are on in the absence of TaskOperations."
/// TaskOperations then prune this, turning OFF types that have been
/// enabled.  This allows users to turn on noncanonicals for design, and
/// then use TaskOperations with the same commutativity rules (turning OFF
/// types only) that are used for canonicals, making mixed design with
/// canonicals and noncanonicals much easier.\nThis was implemented as
/// part of the 2016 Chemical XRW (eXtreme Rosetta Workshop).
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu).

#ifndef INCLUDED_core_pack_palette_PackerPalette_hh
#define INCLUDED_core_pack_palette_PackerPalette_hh

// Unit Headers
#include <core/pack/palette/PackerPalette.fwd.hh>

// Package Headers

// Project Headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// STL Headers
#include <map>
#include <list>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace palette {

enum SpecialPackerPaletteBehaviour {
	KEEP_EXISTING_BASE_TYPE = 1, //At a given position, keep the existing base type in the allowed list.  (Keep this first.)
	FORCE_EXISTING_BASE_TYPE, //At a given position, allow NOTHING but the exsiting base type (i.e. no design).  False by default.
	ONLY_DESIGN_POLYMER_RESIDUES, //If a given positions is NOT a polymer residue type, then (a) keep the existing type, and (b) don't add any additional types.
	ONLY_DESIGN_PROTEIN_PEPTIOID_DNA_AND_SACCHARIDES, //The old default behaviour.  Non-protein, non-peptoid, non-DNA, non-saccharide positions are limited to repacking.
	ONLY_ALPHA_AA_AT_ALPHA_POSITIONS, //If a given position is an alpha-amino acid, only allow alpha-amino acids at that position.
	ONLY_BETA_AA_AT_BETA_POSITIONS, //If a given position is a beta-amino acid, only allow beta-amino acids at that position.
	ONLY_GAMMA_AA_AT_GAMMA_POSITIONS, //If a given position is a gamma-amino acid, only allow gamma-amino acids at that position.
	ONLY_DNA_AT_DNA_POSITIONS, //If a given position is DNA, only allow DNA types at that position.
	ONLY_OLIGOUREA_AT_OLIGOUREA_POSITIONS, //If a given position is an oligourea, only allow oligourea types at that position.
	ONLY_ARAMID_AT_ARAMID_POSITIONS, //If a given position is an aramid, only allow compatible aramid types at that position.
	ONLY_SACCHARIDES_AT_SACCHARIDE_POSITIONS, // If a given position is a sugar, only allow sugar types at that position.
	ONLY_LIGAND_AT_LIGAND_POSITIONS, //If a given position is a ligand, only allow ligand types at that position.
	ONLY_MATCHING_LIGAND_NAMES, //If a position is a ligand, only allow other ligands at that position with the same name.
	ALL_DNA_TYPES_ON, //Irritatingly, the old default was for ALL DNA types to be on by default when designing -- not just canoncal types.
	ONLY_RNA_AT_RNA_POSITIONS, //If a given position is RNA, only allow RNA types at that position.
	EXCLUDE_ADDUCT_VARIANT_AT_DNA_POSITIONS, //If a given position is DNA, when designing it, don't match on the ADDUCT_VARIANT type.  This is a special behaviour that has to be carried over from ResidueLevelTask_ to preserve old functionality.
	STRIP_VIRTUAL_SIDE_CHAIN, //If a given position has a virtual side chain, we convert it automatically to the same without a virtual side chain.
	pH_MODE_EXCEPTIONS, //Are we using exceptions for pH mode, so that PROTONATED/DEPROTONATED variants can be matched?
	KEEP_EXISTING_TERMINAL_VARIANT_TYPES_AT_POSITIONS, //If a position has a modified terminus, only allow types at that position with that modified terminus.
	KEEP_EXISTING_NONTERMINAL_VARIANT_TYPES_FOR_EXISTING_BASE_TYPE, // At a given position, keep the existing variant types in the allowed list (e.g. if phosphorylation is present on a serine, keep phosophoserine as an allowed type at that position).
	KEEP_EXISTING_NONTERMINAL_VARIANT_TYPES_AND_DISALLLOW_INCOMPATIBLE_BASE_TYPES, //At a given position, keep the existing variant types.  Disallow base types that don't have that variant type.
	KEEP_EXISTING_DISULFIDES, // If a residue is involved in a disulfide, only allow the disulfide type at that position.
	NO_METAPATCHES, // Unless the current position is metapatched, do not consider metapatched ResidueType variants.  Needed for speed.
	ALLOW_ALTERNATE_BACKBONE_MATCHING, // Rosetta will attempt to design residues of different backbone families, if their connections can superimpose.
	INVALID_BEHAVIOUR, //Keep this second-to-last.
	END_OF_BEHAVIOUR_LIST = INVALID_BEHAVIOUR //Keep this last.
};


/// @brief  The PackerPalette class gives instructions to the packer about
/// the set of ResidueTypes and VariantTypes to use by default, in the
/// absence of any TaskOperations that prune the list.
class PackerPalette : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< PackerPalette >
{

public:
	/// @brief Default constructor.
	/// @details Defaults to the fa_standard ResidueTypeSet.
	PackerPalette();

	/// @brief Constructor with ResidueTypeSet specifier.
	/// @details Needed for packing with anything but the default (fa_standard) ResidueTypeSet.
	PackerPalette( core::chemical::ResidueTypeSetCOP residue_type_set_in );

	/// @brief Copy constructor.
	PackerPalette( PackerPalette const &src );

	/// @brief Destructor.
	~PackerPalette();

	/// @brief Clone operator: make a copy and return an owning pointer to the copy.
	/// @details Derived classes MUST implement this.
	virtual PackerPaletteOP clone() const = 0;

	/// @brief Self pointers (const).
	inline PackerPaletteCOP get_self_ptr() const { return shared_from_this(); }

	/// @brief Self pointers (non-const).
	inline PackerPaletteOP get_self_ptr() { return shared_from_this(); }

	/// @brief The initialize_residue_level_task ("apply") function -- called during rotamer setup for the packer to get the list of all ResidueTypes that are allowed.
	/// @details Derived classes don't get to implement this to override the default behaviour (which is just to set up all allowed ResidueTypes that result from
	/// the combination of the base types listed and the VariantTypes listed), but specific variations on the default behaviour, provided by enum, are permitted.
	/// Note that the PackerPalette is not meant to be used for position-specific setup, despite having access to the pose residue.  Use TaskOperations for that.
	/// @param [in] existing_residue The existing residue, for reference (though this should be largely unneeded).  It's largely only used for variant type matching.
	/// @param [out] residue_type_list A std::list of ResidueTypeCOPs, which is cleared and populated by this function.
	void //NOT virtual.
	initialize_residue_level_task(
		core::conformation::Residue const & existing_residue,
		std::list< chemical::ResidueTypeCOP > & residue_type_list
	) const;

	/// @brief Function to parse XML tags, implemented by derived classes.
	/// @brief Failure to implement this results in a warning message, but does not prevent compilation or
	/// program execution.
	virtual void
	parse_my_tag(
		utility::tag::TagCOP const &tag,
		basic::datacache::DataMap const &datamap
	);

	/// @brief Get the name of this class.
	/// @details Must be implemented by derived classes.
	virtual std::string const & name() const = 0;

	/// @brief Function to allow a different ResidueTypeSet to be set.
	/// @details Each PackerPalette derived class must implement this.  After setting the new ResidueTypeSet, things need to happen.
	virtual void set_residue_type_set( core::chemical::ResidueTypeSetCOP new_type_set );

	/// @brief Add a base residue type name to the list of base residue type names.
	void add_base_residue_type( std::string const &name );

	/// @brief Add a group of base ResidueTypes and names to the PackerPalette by properties.
	void add_base_residue_types_by_properties( utility::vector1< core::chemical::ResidueProperty > const & properties );

	/// @brief Does the PackerPalette have a base residue type?
	bool has_base_residue_type( std::string const &name );

	/// @brief Clear the base residue type list.
	void clear_base_residue_type_list();

	/// @brief Set up the default base types:
	void set_up_default_base_types();
	void set_up_expanded_base_types();

	/// @brief Set the special behaviours to their default values.
	void set_up_default_special_behaviours();

	/// @brief Get a const owning pointer to the ResidueTypeSet.
	inline core::chemical::ResidueTypeSetCOP residue_type_setCOP() const { return residue_type_set_; }


protected: //Protected functions:

	/// @brief Set whether design is only limited to protein/peptoid/dna/saccharide positions, or can happen
	/// everywhere.
	/// @details The DefaultPackerPalette has this set to true; everything else should set it to false.
	void set_only_design_protein_peptoid_dna_saccharide( bool const setting );

	/// @brief Set whether we're forcing the existing base type.
	/// @details Defaults to false.  The NoDesignPackerPalette sets this to true.
	void set_force_existing_base_type( bool const setting );

	/// @brief Set whether design is only limited to polymer positions, or can happen everywhere.
	/// @details The DefaultPackerPalette has this set to true; everything else should set it to false.
	void set_only_design_polymer_residues( bool const setting );

	/// @brief Return a list of base ResidueType names to be added to the PackerPalette by properties.
	utility::vector1< std::string > get_base_names_by_properties(
		utility::vector1< core::chemical::ResidueProperty > const & properties );


private: //Private class methods:

	/// @brief Given a residue type pointer, get a raw pointer to the base type.
	/// @details Avoids calls to get_self_ptr().
	core::chemical::ResidueType const * get_base_type_raw_ptr( core::chemical::ResidueTypeCOP const & restype ) const;

	/// @brief Get a residue type with a variant removed.
	core::chemical::ResidueTypeCOP get_residue_type_cop_with_variant_removed( core::chemical::ResidueTypeCOP const & restype, core::chemical::VariantType const vartype_to_remove ) const;

	/// @brief Set the special behaviours to be a map of each behaviour type to FALSE.
	/// @details All behaviours default to false, and must be explicitly set to true.
	void initialize_special_behaviours();

	/// @brief Given the existing residue, decide what to do with its base type.
	/// @details Depending on what's defined in the special_behaviours_ list, we might keep the base type,
	/// keep the base type with terminal variants, keep the base type with terminal and side-chain variants,
	/// etc.
	void decide_what_to_do_with_existing_type(
		core::conformation::Residue const & existing_residue,
		std::list< core::chemical::ResidueTypeCOP> & residue_type_list
	) const;

	/// @brief Given the existing residue and a candidate base type in the base types list, decide what to
	/// do with the candidate base type.
	void decide_what_to_do_with_base_type(
		core::conformation::Residue const & existing_residue,
		std::list< core::chemical::ResidueTypeCOP> & residue_type_list,
		core::chemical::ResidueTypeCOP const & candidate_base_type,
		bool &existing_type_processed
	) const;

	/// @brief Set up the list of variant types that don't modify termini.
	void initialize_non_terminal_types();

	/// @brief Given a list of VariantTypes and a residue, populate two new lists with the VariantTypes present on the
	/// residue and the custom VariantTypes present.
	void get_types_on_residue(
		core::conformation::Residue const & residue,
		utility::vector1 < core::chemical::VariantType > const & types,
		utility::vector1 < core::chemical::VariantType > & present_types,
		utility::vector1 < std::string > & present_on_the_fly_types
	) const;

	/// @brief Has the user specified that pre-talaris behaviour should be restored?
	/// @details Having to support this is a royal pain in the neck.
	inline bool restore_pre_talaris_behaviour() const { return restore_pre_talaris_behaviour_; }

	/// @brief has the icoor_05_2009 option been used?
	/// @details More backwarkwards-compatibility malarkey.
	inline bool icoor_05_2009() const { return icoor_05_2009_; }

	/// @brief Has the pH_mode option been used?
	/// @details Still more backwarkwards-compatibility malarkey.
	inline bool pH_mode() const { return pH_mode_; }

	/// @brief Given some ResidueTypes to add to the ResidueType list, add the ones that are not already in the list.
	/// @details Determines equivalency by pointer address comparison.
	void add_residue_types_to_list( utility::vector1< core::chemical::ResidueTypeCOP > const &types_to_add, std::list< core::chemical::ResidueTypeCOP > &residue_type_list, core::Size const seqpos ) const;

	/// @brief Do these two ResidueTypes have compatible backbones?
	/// @author Jason W. Labonte
	bool saccharide_backbones_are_compatible(
		chemical::ResidueType const & existing_type,
		chemical::ResidueType const & candidate_type ) const;

	/// @brief Do these two ResidueTypes have compatible aramid backbones?
	/// @note Returns false if either residue is not an aramid, or if they are both aramids
	/// but of incompatible types.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	bool aramid_backbones_are_compatible(
		core::chemical::ResidueType const &existing_type,
		core::chemical::ResidueType const &candidate_type
	) const;

private: //Private member variables:

	/// @brief Has the user specified that pre-talaris behaviour should be restored?
	/// @details Having to support this is a royal pain in the neck.
	bool restore_pre_talaris_behaviour_;

	/// @brief Has the icoor_05_2009 option been used?
	/// @details More backwarkwards-compatibility malarkey.
	bool icoor_05_2009_;

	/// @brief Has the pH_mode option been used?
	/// @details Still more backwarkwards-compatibility malarkey.
	bool pH_mode_;

	/// @brief List base residue type names and ResidueTypeCOPs.
	/// @details Stores the string and a const owning pointer to the type.
	utility::vector1 < std::pair< std::string, core::chemical::ResidueTypeCOP > > base_residue_type_names_;

	/// @brief A const-access owning pointer to the ResidueTypeSet.
	/// @author Andy Watkins.
	core::chemical::ResidueTypeSetCOP residue_type_set_;

	/// @brief List of variant type names.
	utility::vector1 < std::string > variant_type_names_;

	/// @brief Settings for special behaviours.
	/// @details Uses the SpecialPackerPaletteBehaviour enum, defined in PackerPalette.hh.
	/// Certain enumerated position-specific behaviours are allowed (e.g. preserving special variant types,
	/// like phosphorylation, at a position), but generally,
	/// TaskOperations should be used for position-specific behaviour, not PackerPalettes.
	std::map < SpecialPackerPaletteBehaviour, bool > special_behaviours_;

	/// @brief List of all VariantTypes that exist that modify termini.
	/// @details Loaded from core::chemical::get_terminal_types() on PackerPalette creation.
	utility::vector1 < core::chemical::VariantType > terminal_types_;

	/// @brief List of all VariantTypes that exist that do not modify termini.
	/// @details Generated from terminal_types_ on PackerPalette creation.
	utility::vector1 < core::chemical::VariantType > non_terminal_types_;

#ifdef SERIALIZATION
public:

	/// @brief "Save" function for serialization.
	///
	template< class Archive > void save( Archive & arc ) const;

	/// @brief "Load" function for serialization.
	///
	template< class Archive > void load( Archive & arc );

	/// @brief Given a utility::vector1 < std::pair< std::string, core::chemical::ResidueTypeCOP > >,
	/// serialize it.
	template < class Archive >
	void serialize_base_residue_type_names(
		Archive &archive,
		utility::vector1 < std::pair< std::string, core::chemical::ResidueTypeCOP > > const &vect
	) const;

	/// @brief Given a utility::vector1 < std::pair< std::string, core::chemical::ResidueTypeCOP > >,
	/// deserialize it.
	template < class Archive >
	void deserialize_base_residue_type_names(
		Archive &archive,
		utility::vector1 < std::pair< std::string, core::chemical::ResidueTypeCOP > > &vect
	) const;

	/// @brief Given a vector of VariantTypes, serialize it.
	///
	template < class Archive >
	void serialize_VariantType_vector(
		Archive &archive,
		utility::vector1 < core::chemical::VariantType > const &vect
	) const;

	/// @brief Given a vector of VariantTypes, deserialize it.
	///
	template < class Archive >
	void deserialize_VariantType_vector(
		Archive &archive,
		utility::vector1 < core::chemical::VariantType > &vect
	) const;

	/// @brief Given a std::map of SpecialPackerPaletteBahaviour enum to bool,
	/// serialize it.
	template < class Archive >
	void serialize_behaviours_map(
		Archive &archive,
		std::map < SpecialPackerPaletteBehaviour, bool > const &mymap
	) const;

	/// @brief Given a std::map of SpecialPackerPaletteBahaviour enum to bool,
	/// deserialize it.
	template < class Archive >
	void deserialize_behaviours_map(
		Archive &archive,
		std::map < SpecialPackerPaletteBehaviour, bool > &mymap
	) const;


#endif // SERIALIZATION

};


} //namespace palette
} //namespace pack
} //namespace core

#endif
