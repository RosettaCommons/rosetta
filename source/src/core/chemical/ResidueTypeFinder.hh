// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/ResidueTypeFinder.hh
/// @brief Functions to find residue_type(s) from ResidueTypeSet without requiring instantiation of all rsd_types.
/// @detailed Intended to be super-efficient replacement for aa_map(), name3_map() in ResidueTypeSet.
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_chemical_ResidueTypeFinder_HH
#define INCLUDED_core_chemical_ResidueTypeFinder_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/chemical/ResidueTypeFinder.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/Patch.fwd.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace core {
namespace chemical {

class ResidueTypeFinder: public utility::pointer::ReferenceCount {

public:

	//constructor
	ResidueTypeFinder( core::chemical::ResidueTypeSet const & residue_type_set );

	//destructor
	~ResidueTypeFinder();

public:

	ResidueTypeCOP
	get_representative_type() const;

	ResidueTypeCOP
	get_representative_type_SLOW() const;

	ResidueTypeCOPs
	get_all_possible_residue_types( bool const allow_extra_variants = false ) const;

	ResidueTypeCOPs
	get_possible_base_residue_types() const;

	ResidueTypeCOP
	get_best_match_residue_type_for_atom_names( utility::vector1< std::string > const & atom_names );

	ResidueTypeFinder &
	aa( core::chemical::AA const & setting ) {
		aa_ = setting;
		return *this;
	}

	ResidueTypeFinder &
	name1( char const & setting ) {
		name1_ = setting;
		return *this;
	}

	ResidueTypeFinder &
	name3( std::string const & setting ) {
		name3_ = setting;
		return *this;
	}

	ResidueTypeFinder &
	residue_base_name( std::string const & setting ) {
		residue_type_base_name_ = setting;
		return *this;
	}

	ResidueTypeFinder &
	base_property( ResidueProperty const setting ) {
		base_property_ = setting;
		return *this;
	}

	ResidueTypeFinder &
	atom_names_soft( utility::vector1< std::string > const & setting ) {
		atom_names_soft_ = setting;
		return *this;
	}

	ResidueTypeFinder &
	variants_in_sets( utility::vector1< utility::vector1< VariantType > > const & setting ) {
		variants_in_sets_ = setting;
		return *this;
	}

	ResidueTypeFinder &
	variants( utility::vector1< VariantType > const & setting );

	ResidueTypeFinder &
	variants( utility::vector1< std::string > const & setting );

	ResidueTypeFinder &
	disallow_variants( utility::vector1< VariantType > const & setting ) {
		disallow_variants_ = setting;
		return *this;
	}

	ResidueTypeFinder &
	properties( utility::vector1< ResidueProperty > const & setting ) {
		properties_ = setting;
		return *this;
	}

	ResidueTypeFinder &
	disallow_properties( utility::vector1< ResidueProperty > const & setting ) {
		disallow_properties_ = setting;
		return *this;
	}

	ResidueTypeFinder &
	patch_names( utility::vector1< std::string > const & setting ) {
		patch_names_ = setting;
		return *this;
	}

	ResidueTypeFinder &
	ignore_atom_named_H( bool const setting ) {
		ignore_atom_named_H_ = setting;
		return *this;
	}

	ResidueTypeFinder &
	match_adducts( bool const setting ) {
		match_adducts_ = setting;
		return *this;
	}

	ResidueTypeFinder &
	disallow_carboxyl_conjugation_at_glu_asp( bool const setting ) {
		disallow_carboxyl_conjugation_at_glu_asp_ = setting;
		return *this;
	}

private:

	ResidueTypeCOPs
	apply_filters_after_patches( ResidueTypeCOPs rsd_types,
		bool const allow_extra_variants = false ) const;

	utility::vector1< ResidueTypeCOP >
	apply_patches_recursively( utility::vector1< ResidueTypeCOP > const & rsd_types,
		Size const patch_number,
		bool const get_first_totally_ok_residue_type = false ) const;


	ResidueTypeCOPs
	filter_by_name1( ResidueTypeCOPs const & rsd_types  ) const;

	ResidueTypeCOPs
	filter_by_name3( ResidueTypeCOPs const & rsd_types,
		bool const keep_if_base_type_generates_name3) const;

	ResidueTypeCOPs
	filter_by_aa( ResidueTypeCOPs const & rsd_types  ) const;

	ResidueTypeCOPs
	filter_by_residue_type_base_name( ResidueTypeCOPs const & rsd_types ) const;

	ResidueTypeCOPs
	filter_by_base_property( ResidueTypeCOPs const & base_rsd_types ) const;

	ResidueTypeCOPs
	check_candidate_has_all_variant_sets( ResidueTypeCOPs const & rsd_types ) const;

	ResidueTypeCOPs
	check_variant_sets_have_all_candidate_variants( ResidueTypeCOPs const & rsd_types ) const;

	ResidueTypeCOPs
	filter_all_variants_matched( ResidueTypeCOPs const & rsd_types,
		bool const allow_extra_variants = false ) const;

	ResidueTypeCOPs
	filter_disallow_variants( ResidueTypeCOPs const & rsd_types )  const;

	ResidueTypeCOPs
	filter_all_patch_names( ResidueTypeCOPs const & rsd_types )  const;

	ResidueTypeCOPs
	filter_all_properties( ResidueTypeCOPs const & rsd_types )  const;

	ResidueTypeCOPs
	filter_disallow_properties( ResidueTypeCOPs const & rsd_types )  const;

	ResidueTypeCOPs
	filter_special_cases( ResidueTypeCOPs const & rsd_types )  const;

	bool
	has_disallowed_variant( PatchCOP patch ) const;

	bool
	adds_any_variant( PatchCOP patch ) const;

	bool
	fixes_name3( PatchCOP patch,
		ResidueTypeCOP rsd_type ) const;

	bool
	adds_any_property( PatchCOP patch,
		ResidueTypeCOP rsd_type ) const;

	bool
	deletes_any_property( PatchCOP patch,
		ResidueTypeCOP rsd_type ) const;

	bool
	deletes_any_variant( PatchCOP patch,
		ResidueTypeCOP rsd_type ) const;

	bool
	matches_any_patch_name( PatchCOP patch ) const;

	bool
	matches_any_atom_name( PatchCOP patch,
		ResidueTypeCOP rsd_type ) const;

	ResidueTypeFinder &
	variants( utility::vector1< VariantType > const & setting ) const;


private:

	core::chemical::ResidueTypeSet const & residue_type_set_;

	core::chemical::AA aa_;
	char name1_;
	std::string name3_;
	std::string residue_type_base_name_;
	utility::vector1< std::string > atom_names_soft_;
	utility::vector1< utility::vector1< VariantType > > variants_in_sets_;
	utility::vector1< std::string > custom_variants_;
	utility::vector1< VariantType > disallow_variants_;
	utility::vector1< ResidueProperty > properties_;
	utility::vector1< ResidueProperty > disallow_properties_;
	utility::vector1< std::string > patch_names_;
	ResidueProperty base_property_;
	bool ignore_atom_named_H_;
	bool match_adducts_; // for DNA -- actually not in use, but keep in here if adducts are modernized.
	bool disallow_carboxyl_conjugation_at_glu_asp_; // special case

};

} //chemical
} //core

#endif
