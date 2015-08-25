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

#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/util.hh>

#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/string.functions.hh>
#include <utility/stream_util.hh>

using utility::vector1;
using utility::tools::make_vector1;

static basic::Tracer TR( "core.chemical.ResidueTypeFinder" );

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Intended to be a super-efficient way to get ResidueType's from a ResidueTypeSet, without
//  having to instantiate *everything* in the ResidueTypeSet, which can take an exponentially
//  long time if there are a lot of patches.
//
// Notes:
//  * caches any new ResidueTypes in ResidueTypeSet.
//
//  * does a binary tree search adding patches, one-by-one.
//
//  * for get_representative_type(), best to use as few desired constraints as possible.
//
//  * for get_all_residue_types() or get_best_match_residue_type_for_atom_names(), best to use
//     as many constraints as possible that are specific to the desired residue type(s).
//
// TO DO:
//  * probably should refactor all the checks/filters into functionalities for
//     chemical::ResidueTypeSelector.  Something like:
//             bool check_patch_worth_trying( Patch, ResidueType )
//
//  * need to properly handle replace_residue_type() functionality (those should be required patches...)
//
//  * when grabbing base_residue_types, may want to be smarter about pulling out that list --
//       in the future may want to have 100,000 ligands in Rosetta, which could be kept as bare-bones
//       'placeholder' ResidueTypes with just name3. These would be instantiated only when recognized
//       in say a PDB file or requested explicitly by name3.
//
// -- rhiju, June 2015
//
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace chemical {

//Constructor
ResidueTypeFinder::ResidueTypeFinder( core::chemical::ResidueTypeSet const & residue_type_set ):
	residue_type_set_( residue_type_set ),
	aa_( aa_none ),
	name1_( '?' ),
	base_property_( NO_PROPERTY ),
	ignore_atom_named_H_( false ),
	match_adducts_( true ),
	disallow_carboxyl_conjugation_at_glu_asp_( false )
{}

//Destructor
ResidueTypeFinder::~ResidueTypeFinder()
{}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOP
ResidueTypeFinder::get_representative_type() const
{
	if ( custom_variants_.size() > 0 ) return get_representative_type_SLOW();
	ResidueTypeCOPs rsd_types = get_possible_base_residue_types();
	rsd_types = apply_patches_recursively( rsd_types, 1 /*start with this patch*/, true /*get_first_residue_found*/ );
	rsd_types = apply_filters_after_patches( rsd_types, true /* allow_extra_variants */ );
	if ( rsd_types.size() == 0 ) return 0;
	return rsd_types[ 1 ];
}


ResidueTypeCOP
ResidueTypeFinder::get_representative_type_SLOW() const
{
	// this will instantiate *everything*. Super-slow. The way things were done with the OLD ResidueTypeSet.
	ResidueTypeCOPs rsd_types( residue_type_set_.residue_types_DO_NOT_USE() );

	// this could be made more efficient by exiting after the first rsd_type is found.
	rsd_types = filter_by_residue_type_base_name( rsd_types );
	rsd_types = filter_all_variants_matched( rsd_types );
	if ( rsd_types.size() == 0 ) return 0;
	return rsd_types[ 1 ];
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::get_all_possible_residue_types( bool const allow_extra_variants /* = false */ ) const
{
	// Get all possible basic residues that might match
	ResidueTypeCOPs rsd_types = get_possible_base_residue_types();

	// Go down the binary tree of patches.
	rsd_types = apply_patches_recursively( rsd_types, 1 /*start with this patch*/ );

	// Filter for rsd_types that strictly obey requirements
	rsd_types = apply_filters_after_patches( rsd_types, allow_extra_variants );

	return rsd_types;

}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOP
ResidueTypeFinder::get_best_match_residue_type_for_atom_names( utility::vector1< std::string > const & atom_names )
{

	// clock_t const time_start( clock() );

	// will try to match these ('soft' constraints). Go ahead and strip out whitespace.
	atom_names_soft_.clear();
	for ( Size n = 1; n <= atom_names.size(); n++ ) {
		std::string atom_name_temp = atom_names[ n ];
		atom_names_soft_.push_back( ObjexxFCL::strip_whitespace( atom_name_temp ) );
	}

	ResidueTypeCOPs rsd_types = get_all_possible_residue_types( true /* allow_extra_variants */ );
	ResidueTypeCOP  rsd_type  = find_best_match( rsd_types, atom_names, ignore_atom_named_H_ );

	//  TR << "time to initialize " << rsd_type->name() << " from " << rsd_types.size() << " possible ResidueTypes: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;

	return rsd_type;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::get_possible_base_residue_types() const
{
	ResidueTypeCOPs rsd_types = residue_type_set_.base_residue_types();
	rsd_types = filter_by_aa( rsd_types );
	rsd_types = filter_by_name1( rsd_types );
	rsd_types = filter_by_name3( rsd_types, true /* keep_if_base_type_generates_name3 */ );
	rsd_types = filter_by_residue_type_base_name( rsd_types );
	rsd_types = filter_by_base_property( rsd_types );
	return rsd_types;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::apply_filters_after_patches( ResidueTypeCOPs rsd_types,
	bool const allow_extra_variants  /* = false */ ) const
{
	rsd_types = filter_by_name3( rsd_types, false /* keep_if_base_type_generates_name3 */ );
	rsd_types = filter_disallow_variants( rsd_types );
	rsd_types = filter_all_variants_matched( rsd_types, allow_extra_variants );
	rsd_types = filter_all_properties(  rsd_types );
	rsd_types = filter_disallow_properties(  rsd_types );
	rsd_types = filter_all_patch_names( rsd_types );
	rsd_types = filter_special_cases( rsd_types );
	return rsd_types;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Instantiates ResidueType (and gets exponentially slower), based on the number
//  of desired features that a patch might offer:
//
// (1) xyz_atom_names   Does patch introduce this atom via ADD_ATOM?
// (2) variants         Does patch introduce this VariantType through its TYPES?
// (3) properties       Does patch introduce this ResidueProperty through ADD_PROPERTY?
// (4) patch_names      Used for carbohydrate branching --> look explicitly for name.
//
// This function will take less time and memory if the lists above are *short*.
//
// So you can be smart about tuning performance -- for example, if you are looking for "->2)-branch"
//  carbohydrates, set that through patch_names; and do not also look for BRANCH_POINT under properties (which
//  is also a property of patches applied to amino acids).
//
//                 -- rhiju, 2015
//
////////////////////////////////////////////////////////////////////////////////////////////////////
vector1< ResidueTypeCOP >
ResidueTypeFinder::apply_patches_recursively(
	vector1< ResidueTypeCOP > const & rsd_types,
	Size const patch_number,
	bool const get_first_totally_ok_residue_type /*= false*/
) const {

	ResidueTypeCOPs rsd_types_new = rsd_types;
	PatchCOP patch = residue_type_set_.patches()[ patch_number ];

	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP rsd_type = rsd_types[ n ];

		// absolute no-no's.
		if ( !patch->applies_to( *rsd_type ) )             continue;
		if ( has_disallowed_variant( patch ) )             continue;
		if ( deletes_any_property(   patch, rsd_type ) )   continue;
		if ( deletes_any_variant(    patch, rsd_type ) )   continue;
		// could also add as a no-no: if patch *deletes* an atom in atom_names_.

		// note -- make sure to apply patch if it has a chance of satisfying any of
		// the constraints on variants, branchpoints, or properties.
		bool apply_patch = (  adds_any_variant ( patch ) ||
			adds_any_property( patch, rsd_type ) ||
			matches_any_patch_name( patch )      ||
			matches_any_atom_name( patch, rsd_type ) ||
			fixes_name3( patch, rsd_type ) );

		if ( apply_patch ) {
			//std::cout << "amw patch name was " << patch->name();
			ResidueTypeCOP rsd_type_new_placeholder = patch->apply( *rsd_type, false /*instantiate*/ );
			//std::cout << "amw placeholder RT name was " << rsd_type_new_placeholder->name();
			ResidueTypeCOP rsd_type_new( rsd_type->residue_type_set().name_map( rsd_type_new_placeholder->name() ).get_self_ptr() );
			//std::cout << "amw used name map to get " << rsd_type_new->name();

			rsd_types_new.push_back( rsd_type_new );
		}

	} // end loop

	if ( get_first_totally_ok_residue_type ) { // maybe we're done?
		// note that this repeats some work -- some rsd_types were checked in prior steps in the recursion
		ResidueTypeCOPs rsd_types_filtered = apply_filters_after_patches( rsd_types_new, true /*allow_extra_variants*/ );
		if ( rsd_types_filtered.size() > 0 ) return rsd_types_filtered;
	}

	// end of recursion through patches?
	if ( patch_number == residue_type_set_.patches().size() ) return rsd_types_new;

	return  apply_patches_recursively( rsd_types_new, patch_number + 1, get_first_totally_ok_residue_type );

}


////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_by_name1( ResidueTypeCOPs const & rsd_types  ) const
{
	if ( name1_ == '?' ) return rsd_types;
	ResidueTypeCOPs rsd_types_new;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP rsd_type( rsd_types[ n ] );
		if ( rsd_type->name1() != name1_ ) continue;
		rsd_types_new.push_back( rsd_type );
	}
	return rsd_types_new;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_by_name3( ResidueTypeCOPs const & rsd_types, bool const keep_if_base_type_generates_name3  ) const
{
	if ( name3_.size() == 0 ) return rsd_types;
	ResidueTypeCOPs rsd_types_new;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP rsd_type( rsd_types[ n ] );

		// For specialty amino acids, add them to the name three maps both with their PDB strings and
		// with their specialty string -- the first three letters of the residue name.
		// E.g., CYD will appear in both lists for name3_map_[ "CYS" ] and name3_map_[ "CYD" ]
		if ( rsd_type->name3() == name3_ ||
				rsd_type->name().substr(0,3) == name3_ ||
				( keep_if_base_type_generates_name3 && residue_type_set_.generates_patched_residue_type_with_name3( rsd_type->name(), name3_ ) ) ) {
			rsd_types_new.push_back( rsd_type );
		}
	}
	return rsd_types_new;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_by_aa( ResidueTypeCOPs const & rsd_types  ) const
{
	if ( aa_ == aa_none ) return rsd_types;

	ResidueTypeCOPs rsd_types_new;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP rsd_type( rsd_types[ n ] );
		if ( rsd_type->aa() != aa_ ) continue;
		rsd_types_new.push_back( rsd_type );
	}
	return rsd_types_new;
}

////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_by_residue_type_base_name( ResidueTypeCOPs const & rsd_types ) const
{
	if ( residue_type_base_name_.size() == 0 ) return rsd_types;

	ResidueTypeCOPs filtered_rsd_types;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP const & rsd_type = rsd_types[ n ];
		if ( residue_type_base_name( *rsd_type ) != residue_type_base_name_ ) continue;
		filtered_rsd_types.push_back( rsd_type );
	}
	return filtered_rsd_types;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_by_base_property( ResidueTypeCOPs const & rsd_types ) const
{
	if ( base_property_ == NO_PROPERTY ) return rsd_types;
	ResidueTypeCOPs filtered_rsd_types;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP rsd_type = rsd_types[ n ];
		if ( rsd_type->properties().has_property( base_property_ ) ) {
			filtered_rsd_types.push_back( rsd_type );
		}
	}
	return filtered_rsd_types;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResidueTypeFinder::has_disallowed_variant( PatchCOP patch ) const
{
	vector1< std::string> const & patch_variant_types = patch->types();
	for ( Size k = 1; k <= patch_variant_types.size(); k++ )  {
		VariantType const patch_variant = ResidueProperties::get_variant_from_string( patch_variant_types[ k ] );
		if ( disallow_variants_.has_value( patch_variant ) ) return true;
	}
	return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResidueTypeFinder::adds_any_variant( PatchCOP patch ) const
{
	vector1< std::string> const & patch_variant_types = patch->types();
	for ( Size n = 1; n <= patch_variant_types.size(); n++ ) {
		VariantType const patch_variant = ResidueProperties::get_variant_from_string( patch_variant_types[ n ] );

		for ( Size k = 1; k <= variants_in_sets_.size(); k++ ) {
			if ( variants_in_sets_[ k ].has_value( patch_variant ) ) return true;
			if ( !match_adducts_ && patch_variant == ADDUCT_VARIANT ) return true; // explore all adducts
		}

	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResidueTypeFinder::fixes_name3( PatchCOP patch, ResidueTypeCOP rsd_type ) const
{
	if ( name3_.size() > 0 &&
			rsd_type->name3() != name3_ &&
			patch->apply( *rsd_type, false /*instantiate*/ )->name3() == name3_ ) {
		return true;
	}
	return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResidueTypeFinder::adds_any_property( PatchCOP patch, ResidueTypeCOP rsd_type ) const
{
	vector1< std::string> const & added_properties( patch->adds_properties( *rsd_type ) );
	for ( Size n = 1; n <= added_properties.size(); n++ ) {
		// convert string to ResidueProperty
		ResidueProperty const added_property( rsd_type->properties().get_property_from_string( added_properties[ n ] ) );
		if ( properties_.has_value( added_property ) ) return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResidueTypeFinder::deletes_any_property( PatchCOP patch, ResidueTypeCOP rsd_type ) const
{
	vector1< std::string> const & deleted_properties( patch->deletes_properties( *rsd_type ) );
	for ( Size n = 1; n <= deleted_properties.size(); n++ ) {
		// convert string to ResidueProperty
		ResidueProperty const deleted_property( rsd_type->properties().get_property_from_string( deleted_properties[ n ] ) );
		if ( properties_.has_value( deleted_property ) ) return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResidueTypeFinder::deletes_any_variant( PatchCOP patch, ResidueTypeCOP rsd_type ) const
{
	vector1< std::string > const & deleted_variants( patch->deletes_variants( *rsd_type ) );
	for ( Size n = 1; n <= deleted_variants.size(); n++ ) {
		// convert string to ResidueProperty
		VariantType const deleted_variant( ResidueProperties::get_variant_from_string( deleted_variants[ n ] ) );
		for ( Size k = 1; k <= variants_in_sets_.size(); k++ ) {
			vector1< VariantType > const & variant_set = variants_in_sets_[ k ];
			for ( Size m = 1; m <= variant_set.size(); m++ ) {
				if ( variant_set.has_value( deleted_variant ) ) return true;
			} // loop inside one variant set
		}
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResidueTypeFinder::matches_any_patch_name( PatchCOP patch ) const
{
	return ( patch_names_.has_value( patch->name() ) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// check for any matching atom name -- that's a good sign that we should check out this residue_type.
bool
ResidueTypeFinder::matches_any_atom_name( PatchCOP patch, ResidueTypeCOP rsd_type ) const
{
	vector1< std::string > new_atom_names = patch->adds_atoms( *rsd_type );
	for ( Size m = 1; m <= new_atom_names.size(); m++ ) {
		if ( atom_names_soft_.has_value( ObjexxFCL::strip_whitespace( new_atom_names[ m ] ) ) ) return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This may not work with custom variants. Crap.
ResidueTypeCOPs
ResidueTypeFinder::filter_all_variants_matched( ResidueTypeCOPs const & rsd_types, bool const allow_extra_variants /* = false */ ) const
{
	ResidueTypeCOPs filtered_rsd_types = rsd_types;
	filtered_rsd_types = check_candidate_has_all_variant_sets( filtered_rsd_types );
	if ( !allow_extra_variants ) {
		filtered_rsd_types = check_variant_sets_have_all_candidate_variants( filtered_rsd_types );
	}
	return filtered_rsd_types;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::check_candidate_has_all_variant_sets( ResidueTypeCOPs const & rsd_types ) const
{
	ResidueTypeCOPs filtered_rsd_types;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP rsd_type = rsd_types[ n ];

		vector1< bool > variant_list_found_partner( rsd_type->variant_types().size(), false );

		bool all_variant_sets_matched( true );
		for ( Size k = 1; k <= variants_in_sets_.size(); k++ ) {
			vector1< VariantType > const & variant_set = variants_in_sets_[ k ];
			if ( !match_adducts_ && variant_set.size() == 1 && variant_set[ 1 ] == ADDUCT_VARIANT ) continue;

			bool at_least_one_variant_matched( false );
			for ( Size m = 1; m <= variant_set.size(); m++ ) {
				if ( rsd_type->properties().is_variant_type( variant_set[ m ] ) ) {
					at_least_one_variant_matched = true;
					break;
				}
			} // loop inside one variant set

			if ( !at_least_one_variant_matched ) {
				all_variant_sets_matched = false;
				break;
			}
		}
		if ( all_variant_sets_matched ) {
			filtered_rsd_types.push_back( rsd_type );
		}
	}
	return filtered_rsd_types;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::check_variant_sets_have_all_candidate_variants( ResidueTypeCOPs const & rsd_types ) const
{
	ResidueTypeCOPs filtered_rsd_types;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP rsd_type = rsd_types[ n ];
		utility::vector1< std::string > const & variant_types = rsd_type->variant_types();

		bool all_candidate_variant_types_matched = true;
		for ( Size q = 1; q <= variant_types.size(); q++ ) {

			VariantType candidate_variant_type = ResidueProperties::get_variant_from_string( variant_types[ q ] );
			if ( candidate_variant_type == NO_VARIANT ) { // check in custom variants
				if ( custom_variants_.has_value( variant_types[ q ] ) ) continue;
				all_candidate_variant_types_matched = false;
				break;
			}

			if ( !match_adducts_ && candidate_variant_type == ADDUCT_VARIANT ) continue;
			bool at_least_one_variant_matched( false );
			for ( Size k = 1; k <= variants_in_sets_.size(); k++ ) {
				vector1< VariantType > const & variant_set = variants_in_sets_[ k ];
				for ( Size m = 1; m <= variant_set.size(); m++ ) {
					if ( variant_set[m] == candidate_variant_type ) {
						at_least_one_variant_matched = true;
						break;
					}
				} // loop inside one variant set
				if ( at_least_one_variant_matched ) break;
			}

			if ( !at_least_one_variant_matched ) {
				all_candidate_variant_types_matched = false;
				break;
			}

		}

		if ( all_candidate_variant_types_matched ) {
			filtered_rsd_types.push_back( rsd_type );
		}

	}
	return filtered_rsd_types;
}

////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_all_patch_names( ResidueTypeCOPs const & rsd_types )  const
{
	ResidueTypeCOPs filtered_rsd_types;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP const & rsd_type = rsd_types[ n ];
		bool found_patch_name( true );
		for ( Size m = 1; m <= patch_names_.size(); m++ )  {
			std::string const & patch_name = patch_names_[ m ];
			if ( residue_type_all_patches_name( *rsd_type ).find( patch_name ) == std::string::npos ) {
				found_patch_name = false;
				break;
			}
		}
		if ( !found_patch_name ) continue;
		filtered_rsd_types.push_back( rsd_type );
	}
	return filtered_rsd_types;
}
////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_all_properties( ResidueTypeCOPs const & rsd_types )  const
{
	ResidueTypeCOPs filtered_rsd_types;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP const & rsd_type = rsd_types[ n ];
		bool found_property( true );
		for ( Size m = 1; m <= properties_.size(); m++ )  {
			if ( !rsd_type->has_property( properties_[ m ] ) ) {
				found_property = false;
				break;
			}
		}
		if ( !found_property ) continue;
		filtered_rsd_types.push_back( rsd_type );
	}
	return filtered_rsd_types;
}
////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_disallow_variants( ResidueTypeCOPs const & rsd_types )  const
{
	ResidueTypeCOPs filtered_rsd_types;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP const & rsd_type = rsd_types[ n ];
		bool disallowed( false );
		for ( Size n = 1; n <= disallow_variants_.size(); n++ ) {
			if ( rsd_type->has_variant_type( disallow_variants_[ n ] ) ) {
				disallowed = true; break;
			}
		}
		if ( !disallowed ) filtered_rsd_types.push_back( rsd_type );
	}
	return filtered_rsd_types;
}
////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_disallow_properties( ResidueTypeCOPs const & rsd_types )  const
{
	ResidueTypeCOPs filtered_rsd_types;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP const & rsd_type = rsd_types[ n ];
		bool disallowed( false );
		for ( Size n = 1; n <= disallow_properties_.size(); n++ ) {
			if ( rsd_type->has_property( disallow_properties_[ n ] ) ) {
				disallowed = true; break;
			}
		}
		if ( !disallowed ) filtered_rsd_types.push_back( rsd_type );
	}
	return filtered_rsd_types;
}
////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_special_cases( ResidueTypeCOPs const & rsd_types )  const
{
	ResidueTypeCOPs filtered_rsd_types;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP const & rsd_type = rsd_types[ n ];
		if ( disallow_carboxyl_conjugation_at_glu_asp_ &&
				( rsd_type->aa() == aa_glu || rsd_type->aa() == aa_asp ) &&
				rsd_type->has_variant_type( BRANCH_LOWER_TERMINUS_VARIANT ) ) continue;
		filtered_rsd_types.push_back( rsd_type );
	}
	return filtered_rsd_types;
}
////////////////////////////////////////////////////////////////////
/// @brief   set function for variants
/// @details actually updates variants_in_sets_.
ResidueTypeFinder &
ResidueTypeFinder::variants( utility::vector1< VariantType > const & setting )
{
	variants_in_sets_.clear();
	for ( Size n = 1; n <= setting.size(); n++ ) {
		variants_in_sets_.push_back( make_vector1( setting[ n ] ) );
	}

	return *this;
}

////////////////////////////////////////////////////////////////////
/// @brief   set function for variants
/// @details actually updates variants_in_sets_.
ResidueTypeFinder &
ResidueTypeFinder::variants( utility::vector1< std::string > const & setting )
{
	vector1< VariantType > variants_;
	for ( Size n = 1; n <= setting.size(); n++ ) {
		VariantType variant_type = ResidueProperties::get_variant_from_string( setting[ n ] );
		if ( variant_type != NO_VARIANT ) {
			variants_.push_back( variant_type );
		} else {
			custom_variants_.push_back( setting[ n ] );
		}
	}

	return variants( variants_ );
}


} //chemical
} //core
