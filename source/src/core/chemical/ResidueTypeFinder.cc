// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/ResidueTypeFinder.hh
/// @brief Functions to find residue_type(s) from ResidueTypeSet without requiring instantiation of all rsd_types.
/// @details Intended to be super-efficient replacement for aa_map(), name3_map() in ResidueTypeSet.
/// @author Rhiju Das, rhiju@stanford.edu
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- cleaned this up a bit and made it more efficient as part of the PackerPalette implementation.

#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/Metapatch.hh>
#include <core/chemical/util.hh>

#include <basic/options/option.hh>

#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/string.functions.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>

using namespace basic::options;
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
	name3_(""),
	residue_type_base_name_(""),
	base_type_(),
	interchangeability_group_(""),
	base_property_( NO_PROPERTY ),
	ignore_atom_named_H_( false ),
	disallow_carboxyl_conjugation_at_glu_asp_( false ),
	check_nucleic_acid_virtual_phosphates_( false ),
	no_metapatches_(false),
	no_CCD_on_name3_match_( ! option[ OptionKeys::in::file::check_all_PDB_components ]() )
{}

//Destructor
ResidueTypeFinder::~ResidueTypeFinder() = default;

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOP
ResidueTypeFinder::get_representative_type( bool const metapatches ) const
{
	ResidueTypeCOPs rsd_types;
	rsd_types = get_possible_base_residue_types( false /* include_unpatchable */ );

	rsd_types = apply_patches_recursively( rsd_types, 1 /*start with this patch*/, true /*get_first_residue_found*/ );
	if ( metapatches && !no_metapatches() ) {
		rsd_types = apply_metapatches_recursively( rsd_types, 1 /*start with this patch*/ );
		// We need to apply metapatches again just in case there are some double variants.
		// Only needed for packing metapatched residues.
		// TODO: this is atrocious.
		rsd_types = apply_metapatches_recursively( rsd_types, 1 /*start with this patch*/ );
	}
	if ( rsd_types.size() == 0 ) rsd_types =  get_possible_unpatchable_residue_types();

	rsd_types = apply_filters_after_patches( rsd_types, true /* allow_extra_variants */ );

	if ( rsd_types.size() == 0 ) return nullptr;
	return rsd_types[ 1 ];
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::get_all_possible_residue_types( bool const allow_extra_variants /* = false */ ) const
{
	// Get all possible basic residues that might match
	ResidueTypeCOPs rsd_types = get_possible_base_residue_types( false /* include_unpatchable*/ );

	TR.Debug << "Found " << rsd_types.size() << " base ResidueTypes." << std::endl;

	// Go down the binary tree of patches.
	rsd_types = apply_patches_recursively( rsd_types, 1 /*start with this patch*/ );

	if ( !no_metapatches() ) {
		rsd_types = apply_metapatches_recursively( rsd_types, 1 /*start with this patch*/ );
		// We need to apply metapatches again just in case there are some double variants.
		// Only needed for packing metapatched residues.
		// TODO: this is atrocious.
		rsd_types = apply_metapatches_recursively( rsd_types, 1 /*start with this patch*/ );
	}

	// add in any unpatchable residues.
	rsd_types.append( get_possible_unpatchable_residue_types() );

	TR.Debug << "Patched up to " << rsd_types.size() << " ResidueTypes." << std::endl;

	// Filter for rsd_types that strictly obey requirements
	rsd_types = apply_filters_after_patches( rsd_types, allow_extra_variants );

	rsd_types = apply_preferences_and_discouragements( rsd_types );

	return rsd_types;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOP
ResidueTypeFinder::get_best_match_residue_type_for_atom_names( utility::vector1< std::string > const & atom_names )
{
	//clock_t const time_start( clock() );

	// will try to match these ('soft' constraints). Go ahead and strip out whitespace.
	atom_names_soft_.clear();
	for ( Size n = 1; n <= atom_names.size(); ++n ) {
		std::string atom_name_temp = atom_names[ n ];
		atom_names_soft_.push_back( ObjexxFCL::strip_whitespace( atom_name_temp ) );
	}

	ResidueTypeCOPs rsd_types = get_all_possible_residue_types( true /* allow_extra_variants */ );

	Size const n_types( rsd_types.size() );
	if ( n_types == 0 ) {
		return nullptr;
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Finding best match from among " << n_types << " ResidueTypes." << std::endl;
		for ( uint i( 1 ); i <= n_types; ++i ) {
			TR.Trace << ' ' << rsd_types[ i ]->name() << std::endl;
		}
	}

	ResidueTypeCOP rsd_type = find_best_match( rsd_types, atom_names, ignore_atom_named_H_ );

	//TR << "time to initialize " << rsd_type->name() << " from " << rsd_types.size() << " possible ResidueTypes: " <<
	// static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;

	return rsd_type;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::get_possible_base_residue_types(
	bool const include_unpatchable /*=true*/, bool const apply_all_filters /*=false*/ ) const
{
	if ( base_type_ ) { //If a base type has already been specified, there's no need to bother with a lot of other rigamarole.
		ResidueTypeCOPs rsd_types;
		rsd_types.push_back( base_type_ );
		return rsd_types;
	}

	//Otherwise, load the whole set of base types and start pruning:
	initialize_relevant_pdb_components();

	ResidueTypeCOPs rsd_types = residue_type_set_.base_residue_types();


	if ( include_unpatchable ) {
		rsd_types.append( get_possible_base_unpatchable_residue_types() );
	}
	rsd_types = apply_basic_filters( rsd_types );
	if ( apply_all_filters ) {
		rsd_types = apply_filters_after_patches( rsd_types );
		rsd_types = apply_preferences_and_discouragements( rsd_types );
	}
	return rsd_types;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// @brief Load the PDB component into the RTS's base residue type list, if we're specifying a name3 search
void
ResidueTypeFinder::initialize_relevant_pdb_components() const
{
	if ( name3_.empty() ) { return; } // Only take components for name3 searches

	if ( no_CCD_on_name3_match_ ) {
		// if there's already a Rosetta type with the name3 in the base residue types, don't bother loading the component.
		for ( auto rsd_type : residue_type_set_.base_residue_types() ) {
			if ( rsd_type->name3() == name3_ ||
					residue_type_set_.generates_patched_residue_type_with_name3( residue_type_base_name( *rsd_type ), name3_ ) ) {
				return; // Don't bother with component loading.
			}
		}
	}

	// Will cause the RTS to load the PDB component into its base residue types
	residue_type_set_.name_mapOP( "pdb_" + utility::strip(name3_) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::get_possible_unpatchable_residue_types() const
{
	ResidueTypeCOPs rsd_types = residue_type_set_.unpatchable_residue_types();
	rsd_types = apply_basic_filters( rsd_types );
	return rsd_types;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::get_possible_base_unpatchable_residue_types() const
{
	ResidueTypeCOPs filtered_rsd_types;
	ResidueTypeCOPs rsd_types = residue_type_set_.unpatchable_residue_types();
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP rsd_type = rsd_types[ n ];
		if ( residue_type_base_name( *rsd_type ) == rsd_type->name() ) {
			filtered_rsd_types.push_back( rsd_type );
		}
	}
	return filtered_rsd_types;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::apply_basic_filters( ResidueTypeCOPs rsd_types ) const
{
	if ( rsd_types.empty() ) {
		TR.Debug << "Going into apply_basic_filters() ResidueType filtering, no residue types were passed." << std::endl;
		return rsd_types;
	}
	rsd_types = filter_by_aa( rsd_types );
	rsd_types = filter_by_name1( rsd_types );
	rsd_types = filter_by_name3( rsd_types, true /* keep_if_base_type_generates_name3 */ );
	rsd_types = filter_by_residue_type_base_name( rsd_types );
	rsd_types = filter_by_interchangeability_group( rsd_types, true /* keep_if_base_type_generates_interchangeability_group */ );
	rsd_types = filter_by_base_property( rsd_types );
	return rsd_types;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::apply_filters_after_patches( ResidueTypeCOPs rsd_types,
	bool const allow_extra_variants  /* = false */ ) const
{
	if ( rsd_types.empty() ) {
		TR.Debug << "Going into apply_filters_after_patches() ResidueType filtering, no residue types were passed." << std::endl;
		return rsd_types;
	}
	rsd_types = filter_by_name3( rsd_types, false /* keep_if_base_type_generates_name3 */ );
	rsd_types = filter_by_interchangeability_group( rsd_types, false /* keep_if_base_type_generates_interchangeability */ );
	rsd_types = filter_disallow_variants( rsd_types );
	rsd_types = filter_all_variants_matched( rsd_types, allow_extra_variants );
	rsd_types = filter_all_properties( rsd_types );
	rsd_types = filter_disallow_properties( rsd_types );
	rsd_types = filter_all_patch_names( rsd_types );
	rsd_types = filter_connections( rsd_types );
	rsd_types = filter_special_cases( rsd_types );

	return rsd_types;
}

////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::apply_preferences_and_discouragements( ResidueTypeCOPs const & rsd_types ) const {
	if ( rsd_types.empty() ) return rsd_types;

	if ( preferred_properties_.empty() && discouraged_properties_.empty() && ! no_CCD_on_name3_match_ ) {
		return rsd_types; // nothing to do
	}

	ResidueTypeCOPs current_type_list( rsd_types );
	ResidueTypeCOPs new_type_list;

	// Discouragement before encouragment

	// keep the ResidueTypes in order
	if ( ! discouraged_properties_.empty() ) {
		utility::vector1< core::Size > property_counts;
		for ( ResidueTypeCOP const & rsd_type: current_type_list ) {
			core::Size prop_count = 0;
			for ( ResidueProperty prop: discouraged_properties_ ) {
				if ( rsd_type->has_property( prop ) ) {
					++prop_count;
				}
			}
			property_counts.push_back( prop_count );
		}
		debug_assert( ! property_counts.empty() );
		core::Size min_count = *std::min_element( property_counts.begin(), property_counts.end() );
		new_type_list.clear();
		for ( core::Size ii(1); ii <= current_type_list.size(); ++ii ) {
			if ( property_counts[ ii ] == min_count ) {
				new_type_list.push_back( current_type_list[ ii ] );
			}
		}
		if ( TR.Debug.visible() ) {
			TR.Debug << "Discouraging " << discouraged_properties_.size() << " properties, "
				<< "going from " << current_type_list.size() << " types to " << new_type_list.size() << " types." << std::endl;
			TR.Debug << "Discouraged: " <<  discouraged_properties_ << std::endl;
			TR.Debug << "Going from ";
			for ( auto rt: current_type_list ) { TR.Debug << " " << rt->name(); }
			TR.Debug << std::endl;
			TR.Debug << "To ";
			for ( auto rt: new_type_list ) { TR.Debug << " " << rt->name(); }
			TR.Debug << std::endl;
		}
		current_type_list = new_type_list;
	}

	if ( ! preferred_properties_.empty() ) {
		utility::vector1< core::Size > property_counts;
		for ( ResidueTypeCOP const & rsd_type: current_type_list ) {
			core::Size prop_count = 0;
			for ( ResidueProperty prop: preferred_properties_ ) {
				if ( rsd_type->has_property( prop ) ) {
					++prop_count;
				}
			}
			property_counts.push_back( prop_count );
		}
		debug_assert( ! property_counts.empty() );
		core::Size max_count = *std::max_element( property_counts.begin(), property_counts.end() );
		new_type_list.clear();
		for ( core::Size ii(1); ii <= current_type_list.size(); ++ii ) {
			if ( property_counts[ ii ] == max_count ) {
				new_type_list.push_back( current_type_list[ ii ] );
			}
		}
		if ( TR.Debug.visible() ) {
			TR.Debug << "Encouraging " << preferred_properties_.size() << " properties, "
				<< "going from " << current_type_list.size() << " types to " << new_type_list.size() << " types." << std::endl;
			TR.Debug << "Encouraged: " << preferred_properties_ << std::endl;
			TR.Debug<< "Going from ";
			for ( auto rt: current_type_list ) { TR.Debug << " " << rt->name(); }
			TR.Debug << std::endl;
			TR.Debug << "To ";
			for ( auto rt: new_type_list ) { TR.Debug << " " << rt->name(); }
			TR.Debug << std::endl;
		}

		current_type_list = new_type_list;
	}

	current_type_list = prioritize_rosetta_types_over_pdb_components( current_type_list );

	return current_type_list;
}

////////////////////////////////////////////////////////
// @brief if no_CCD_on_name3_match_ is on, then any pdb components should be ignored if we have an otherwise qualifying
// Rosetta database file, even if there's no chemical match.
// (Note that the components which are chemically equivalent to Rosetta types should be handled by the
// exclude_pdb_component_list.txt file in the database, along with the exclude_pdb_component_ids_ functionality in
// the GlobalResidueType set.)
ResidueTypeCOPs
ResidueTypeFinder::prioritize_rosetta_types_over_pdb_components( ResidueTypeCOPs const & rsd_types ) const
{
	if ( ! no_CCD_on_name3_match_ ) { return rsd_types; } // Shouldn't do any filtering.

	ResidueTypeCOPs filtered_rsd_types;
	for ( auto const & rsd_type : rsd_types ) {
		if ( rsd_type->name().size() < 4 || !( rsd_type->name().substr(0,4) == "pdb_" ) ) {
			filtered_rsd_types.push_back( rsd_type );
		}
	}

	if ( filtered_rsd_types.empty() ) {
		return rsd_types; // Only have the components.
	} else {
		return filtered_rsd_types;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Instantiates ResidueType (and gets exponentially slower), based on the number
///  of desired features that a patch might offer:
///
/// (1) xyz_atom_names   Does patch introduce this atom via ADD_ATOM?
/// (2) variants         Does patch introduce this VariantType through its TYPES?
/// (3) properties       Does patch introduce this ResidueProperty through ADD_PROPERTY?
/// (4) patch_names      Used for carbohydrate branching --> look explicitly for name.
///
/// This function will take less time and memory if the lists above are *short*.
///
/// So you can be smart about tuning performance -- for example, if you are looking for "->2)-branch"
///  carbohydrates, set that through patch_names; and do not also look for BRANCH_POINT under properties (which
///  is also a property of patches applied to amino acids).
///
///                 -- rhiju, 2015
///
/// @details This function made extensive use of string parsing, which is very inefficient,
/// particularly in a recursive context.  Cleaned up considerably by V. Mulligan on 18 Aug 2016.
///
////////////////////////////////////////////////////////////////////////////////////////////////////
vector1< ResidueTypeCOP >
ResidueTypeFinder::apply_patches_recursively(
	vector1< ResidueTypeCOP > const & rsd_types,
	Size const patch_number,
	bool const get_first_totally_ok_residue_type /*= false*/
) const {

	// Pointless to apply patches if we don't have any residue types to apply them to.
	if ( rsd_types.empty() ) {
		return rsd_types;
	}

	ResidueTypeCOPs rsd_types_new = rsd_types;
	PatchCOP patch = residue_type_set_.patches()[ patch_number ];

	for ( auto const & rsd_type : rsd_types ) {

		// absolute no-no's.
		if ( !patch->applies_to( *rsd_type ) )             continue;
		if ( has_disallowed_variant( patch ) )             continue;
		if ( deletes_any_property(   patch, rsd_type ) )   continue;
		if ( deletes_any_variant(    patch, rsd_type ) )   continue;
		if ( changes_to_wrong_aa(    patch, rsd_type ) )   continue;
		// could also add as a no-no: if patch *deletes* an atom in atom_names_.

		// note -- make sure to apply patch if it has a chance of satisfying any of
		// the constraints on variants, branchpoints, or properties.
		bool apply_patch = (  adds_any_variant( patch ) ||
			adds_any_property( patch, rsd_type ) ||
			matches_any_patch_name( patch )      ||
			matches_any_atom_name( patch, rsd_type ) ||
			fixes_name3( patch, rsd_type ) ||
			fixes_interchangeability_group( patch, rsd_type ) ||
			fixes_connects( patch, rsd_type ) );

		if ( apply_patch ) {
			// following just gets the right name of the patched residue
			std::string const & patched_name( patch->patched_name( *rsd_type ) );
			// by using name_map, forces residue_type_set to generate the real residue_type, cache it, and return the COP:
			ResidueTypeCOP rsd_type_new( residue_type_set_.name_mapOP( patched_name ) );
			if ( rsd_type_new ) {
				rsd_types_new.push_back( rsd_type_new );
			}
		}

	} // end loop

	if ( get_first_totally_ok_residue_type && ! rsd_types_new.empty() ) { // maybe we're done?
		// note that this repeats some work -- some rsd_types were checked in prior steps in the recursion
		ResidueTypeCOPs rsd_types_filtered = apply_filters_after_patches( rsd_types_new, true /*allow_extra_variants*/ );
		if ( ! rsd_types_filtered.empty() ) return rsd_types_filtered;
	}

	// end of recursion through patches?
	if ( patch_number == residue_type_set_.patches().size() ) return rsd_types_new;

	return  apply_patches_recursively( rsd_types_new, patch_number + 1, get_first_totally_ok_residue_type );

}

vector1< ResidueTypeCOP >
ResidueTypeFinder::apply_metapatches_recursively(
	vector1< ResidueTypeCOP > const & rsd_types,
	Size const metapatch_number,
	bool const get_first_totally_ok_residue_type /*= false*/
) const {
	if ( no_metapatches() ) return rsd_types;
	ResidueTypeCOPs rsd_types_new = rsd_types;
	utility::vector1< MetapatchCOP > metapatch_list( residue_type_set_.metapatches() ); // Returned by value
	if ( metapatch_number == 0 || metapatch_number > metapatch_list.size() ) return rsd_types_new;

	MetapatchCOP metapatch = metapatch_list[ metapatch_number ];

	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP rsd_type = rsd_types[ n ];

		utility::vector1< std::string > atoms = metapatch->atoms( *rsd_type );
		for ( Size i = 1; i <= atoms.size(); ++i ) {
			PatchCOP patch = metapatch->get_one_patch( /**rsd_type,*/ atoms[ i ] );

			//TR << "Considering " << rsd_type->name() << " for " << patch->name() << std::endl;

			// absolute no-no's.
			if ( !patch->applies_to( *rsd_type ) )             continue;
			if ( has_disallowed_variant( patch ) )             continue;
			if ( deletes_any_property(   patch, rsd_type ) )   continue;
			if ( deletes_any_variant(    patch, rsd_type ) )   continue;
			if ( changes_to_wrong_aa(    patch, rsd_type ) )   continue;
			// could also add as a no-no: if patch *deletes* an atom in atom_names_.
			//TR << "Passed continue " << rsd_type->name() << " for " << patch->name() << std::endl;

			// note -- make sure to apply patch if it has a chance of satisfying any of
			// the constraints on variants, branchpoints, or properties.
			bool apply_patch = (  adds_any_variant( patch ) ||
				adds_any_property( patch, rsd_type ) ||
				matches_any_patch_name( patch )      ||
				matches_any_atom_name( patch, rsd_type ) ||
				fixes_name3( patch, rsd_type ) ||
				fixes_interchangeability_group( patch, rsd_type ) // ||
				// fixes_connects( patch, rsd_type ) // Currently an issue, as the connect metapatch plays havoc
			);

			if ( apply_patch ) {
				//TR << "Considering " << rsd_type->name() << " for " << patch->name() << std::endl;

				// following just gets the right name of the patched residue
				ResidueTypeCOP rsd_type_new_placeholder = patch->apply( *rsd_type, false /*instantiate*/ );
				if ( rsd_type_new_placeholder ) {
					// by using name_map, forces residue_type_set to generate the real residue_type, cache it, and return the COP:
					ResidueTypeCOP rsd_type_new( residue_type_set_.name_mapOP( rsd_type_new_placeholder->name() ) );
					if ( rsd_type_new ) {
						rsd_types_new.push_back( rsd_type_new );
					}
				}
			}
		}

	} // end loop

	if ( get_first_totally_ok_residue_type ) { // maybe we're done?
		// note that this repeats some work -- some rsd_types were checked in prior steps in the recursion
		ResidueTypeCOPs rsd_types_filtered = apply_filters_after_patches( rsd_types_new, true /*allow_extra_variants*/ );
		if ( rsd_types_filtered.size() > 0 ) return rsd_types_filtered;
	}

	// end of recursion through patches?
	if ( metapatch_number == residue_type_set_.metapatches().size() ) return rsd_types_new;

	return  apply_metapatches_recursively( rsd_types_new, metapatch_number + 1, get_first_totally_ok_residue_type );
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
	if ( rsd_types_new.empty() && ! rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering by one letter code: '" << name1_ << "'" << std::endl;
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
				utility::strip( rsd_type->name3() ) == utility::strip( name3_ ) || // name3 may be whitespace padded
				( keep_if_base_type_generates_name3 && residue_type_set_.generates_patched_residue_type_with_name3( rsd_type->name(), name3_ ) ) ) {
			rsd_types_new.push_back( rsd_type );
		}
	}
	if ( rsd_types_new.empty() && ! rsd_types.empty() ) {
		// RM: This seems to occur frequently in the integration tests.
		// We may want to examine why we're doing so much filtering by name3 on sets that don't contain what we want.
		TR.Debug << "No ResidueTypes remain after filtering by three letter code: '" << name3_ << "'" << std::endl;
	}
	return rsd_types_new;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_by_interchangeability_group( ResidueTypeCOPs const & rsd_types, bool const keep_if_base_type_generates_interchangeability_group  ) const
{
	if ( interchangeability_group_.size() == 0 ) return rsd_types;
	ResidueTypeCOPs rsd_types_new;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP rsd_type( rsd_types[ n ] );
		if ( rsd_type->interchangeability_group() == interchangeability_group_ ||
				( keep_if_base_type_generates_interchangeability_group &&
				residue_type_set_.generates_patched_residue_type_with_interchangeability_group( rsd_type->name(), interchangeability_group_ ) ) ) {
			rsd_types_new.push_back( rsd_type );
		}
	}
	if ( rsd_types_new.empty() && ! rsd_types.empty() ) {
		// RM: This seems to occur frequently in the integration tests.
		// We may want to examine why we're doing so much filtering on sets that don't contain what we want.
		TR.Debug << "No ResidueTypes remain after filtering by interchangeability group: '" << interchangeability_group_ << "'" << std::endl;
	}
	return rsd_types_new;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_by_aa( ResidueTypeCOPs const & rsd_types ) const
{
	if ( aa_ == aa_none ) return rsd_types;

	ResidueTypeCOPs rsd_types_new;
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP rsd_type( rsd_types[ n ] );
		if ( rsd_type->aa() == aa_ ) {
			rsd_types_new.push_back( rsd_type );
		}
	}
	if ( rsd_types_new.empty() && ! rsd_types.empty() ) {
		// RM: This seems to occur frequently in the integration tests.
		// We may want to examine why we're doing so much filtering by aa_type on sets that don't contain what we want.
		TR.Debug << "No ResidueTypes remain after filtering by amino acid designation: '" << aa_ << "'" << std::endl;
	}
	return rsd_types_new;
}

////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_by_residue_type_base_name( ResidueTypeCOPs const & rsd_types ) const
{
	if ( residue_type_base_name_.size() == 0 ) return rsd_types;

	ResidueTypeCOPs filtered_rsd_types;
	for ( ResidueTypeCOP const & rsd_type : rsd_types ) {
		std::string base_name( residue_type_base_name( *rsd_type ) );
		if ( base_name == residue_type_base_name_ || rsd_type->base_name() == residue_type_base_name_ ) {
			filtered_rsd_types.push_back( rsd_type );
		}
	}
	if ( filtered_rsd_types.empty() && ! rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering by ResidueType base name: '" << residue_type_base_name_ << "'" << std::endl;
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
	if ( filtered_rsd_types.empty() && ! rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering by ResidueType base property: '" << base_property_ << "'" << std::endl;
	}
	return filtered_rsd_types;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResidueTypeFinder::has_disallowed_variant( PatchCOP patch ) const
{
	vector1< core::chemical::VariantType > const & patch_variant_types = patch->types();
	for ( Size k(1), kmax(patch_variant_types.size()); k <= kmax; ++k )  {
		if ( disallow_variants_.has_value( patch_variant_types[ k ] ) ) return true;
		//Does not currently support on-the-fly types.  --VKM.
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResidueTypeFinder::adds_any_variant( PatchCOP patch ) const
{
	if ( !patch->types().empty() ) {
		vector1< core::chemical::VariantType > const & patch_variant_types = patch->types();
		for ( auto const & patch_variant_type : patch_variant_types ) {
			for ( Size k = 1; k <= variants_in_sets_.size(); k++ ) {
				if ( variants_in_sets_[ k ].has_value( patch_variant_type ) ) return true;
				if ( variant_exceptions_.has_value( patch_variant_type ) ) return true; // explore all of these 'exceptions' (used for adducts)
			}

			// Since PDB-input matching tends to aggressively favor justifying patches
			// with the EXISTENCE of atoms, we could get around the below by checking
			// if a patch virtualizes any atoms.
			if ( check_nucleic_acid_virtual_phosphates_ &&
					( patch_variant_type == VIRTUAL_DNA_PHOSPHATE || patch_variant_type == VIRTUAL_PHOSPHATE ) &&
					!atom_names_soft_.has_value( "P" ) ) return true;
		}
	}

	if ( !patch->custom_types().empty() ) {
		vector1< std::string> const & patch_custom_variant_types( patch->custom_types() );
		for ( auto const & custom_variant : custom_variants_ ) {
			for ( auto const & patch_custom_variant : patch_custom_variant_types ) {
				if ( custom_variant == patch_custom_variant ) return true;
				//Andy, you had the following, but it was almost certainly not doing what you thought it was doing. --VKM
				//if ( custom_variants_[ n ].substr(0, custom_variants_.size()-1) == patch_variant_types[ k ] ) return true;
				//if ( custom_variants_[ n ].substr(0, custom_variants_.size()-2) == patch_variant_types[ k ] ) return true;
			}
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
			residue_type_set_.generates_patched_residue_type_with_name3( residue_type_base_name( *rsd_type ), name3_ ) ) {
		ResidueTypeCOP new_type( patch->apply( *rsd_type, false /*instantiate*/ ) );
		if ( new_type && new_type->name3() == name3_ ) {
			return true;
		}
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResidueTypeFinder::fixes_interchangeability_group( PatchCOP patch, ResidueTypeCOP rsd_type ) const
{
	if ( interchangeability_group_.size() > 0 &&
			rsd_type->interchangeability_group() != interchangeability_group_ &&
			residue_type_set_.generates_patched_residue_type_with_interchangeability_group( residue_type_base_name( *rsd_type ),
			interchangeability_group_ ) ) {
		ResidueTypeCOP new_type( patch->apply( *rsd_type, false /*instantiate*/ ) );
		if ( new_type && new_type->interchangeability_group() == interchangeability_group_ ) {
			return true;
		}
		return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResidueTypeFinder::fixes_connects( PatchCOP patch, ResidueTypeCOP rsd_type ) const {
	if ( connect_atoms_.empty() ) return false; // Can't fix what isn't broken.
	for ( std::string const & atom: connect_atoms_ ) {
		if ( rsd_type->has(atom) ) {
			// patch->changes_connections_on() should be whitespace padding insensitive.
			if ( rsd_type->residue_connections_for_atom( rsd_type->atom_index(atom) ).empty() &&
					patch->changes_connections_on( *rsd_type, atom ) ) {
				return true;
			}
		} else {
			// Don't have the atom -- get patches which may add the atom.
			if ( patch->adds_atoms( *rsd_type ).has_value( atom ) ) {
				return true;
			}
		}
	}
	return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Returns true if this patch adds any of the properties that the ResidueTypeFinder is seeking.
/// @details ONLY works for canonical properties (not on-the-fly properties) at present.  Modified on
/// 24 Aug 2016 by VKM to remove string parsing.
bool
ResidueTypeFinder::adds_any_property( PatchCOP patch, ResidueTypeCOP rsd_type ) const
{
	utility::vector1< ResidueProperty > const added_properties( patch->adds_properties_enums( *rsd_type ) );

	for ( auto const & prop : added_properties ) {
		if ( properties_.has_value( prop ) ) return true;
	}

	// Do we need to check for preferred properties now? Why didn't we previously?
	for ( auto const & prop : added_properties ) {
		if ( preferred_properties_.has_value( prop ) ) return true;
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Returns ture if this patch deletes any of the properties that the ResidueTypeFinder is seeking.
/// @details ONLY works for canonical properties (not on-the-fly properties) at present.  Modified on
/// 25 Aug 2016 by VKM to remove string parsing.
bool
ResidueTypeFinder::deletes_any_property( PatchCOP patch, ResidueTypeCOP rsd_type ) const
{
	utility::vector1< ResidueProperty > const deleted_properties( patch->deletes_properties_enums( *rsd_type ) );
	for ( auto const & prop : deleted_properties ) {
		if ( properties_.has_value( prop ) ) return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResidueTypeFinder::changes_to_wrong_aa( PatchCOP patch, ResidueTypeCOP rsd_type ) const
{
	if ( aa_ == aa_none ) return false;
	AA new_aa = patch->generates_aa( *rsd_type );
	if ( new_aa != aa_none && new_aa != aa_  ) return true;
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResidueTypeFinder::deletes_any_variant( PatchCOP patch, ResidueTypeCOP rsd_type ) const
{
	vector1< VariantType > const & deleted_variants_by_enum( patch->deletes_variants_by_enum( *rsd_type ) );
	for ( auto const & deleted_variant : deleted_variants_by_enum ) {
		for ( auto const & variant_set : variants_in_sets_ ) {
			if ( variant_set.has_value( deleted_variant ) ) return true;
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
	for ( auto & new_atom_name : new_atom_names ) {
		if ( atom_names_soft_.has_value( ObjexxFCL::strip_whitespace( new_atom_name ) ) ) return true;
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
	if ( filtered_rsd_types.empty() && ! rsd_types.empty() ) {
		// RM: This seems to occur frequently in the integration tests.
		// We may want to examine why we're doing so much filtering by matched variants on sets that don't contain what we want.
		TR.Debug << "No ResidueTypes remain after filtering for matched variants." << std::endl;
		//TODO: Print what those matched variants should have been.
	}
	return filtered_rsd_types;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::check_candidate_has_all_variant_sets( ResidueTypeCOPs const & rsd_types ) const
{
	ResidueTypeCOPs filtered_rsd_types;
	for ( auto const & rsd_type : rsd_types ) {
		vector1< bool > variant_list_found_partner( rsd_type->variant_types().size(), false );

		bool all_variant_sets_matched( true );
		for ( auto const & variant_set : variants_in_sets_ ) {
			if ( variant_set.size() == 1 && variant_exceptions_.has_value( variant_set[ 1 ] ) ) continue;

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

		for ( auto const & custom_variant : custom_variants_ ) {
			if ( !rsd_type->properties().is_variant_type( custom_variant ) ) {
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
	for ( auto const & rsd_type : rsd_types ) {

		utility::vector1< core::chemical::VariantType > const & variant_types( rsd_type->variant_type_enums() );

		bool all_candidate_variant_types_matched( true );
		for ( auto const & variant_type : variant_types ) {
			if ( variant_exceptions_.has_value( variant_type ) ) continue;

			bool at_least_one_variant_matched( false );

			for ( auto const & variant_set : variants_in_sets_ ) {
				for ( auto const & elem : variant_set ) {
					if ( elem == variant_type ) {
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

		if ( !all_candidate_variant_types_matched ) continue; //Go on to the next if not everything matched

		//Now, we need to check on-the-fly types.  This and only this requires string parsing:
		utility::vector1 < std::string > const & custom_variant_types( rsd_type->custom_variant_types() );
		for ( auto const & custom_variant : custom_variant_types ) {
			if ( !custom_variants_.has_value( custom_variant ) ) {
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
	if ( TR.Debug.visible() && filtered_rsd_types.empty() && ! rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering for the presence of all patch names." << std::endl;
		TR.Debug << "    Patches sought: ";
		for ( Size m = 1; m <= patch_names_.size(); m++ )  { TR.Debug << patch_names_[ m ] << "   "; }
		TR.Debug << std::endl;
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
	if ( TR.Debug.visible() && filtered_rsd_types.empty() && ! rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering for the presence of all properties." << std::endl;
		TR.Debug << "    Properties sought: ";
		for ( Size m = 1; m <= properties_.size(); m++ )  { TR.Debug << ResidueProperties::get_string_from_property( properties_[ m ] ) << "   "; }
		TR.Debug << std::endl;
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
		for ( Size i = 1; i <= disallow_variants_.size(); ++i ) {
			if ( rsd_type->has_variant_type( disallow_variants_[ i ] ) ) {
				disallowed = true;
				break;
			}
		}
		if ( !disallowed ) filtered_rsd_types.push_back( rsd_type );
	}
	if ( TR.Debug.visible() && filtered_rsd_types.empty() && ! rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering against the presence of variants." << std::endl;
		TR.Debug << "    Variants prohibited: ";
		for ( Size m = 1; m <= disallow_variants_.size(); m++ )  { TR.Debug << ResidueProperties::get_string_from_variant( disallow_variants_[ m ] ) << "   "; }
		TR.Debug << std::endl;
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
		for ( Size i = 1; i <= disallow_properties_.size(); ++i ) {
			if ( rsd_type->has_property( disallow_properties_[ i ] ) ) {
				disallowed = true;
				break;
			}
		}
		if ( !disallowed ) filtered_rsd_types.push_back( rsd_type );
	}
	if ( TR.Debug.visible() && filtered_rsd_types.empty() && ! rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering against the presence of properties." << std::endl;
		TR.Debug << "    Properties prohibited: ";
		for ( Size m = 1; m <= disallow_properties_.size(); m++ )  { TR.Debug << ResidueProperties::get_string_from_property( disallow_properties_[ m ] ) << "   "; }
		TR.Debug << std::endl;
	}
	return filtered_rsd_types;
}

////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_connections( ResidueTypeCOPs const & rsd_types )  const
{
	if ( connect_atoms_.empty() || rsd_types.empty() ) { return rsd_types; }

	ResidueTypeCOPs filtered_rsd_types;
	for ( ResidueTypeCOP rsd: rsd_types ) {
		bool has_all_required_connections = true;
		for ( std::string const & atom: connect_atoms_ ) {
			if ( !rsd->has( atom ) ||
					rsd->residue_connections_for_atom( rsd->atom_index(atom) ).empty() ) {
				has_all_required_connections = false;
				break;
			}
		}
		if ( has_all_required_connections ) {
			filtered_rsd_types.push_back( rsd );
		}
	}
	return filtered_rsd_types;
}

////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::filter_special_cases( ResidueTypeCOPs const & rsd_types )  const
{
	ResidueTypeCOPs filtered_rsd_types;

	bool const actually_check_nucleic_acid_virtual_phosphates =
		check_nucleic_acid_virtual_phosphates_ && !atom_names_soft_.has_value( "P" );

	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		ResidueTypeCOP const & rsd_type = rsd_types[ n ];

		if ( actually_check_nucleic_acid_virtual_phosphates ) {
			if ( rsd_type->is_DNA() && !rsd_type->has_variant_type( VIRTUAL_DNA_PHOSPHATE ) ) continue;

			// triphosphate termination does not require virtual phosphates!
			if ( rsd_type->is_RNA() && !rsd_type->has_variant_type( VIRTUAL_PHOSPHATE )
					&& !rsd_type->has_variant_type( FIVE_PRIME_PACKABLE_TRIPHOSPHATE ) ) continue;
		}

		filtered_rsd_types.push_back( rsd_type );
	}
	if ( TR.Debug.visible() && filtered_rsd_types.empty() && ! rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering for special cases:" << std::endl;
		if ( actually_check_nucleic_acid_virtual_phosphates ) {
			TR.Debug << "    * Nucleic acid has virtual phosphates." << std::endl;
		}
	}
	return filtered_rsd_types;
}

////////////////////////////////////////////////////////////////////
/// @brief   set function for variants
/// @details actually updates variants_in_sets_.
ResidueTypeFinder &
ResidueTypeFinder::variants(
	utility::vector1< VariantType > const & setting,
	bool const clear_existing /*=true*/
) {
	if ( clear_existing ) {
		variants_in_sets_.clear();
	}
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

/// @brief Specify a list of standard variant types (by enum) and custom variant types (by string).
/// @details This is the most efficient way to handle variants, since it minimizes the string handling.  Everything that
/// can be handled by enum is handled by enum.
/// @param[in] std_variants A vector of enums of standard variants that the ResidueTypeFinder should match.
/// @param[in] custom_variants A vector of strings of custom variant types that the ResidueTypeFinder should match.  Note that
/// standard types should NOT be included in this list.  There is no check for this!
/// @param[in] clear_existing If true (default), the existing VariantType lists are cleared.  If false, this just appends to those lists.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
ResidueTypeFinder &
ResidueTypeFinder::variants(
	utility::vector1< VariantType > const & std_variants,
	utility::vector1< std::string > const & custom_variants,
	bool const clear_existing /*=true*/
) {
	if ( clear_existing ) {
		custom_variants_ = custom_variants;
	} else {
		for ( core::Size i=1, imax=custom_variants.size(); i<=imax; ++i ) {
			custom_variants_.push_back( custom_variants[i] );
		}
	}
	return variants( std_variants, clear_existing );
}

/// @brief Provide a list of VariantTypes that a matched ResidueType must NOT have.
/// @details By default, this overwrites the existing list.  To append to the existing list,
/// set clear_existing to false.
ResidueTypeFinder &
ResidueTypeFinder::disallow_variants(
	utility::vector1< VariantType > const & setting,
	bool const clear_existing /*=true*/
) {
	if ( clear_existing ) {
		disallow_variants_ = setting;
	} else {
		for ( core::Size i=1, imax=setting.size(); i<=imax; ++i ) {
			disallow_variants_.push_back( setting[i] );
		}
	}
	return *this;
}


////////////////////////////////////////////////////////////////////
/// @brief   set function for variant exceptions
ResidueTypeFinder &
ResidueTypeFinder::variant_exceptions( utility::vector1< std::string > const & setting, bool const clear_existing/*=true*/ )
{
	if ( clear_existing ) variant_exceptions_.clear();

	for ( Size n = 1; n <= setting.size(); n++ ) {
		VariantType variant_type = ResidueProperties::get_variant_from_string( setting[ n ] );
		if ( variant_type != NO_VARIANT ) {
			if ( !variant_exceptions_.has_value( variant_type ) ) {
				variant_exceptions_.push_back( variant_type );
			}
		} else {
			utility_exit_with_message( "not currently handling custom variants within variant_exceptions" );
		}
	}
	return *this;
}


/// @brief Provide a list of VariantTypes that will be ignored when matching.
/// @details By default, this overwritest the existing list.  To append to the existing list,
/// set clear_existing=false.
ResidueTypeFinder &
ResidueTypeFinder::variant_exceptions(
	utility::vector1< VariantType > const & setting,
	bool const clear_existing/*=true*/
) {
	if ( clear_existing ) {
		variant_exceptions_ = setting;
	} else {
		for ( core::Size i=1, imax=setting.size(); i<=imax; ++i ) {
			if ( !variant_exceptions_.has_value( setting[i] ) ) variant_exceptions_.push_back( setting[i] );
		}
	}
	return *this;
}

///// Saving this sketch of a class for only creating RTs if they would join the
///// pareto optimal set, i.e. they add atoms, or variants that are not covered by
///// any of the RTs that have been identified so far.
// class ResidueTypeQualifications
// {
// public:
//  ResidueTypeQualifications();
//  ResidueTypeQualifications( ResidueTypeQualifications const & );
//  ~ResidueTypeQualifications();
//
// private:
//  ResidueTypeCOP restype_;
//  utility::vector1< bool > variant_sets_satisfied_;
//  utility::vector1< bool > disallowed_variants_avoided_;
//  utility::vector1< bool > named_atoms_covered_;
//  utility::vector1< bool > properties_satisfied_;
//  utility::vector1< bool > disallowed_properties_avoided_;
//  utility::vector1< bool > patch_names_satified_;
//
// }


} //chemical
} //core
