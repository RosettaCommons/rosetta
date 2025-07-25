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
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/Metapatch.hh>
#include <core/chemical/util.hh>

#include <basic/options/option.hh>

#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/string.functions.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <algorithm>

#include <utility/stream_util.hh> // AUTO IWYU For operator<<

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
	no_metapatches_(false), // true would disable consideration of metapatches (relevant for speed)
	no_CCD_on_name3_match_( ! option[ OptionKeys::in::file::check_all_PDB_components ]() )
{}

//Destructor
ResidueTypeFinder::~ResidueTypeFinder() = default;

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOP
ResidueTypeFinder::get_representative_type( bool const metapatches ) const
{
	ResidueTypeCOPs base_types = get_possible_base_residue_types( false /* include_unpatchable */ );
	ResidueTypeCOPs filtered_types = base_types;
	apply_filters_after_patches( filtered_types, true /* allow_extra_variants */ );
	if ( ! filtered_types.empty() ) {
		return filtered_types[1];
	}

	ResidueTypeCOPs patched_types = get_patched_types( base_types, true );
	if ( !patched_types.empty()  ) {
		filtered_types = patched_types;
		apply_filters_after_patches( filtered_types, true /* allow_extra_variants */ );
		if ( ! filtered_types.empty() ) {
			return filtered_types[1];
		}
	}

	// If there are metapatches to apply and they are to be considered
	// disable_metapatches() would disable consideration of metapatches (relevant for speed)
	if ( metapatches && !no_metapatches() ) {
		base_types.append( patched_types ); // Because we want to metapatch all of them

		ResidueTypeCOPs metapatched_types = get_singly_metapatched_types( base_types, true );
		while ( ! metapatched_types.empty() ) {
			filtered_types = metapatched_types;
			apply_filters_after_patches( filtered_types, true /* allow_extra_variants */ );
			if ( ! filtered_types.empty() ) {
				return filtered_types[1];
			}

			metapatched_types = get_singly_metapatched_types( metapatched_types, true ); // Go for N+1 metapatches
		}
	}

	// No base or (meta)patched types fit -- try the unpatchable types.
	filtered_types = get_possible_unpatchable_residue_types();
	apply_filters_after_patches( filtered_types, true /* allow_extra_variants */ );
	if ( ! filtered_types.empty() ) {
		return filtered_types[1];
	}

	// No residue types left after patching and filtering
	return nullptr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::get_all_possible_residue_types( bool const allow_extra_variants /* = false */ ) const
{
	// Get all possible basic residues that might match
	ResidueTypeCOPs rsd_types = get_possible_base_residue_types( false /* include_unpatchable*/ );

	TR.Debug << "Found " << rsd_types.size() << " base ResidueType(s) that might match." << std::endl;

	rsd_types.append( get_patched_types( rsd_types ) );

	// If metapatches to be considered
	// disable_metapatches() would disable consideration of metapatches (relevant for speed)
	if ( !no_metapatches() ) {
		ResidueTypeCOPs metapatched_types = get_singly_metapatched_types( rsd_types );
		while ( ! metapatched_types.empty() ) {
			rsd_types.append( metapatched_types );

			metapatched_types = get_singly_metapatched_types( metapatched_types ); // Go for N+1 metapatches
		}
	}

	// add in any unpatchable residues.
	rsd_types.append( get_possible_unpatchable_residue_types() );

	TR.Debug << "Patched up to " << rsd_types.size() << " ResidueType(s)." << std::endl;

	// Filter for rsd_types that strictly obey requirements
	apply_filters_after_patches( rsd_types, allow_extra_variants );

	rsd_types = apply_preferences_and_discouragements( rsd_types );

	TR.Debug << "Keeping up to " << rsd_types.size() << " ResidueType(s) after filtering." << std::endl;

	return rsd_types;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOP
ResidueTypeFinder::get_best_match_residue_type_for_atom_names( utility::vector1< std::string > const & atom_names )
{
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

	// If only one rsd_type to choose from, return it
	if ( n_types == 1 ) {
		return rsd_types[ 1 ];
	}

	// Otherwise, find the rsd_type best match
	if ( TR.Debug.visible() ) {
		TR.Debug << "Finding best match from among " << n_types << " ResidueTypes." << std::endl;
		for ( uint i( 1 ); i <= n_types; ++i ) {
			TR.Trace << ' ' << rsd_types[ i ]->name() << std::endl;
		}
	}
	return find_best_match( rsd_types, atom_names, ignore_atom_named_H_ );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::get_possible_base_residue_types(
	bool const include_unpatchable /*=true*/, bool const apply_all_filters /*=false*/ ) const
{
	//If a base type has already been specified,
	//there's no need to bother with a lot of other rigamarole.
	if ( base_type_ ) {
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
	apply_basic_filters( rsd_types );
	if ( apply_all_filters ) {
		apply_filters_after_patches( rsd_types );
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
	// No need to filter if there are no unpatchable residue types
	if ( rsd_types.size() == 0 ) { return rsd_types; }
	apply_basic_filters( rsd_types );
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
void
ResidueTypeFinder::apply_basic_filters( ResidueTypeCOPs & rsd_types ) const
{
	if ( rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes were passed for applying basic filters." << std::endl;
		return;
	}
	filter_by_aa( rsd_types );
	filter_by_name1( rsd_types );
	filter_by_name3( rsd_types, true /* keep_if_base_type_generates_name3 */ );
	filter_by_residue_type_base_name( rsd_types );
	filter_by_interchangeability_group( rsd_types, true /* keep_if_base_type_generates_interchangeability_group */ );
	filter_by_base_property( rsd_types );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::apply_filters_after_patches( ResidueTypeCOPs & rsd_types,
	bool const allow_extra_variants  /* = false */ ) const
{
	if ( rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes were passed for applying filters after patching." << std::endl;
		return;
	}
	filter_by_name3( rsd_types, false /* keep_if_base_type_generates_name3 */ );
	filter_by_interchangeability_group( rsd_types, false /* keep_if_base_type_generates_interchangeability */ );
	filter_disallow_variants( rsd_types );
	filter_all_variants_matched( rsd_types, allow_extra_variants );
	filter_all_properties( rsd_types );
	filter_disallow_properties( rsd_types );
	filter_all_patch_names( rsd_types );
	filter_connections( rsd_types );
	filter_special_cases( rsd_types );
}

////////////////////////////////////////////////////////
ResidueTypeCOPs
ResidueTypeFinder::apply_preferences_and_discouragements( ResidueTypeCOPs const & rsd_types ) const {
	if ( rsd_types.empty() ) return rsd_types;

	if ( preferred_properties_.empty() && discouraged_properties_.empty() && preferred_connects_.empty() && discouraged_connects_.empty() && ! no_CCD_on_name3_match_ ) {
		return rsd_types; // nothing to do
	}

	if ( rsd_types.size() == 1 ) return rsd_types; // If there's only one possibility, we're going to be using it.

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
			TR.Debug << "Discouraging " << discouraged_properties_.size() << " properties, " <<
				"going from " << current_type_list.size() << " types to " <<
				new_type_list.size() << " types." << std::endl;
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

	if ( ! discouraged_connects_.empty() ) {
		utility::vector1< core::Size > connect_counts;
		for ( ResidueTypeCOP const & type: current_type_list ) {
			core::Size count = 0;
			for ( std::string const & connect_point: discouraged_connects_ ) {
				if ( connect_point == "UPPER" && type->upper_connect_id() != 0 ) {
					++count;
				} else if ( connect_point == "LOWER" && type->lower_connect_id() != 0 ) {
					++count;
				} else if ( type->has(connect_point) && type->atom_forms_residue_connection( type->atom_index(connect_point) ) ) {
					++count;
				}
			}
			connect_counts.push_back(count);
		}
		debug_assert( ! connect_counts.empty() );
		core::Size min_count = *std::min_element( connect_counts.begin(), connect_counts.end() );
		new_type_list.clear();
		for ( core::Size ii(1); ii <= current_type_list.size(); ++ii ) {
			if ( connect_counts[ ii ] == min_count ) {
				new_type_list.push_back( current_type_list[ ii ] );
			}
		}
		if ( TR.Debug.visible() ) {
			TR.Debug << "Discouraging " << discouraged_connects_.size() << " connection points, " <<
				"going from " << current_type_list.size() << " types to " <<
				new_type_list.size() << " types." << std::endl;
			TR.Debug<< "Discouraged connections: " << discouraged_connects_ << std::endl;
			TR.Debug<< "Going from ";
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
			TR.Debug << "Encouraging " << preferred_properties_.size() << " properties, " <<
				"going from " << current_type_list.size() << " types to " <<
				new_type_list.size() << " types." << std::endl;
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

	if ( ! preferred_connects_.empty() ) {
		utility::vector1< core::Size > connect_counts;
		for ( ResidueTypeCOP const & type: current_type_list ) {
			core::Size count = 0;
			for ( std::string const & connect_point: preferred_connects_ ) {
				if ( connect_point == "UPPER" && type->upper_connect_id() != 0 ) {
					++count;
				} else if ( connect_point == "LOWER" && type->lower_connect_id() != 0 ) {
					++count;
				} else if ( type->has(connect_point) && type->atom_forms_residue_connection( type->atom_index(connect_point) ) ) {
					++count;
				}
			}
			connect_counts.push_back(count);
		}
		debug_assert( ! connect_counts.empty() );
		core::Size max_count = *std::max_element( connect_counts.begin(), connect_counts.end() );
		new_type_list.clear();
		for ( core::Size ii(1); ii <= current_type_list.size(); ++ii ) {
			if ( connect_counts[ ii ] == max_count ) {
				new_type_list.push_back( current_type_list[ ii ] );
			}
		}
		if ( TR.Debug.visible() ) {
			TR.Debug << "Encouraging " << preferred_connects_.size() << " connection points, " <<
				"going from " << current_type_list.size() << " types to " <<
				new_type_list.size() << " types." << std::endl;
			TR.Debug<< "Encouraged connections: " << preferred_connects_ << std::endl;
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


/// @details This does not call the filters on the returned types -- that's for the caller to handle.
utility::vector1< ResidueTypeCOP >
ResidueTypeFinder::get_patched_types( utility::vector1< ResidueTypeCOP > const & rsd_types,
	bool const get_first_totally_ok_residue_type /* false */ ) const
{
	utility::vector1< ResidueTypeCOP > patched_types; // The return value. Always return this for NRVO
	if ( rsd_types.empty() ) { return patched_types; }

	for ( PatchCOP patch: residue_type_set_.patches() ) {
		// Make a new vector to iterate over, consisting of the base types, as well as types patched in previous rounds.
		// Make a new vector to avoid modifying the vector we're iterating over
		utility::vector1< ResidueTypeCOP > to_patch( rsd_types );
		to_patch.append( patched_types );
		for ( ResidueTypeCOP const & rsd_type: to_patch ) {
			if ( TR.Trace.visible() ) {
				TR.Trace << "Checking if patch " << patch->name() << " applies to residue " << rsd_type->name() << std::endl;
			}

			// absolute no-no's.
			if ( !patch->applies_to( *rsd_type ) ) { continue; }
			if ( has_disallowed_variant( patch ) ) { continue; }
			if ( deletes_any_property( patch, rsd_type ) ) { continue; }
			if ( deletes_any_variant( patch, rsd_type ) ) { continue; }
			if ( changes_to_wrong_aa( patch, rsd_type ) ) { continue; }
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
					patched_types.push_back( rsd_type_new );

					// If we've found a reasonable type, we can stop early.
					if ( get_first_totally_ok_residue_type ) {
						ResidueTypeCOPs filtered{ {rsd_type_new} }; // Just the most recent entry.
						apply_filters_after_patches( filtered, true /*allow_extra_variants*/ );
						if ( ! filtered.empty() ) {
							patched_types = filtered; // This dance is to enable NRVO
							return patched_types;
						}
					}
				}
			} // if apply_patch
		} // for residue types
	} // for patches

	return patched_types;
}

/// @details This does not call the filters on the returned types -- that's for the caller to handle.
utility::vector1< ResidueTypeCOP >
ResidueTypeFinder::get_singly_metapatched_types( utility::vector1< ResidueTypeCOP > const & rsd_types,
	bool const get_first_totally_ok_residue_type /* false */ ) const
{
	utility::vector1< ResidueTypeCOP > patched_types; // The return value. Always return this for NRVO
	if ( rsd_types.empty() ) { return patched_types; }
	// If metapatches are not to be considered, return an empty list
	// no_metapatches() == true means do not consider metapatches
	if ( no_metapatches() ) { return patched_types; }

	for ( MetapatchCOP const & metapatch: residue_type_set_.metapatches() ) {
		for ( ResidueTypeCOP const & rsd_type: rsd_types ) {
			// Don't attempt to run through the atoms if the metapatch not applicable to the residue type.
			if ( ! metapatch->applies_to( *rsd_type ) ) { continue; }

			for ( std::string const & atom_name: metapatch->atoms( *rsd_type ) ) {
				PatchCOP patch = metapatch->get_one_patch( atom_name );

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
					fixes_interchangeability_group( patch, rsd_type ) // ||
					// fixes_connects( patch, rsd_type ) // Currently an issue, as the connect metapatch plays havoc
				);

				if ( apply_patch ) {
					// following just gets the right name of the patched residue
					std::string const & patched_name( patch->patched_name( *rsd_type ) );
					ResidueTypeCOP rsd_type_new( residue_type_set_.name_mapOP( patched_name ) );
					if ( rsd_type_new ) {
						patched_types.push_back( rsd_type_new );

						// If we've found a reasonable type, we can stop early.
						if ( get_first_totally_ok_residue_type ) {
							ResidueTypeCOPs filtered{ {rsd_type_new} }; // Just the most recent entry.
							apply_filters_after_patches( filtered, true /*allow_extra_variants*/ );
							if ( ! filtered.empty() ) {
								patched_types = filtered; // This dance is to enable NRVO
								return patched_types;
							}
						}
					}
				} // if apply_patch
			} // for atoms
		} // for residue types
	} // for metapatches

	return patched_types;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_by_name1( ResidueTypeCOPs & rsd_types  ) const
{
	if ( name1_ == '?' || rsd_types.empty() ) { return; }

	rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
		[&](ResidueTypeCOP const & rsd_type) {
			return rsd_type->name1() != name1_;
		}),
		rsd_types.end()
	);

	if ( rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering by one-letter code: '" << name1_ << "'" << std::endl;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_by_name3( ResidueTypeCOPs & rsd_types, bool const keep_if_base_type_generates_name3  ) const
{
	if ( name3_.empty() || rsd_types.empty() ) { return; }

	std::string const & name3_stripped = utility::strip( name3_ );
	rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
		// For specialty amino acids, add them to the name three maps both with their PDB strings and
		// with their specialty string -- the first three letters of the residue name.
		// E.g., CYD will appear in both lists for name3_map_[ "CYS" ] and name3_map_[ "CYD" ]
		[&](ResidueTypeCOP const & rsd_type) {
			std::string const & name3 = rsd_type->name3();
			return name3 != name3_ && utility::strip( name3 ) != name3_stripped &&
			!( keep_if_base_type_generates_name3 && residue_type_set_.generates_patched_residue_type_with_name3( name3, name3_ ) );
		}),
		rsd_types.end()
	);

	if ( rsd_types.empty() ) {
		// RM: This seems to occur frequently in the integration tests.
		// We may want to examine why we're doing so much filtering by name3 on sets that don't contain what we want.
		TR.Debug << "No ResidueTypes remain after filtering by three-letter code: '" << name3_ << "'" << std::endl;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_by_interchangeability_group( ResidueTypeCOPs & rsd_types, bool const keep_if_base_type_generates_interchangeability_group  ) const
{
	if ( interchangeability_group_.empty() || rsd_types.empty() ) { return; }

	rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
		[&](ResidueTypeCOP const & rsd_type) {
			return !(
			rsd_type->interchangeability_group() == interchangeability_group_ ||
			(keep_if_base_type_generates_interchangeability_group && residue_type_set_.generates_patched_residue_type_with_interchangeability_group( rsd_type->name(), interchangeability_group_ ) )
			);
		}),
		rsd_types.end()
	);

	if ( rsd_types.empty() ) {
		// RM: This seems to occur frequently in the integration tests.
		// We may want to examine why we're doing so much filtering on sets that don't contain what we want.
		TR.Debug << "No ResidueTypes remain after filtering by interchangeability group: '" <<
			interchangeability_group_ << "'" << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_by_aa( ResidueTypeCOPs & rsd_types ) const
{
	if ( aa_ == aa_none || rsd_types.empty() ) { return; }

	rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
		[&](ResidueTypeCOP const & rsd_type) {
			return rsd_type->aa() != aa_;
		}),
		rsd_types.end()
	);

	if ( rsd_types.empty() ) {
		// RM: This seems to occur frequently in the integration tests.
		// We may want to examine why we're doing so much filtering by aa_type on sets that don't contain what we want.
		TR.Debug << "No ResidueTypes remain after filtering by amino acid designation: '" << aa_ << "'" << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_by_residue_type_base_name( ResidueTypeCOPs & rsd_types ) const
{
	if ( residue_type_base_name_.size() == 0 || rsd_types.empty() ) { return; }

	rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
		[&](ResidueTypeCOP const & rsd_type) {
			return rsd_type->base_name() != residue_type_base_name_ && residue_type_base_name( *rsd_type ) != residue_type_base_name_;
		}),
		rsd_types.end()
	);

	if ( rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering by ResidueType base name: '" <<
			residue_type_base_name_ << "'" << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_by_base_property( ResidueTypeCOPs & rsd_types ) const
{
	if ( base_property_ == NO_PROPERTY || rsd_types.empty() ) { return; }

	rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
		[&](ResidueTypeCOP const & rsd_type) {
			return !rsd_type->properties().has_property( base_property_ );
		}),
		rsd_types.end()
	);

	if ( rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering by ResidueType base property: '" <<
			base_property_ << "'" << std::endl;
	}
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
		MutableResidueTypeCOP new_type( patch->apply( *rsd_type, false /*instantiate*/ ) );
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
		MutableResidueTypeCOP new_type( patch->apply( *rsd_type, false /*instantiate*/ ) );
		if ( new_type && new_type->interchangeability_group() == interchangeability_group_ ) {
			return true;
		}
		// @mlnance: why is this return true here? seems like the above if-statement doesn't matter then
		return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
adds_connects_helper( PatchCOP patch, ResidueTypeCOP rsd_type, std::string const & atom ) {
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
	return false;
}

bool
ResidueTypeFinder::fixes_connects( PatchCOP patch, ResidueTypeCOP rsd_type ) const {

	for ( std::string const & atom: connect_atoms_ ) {
		if ( adds_connects_helper( patch, rsd_type, atom ) ) { return true; }
	}

	for ( std::string const & atom: preferred_connects_ ) {
		if ( atom == "LOWER" && rsd_type->lower_connect_id() == 0 ) {
			if ( patch->changes_connections_on( *rsd_type, atom ) ) { return true; }
		} else if ( atom == "UPPER" && rsd_type->upper_connect_id() == 0 ) {
			if ( patch->changes_connections_on( *rsd_type, atom ) ) { return true; }
		} else {
			if ( adds_connects_helper( patch, rsd_type, atom ) ) { return true; }
		}
	}

	for ( std::string const & atom: discouraged_connects_ ) {
		if ( atom == "LOWER" && rsd_type->lower_connect_id() != 0 ) {
			if ( patch->changes_connections_on( *rsd_type, atom ) ) { return true; }
		} else if ( atom == "UPPER" && rsd_type->upper_connect_id() != 0 ) {
			if ( patch->changes_connections_on( *rsd_type, atom ) ) { return true; }
		} else if ( rsd_type->has( atom ) && rsd_type->atom_forms_residue_connection( rsd_type->atom_index(atom) ) ) {
			if ( patch->changes_connections_on( *rsd_type, atom ) ) { return true; }
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
/// @brief Returns true if this patch deletes any of the properties that the ResidueTypeFinder is seeking.
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
		if ( atom_names_soft_.has_value( ObjexxFCL::strip_whitespace( new_atom_name ) ) ) {
			// @mlnance a verbose debug here is helpful for carbohydrate debugging
			TR.Debug << "Patch '" << patch->name() <<
				"' is relevant because this patch has atom name " << new_atom_name <<
				" which is an atom name found in the PDB file" << std::endl;
			return true;
		}
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This may not work with custom variants. Crap.
void
ResidueTypeFinder::filter_all_variants_matched( ResidueTypeCOPs & rsd_types, bool const allow_extra_variants /* = false */ ) const
{
	if ( rsd_types.empty() ) { return; }

	filter_for_all_variant_sets( rsd_types );
	if ( !allow_extra_variants ) {
		filter_variant_sets_have_all_candidate_variants( rsd_types );
	}

	if ( TR.Debug.visible() && rsd_types.empty() ) {
		// RM: This seems to occur frequently in the integration tests.
		// We may want to examine why we're doing so much filtering by matched variants on sets that don't contain what we want.
		TR.Debug << "No ResidueTypes remain after filtering for matched variants." << std::endl;
		//TODO: Print what those matched variants should have been.
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_for_all_variant_sets( ResidueTypeCOPs & rsd_types ) const
{
	if ( rsd_types.empty() ) { return; }

	rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
		[&](ResidueTypeCOP const & rsd_type) {
			for ( auto const & variant_set : variants_in_sets_ ) {
				if ( variant_set.size() == 1 && variant_exceptions_.has_value( variant_set[ 1 ] ) ) continue;

				bool at_least_one_variant_matched( false );
				for ( auto const & variant: variant_set ) {
					if ( rsd_type->properties().is_variant_type( variant ) ) {
						at_least_one_variant_matched = true;
						break;
					}
				} // loop inside one variant set

				if ( !at_least_one_variant_matched ) {
					return true;
				}
			}

			for ( auto const & custom_variant : custom_variants_ ) {
				if ( !rsd_type->properties().is_variant_type( custom_variant ) ) {
					return true;
				}
			}

			return false;
		}),
		rsd_types.end()
	);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_variant_sets_have_all_candidate_variants( ResidueTypeCOPs & rsd_types ) const
{
	if ( rsd_types.empty() ) { return; }

	rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
		[&](ResidueTypeCOP const & rsd_type) {

			for ( auto const & variant_type : rsd_type->variant_type_enums() ) {
				if ( variant_exceptions_.has_value( variant_type ) ) continue;

				bool at_least_one_variant_matched( false );
				for ( auto const & variant_set : variants_in_sets_ ) {
					if ( variant_set.contains( variant_type ) ) {
						at_least_one_variant_matched = true;
						break;
					}
				}

				if ( !at_least_one_variant_matched ) {
					return true;
				}
			}

			//Now, we need to check on-the-fly types.  This and only this requires string parsing:
			for ( std::string const & custom_variant : rsd_type->custom_variant_types() ) {
				if ( !custom_variants_.has_value( custom_variant ) ) {
					return true;
				}
			}
			return false;

		}),
		rsd_types.end()
	);

}

////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_all_patch_names( ResidueTypeCOPs & rsd_types )  const
{
	if ( patch_names_.empty() || rsd_types.empty() ) { return; }

	rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
		[&](ResidueTypeCOP const & rsd_type) {
			auto const & all_patches_names = residue_type_all_patches_name( *rsd_type );
			for ( std::string const & patch_name: patch_names_ ) {
				if ( all_patches_names.find( patch_name ) == std::string::npos ) {
					return true;
				}
			}
			return false;
		}),
		rsd_types.end()
	);

	if ( TR.Debug.visible() && rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering for the presence of all patch names." << std::endl;
		TR.Debug << "    Patches sought: ";
		for ( std::string const & patch_name: patch_names_ ) {
			TR.Debug << patch_name << "   ";
		}
		TR.Debug << std::endl;
	}
}
////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_all_properties( ResidueTypeCOPs & rsd_types )  const
{
	if ( properties_.empty() || rsd_types.empty() ) { return; }

	rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
		[&](ResidueTypeCOP const & rsd_type) {
			for ( auto const & property: properties_ ) {
				if ( !rsd_type->has_property( property ) ) {
					return true;
				}
			}
			return false;
		}),
		rsd_types.end()
	);

	if ( TR.Debug.visible() && rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering for the presence of all properties." << std::endl;
		TR.Debug << "    Properties sought: ";
		for ( auto const & property: properties_ ) {
			TR.Debug << ResidueProperties::get_string_from_property( property ) << "   ";
		}
		TR.Debug << std::endl;
	}
}
////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_disallow_variants( ResidueTypeCOPs & rsd_types )  const
{
	if ( disallow_variants_.empty() || rsd_types.empty() ) { return; }

	rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
		[&](ResidueTypeCOP const & rsd_type) {
			for ( auto const & variant: disallow_variants_ ) {
				if ( rsd_type->has_variant_type( variant ) ) {
					return true;
				}
			}
			return false;
		}),
		rsd_types.end()
	);

	if ( TR.Debug.visible() && rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering against the presence of variants." << std::endl;
		TR.Debug << "    Variants prohibited: ";
		for ( auto const & variant: disallow_variants_ ) {
			TR.Debug << ResidueProperties::get_string_from_variant( variant ) << "   ";
		}
		TR.Debug << std::endl;
	}
}
////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_disallow_properties( ResidueTypeCOPs & rsd_types )  const
{
	if ( disallow_properties_.empty() || rsd_types.empty() ) { return; }

	rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
		[&](ResidueTypeCOP const & rsd_type) {
			for ( auto const & property: disallow_properties_ ) {
				if ( rsd_type->has_property( property ) ) {
					return true;
				}
			}
			return false;
		}),
		rsd_types.end()
	);

	if ( TR.Debug.visible() && rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering against the presence of properties." << std::endl;
		TR.Debug << "    Properties prohibited: ";
		for ( auto const & property: disallow_properties_ )  {
			TR.Debug << ResidueProperties::get_string_from_property( property ) << "   ";
		}
		TR.Debug << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_connections( ResidueTypeCOPs & rsd_types )  const
{
	if ( connect_atoms_.empty() || rsd_types.empty() ) { return; }

	rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
		[&](ResidueTypeCOP const & rsd_type) {
			for ( std::string const & atom: connect_atoms_ ) {
				if ( !rsd_type->has( atom ) ||
						rsd_type->residue_connections_for_atom( rsd_type->atom_index(atom) ).empty() ) {
					return true;
				}
			}
			return false;
		}),
		rsd_types.end()
	);

	if ( TR.Debug.visible() && rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering against the presence of connections." << std::endl;
		TR.Debug << "    Connection atoms required: ";
		for ( std::string const & atom: connect_atoms_ )  {
			TR.Debug << atom << "   ";
		}
		TR.Debug << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////
void
ResidueTypeFinder::filter_special_cases( ResidueTypeCOPs & rsd_types )  const
{
	if ( rsd_types.empty() ) { return; }

	bool const actually_check_nucleic_acid_virtual_phosphates =
		check_nucleic_acid_virtual_phosphates_ && !atom_names_soft_.has_value( "P" );

	if ( actually_check_nucleic_acid_virtual_phosphates ) {
		rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
			[&](ResidueTypeCOP const & rsd_type) {
				if ( rsd_type->is_DNA() && !rsd_type->has_variant_type( VIRTUAL_DNA_PHOSPHATE ) ) { return true; }

				// triphosphate termination does not require virtual phosphates!
				if ( rsd_type->is_RNA() && !rsd_type->has_variant_type( VIRTUAL_PHOSPHATE )
						&& !rsd_type->has_variant_type( FIVE_PRIME_PACKABLE_TRIPHOSPHATE ) ) { return true; }

				return false;
			}),
			rsd_types.end()
		);
	}

	bool const is_water_and_user_specified = ( aa_ == AA::aa_h2o || name1_ == 'w' ) &&
		option[ OptionKeys::in::water_type_if_unspecified ].user();
	std::string const & user_water_name = option[ OptionKeys::in::water_type_if_unspecified ]();

	if ( is_water_and_user_specified ) {
		rsd_types.erase( std::remove_if( rsd_types.begin(), rsd_types.end(),
			[&](ResidueTypeCOP const & rsd_type) {
				return rsd_type->name() != user_water_name;
			}),
			rsd_types.end()
		);
	}

	if ( TR.Debug.visible() && rsd_types.empty() ) {
		TR.Debug << "No ResidueTypes remain after filtering for special cases:" << std::endl;
		if ( actually_check_nucleic_acid_virtual_phosphates ) {
			TR.Debug << "    * Nucleic acid has virtual phosphates." << std::endl;
		}
		if ( is_water_and_user_specified ) {
			TR << "    * Water residue: No match for user specified name: \"" << user_water_name << "\"" << std::endl;
		}
	}
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
			utility_exit_with_message( "Not currently handling custom variants within variant_exceptions!" );
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
