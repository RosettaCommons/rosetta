// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/adduct_util.cc
/// @author Jim Havranek


// Unit header
#include <core/chemical/adduct_util.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>

// Project headers
#include <core/chemical/VariantType.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/conversions.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/options/StringVectorOption.hh>

// Boost Headers
#include <boost/foreach.hpp>


static THREAD_LOCAL basic::Tracer TR( "core.chemical.adduct_util" );

namespace core {
namespace chemical {

/// @brief Convert input string to map of adducts->max usage
std::map< std::string, int >
parse_adduct_string(
	utility::options::StringVectorOption & add_vec
)
{
	std::map< std::string, int > add_map;
	// walk via an index, plucking a second string each time if
	// it is an integer
	Size index( 1 );
	Size over_index( add_vec.size() + 1 );
	while ( index != over_index ) {
		// Always store as lower for ease of comparison
		std::string this_adduct( add_vec[ index ]  );
		ObjexxFCL::lowercase( this_adduct );
		index++;
		int max_uses_for_this_adduct( 9999 );
		if ( index != over_index && ObjexxFCL::is_int( add_vec[ index ] ) ) {
			max_uses_for_this_adduct = ObjexxFCL::int_of( add_vec[ index ] );
			index++;
		}
		// Add to the map
		add_map[ this_adduct ] = max_uses_for_this_adduct;
	}

	// Debug - remove later
	// for( std::map< std::string, int >::iterator iter = add_map.begin() ;
	//    iter != add_map.end() ; iter++ ) {
	//  TR << "Adduct map " << iter->first << "\t" << iter->second << std::endl;
	// }


	return add_map;
}


/// @brief Make sure any adducts requested actually exist
void
error_check_requested_adducts( std::map< std::string, int > const & add_map,
	ResidueTypeCOPs const & rsd_types ) {

	for ( std::map< std::string, int >::const_iterator this_add = add_map.begin() ;
			this_add != add_map.end() ; ++this_add ) {
		bool not_found( true );

		for ( ResidueTypeCOPs::const_iterator this_rsd = rsd_types.begin(),
				end_rsd = rsd_types.end(); this_rsd != end_rsd ; ++this_rsd ) {
			// shortcircuit if we've already found an instance of the adduct
			ResidueType const & rsd( **this_rsd );
			if ( not_found == false ) break;

			for ( utility::vector1< Adduct >::const_iterator rsd_add = rsd.defined_adducts().begin(),
					end_rsd_add = rsd.defined_adducts().end() ;
					rsd_add != end_rsd_add ; ++rsd_add ) {
				std::string check_name( rsd_add->adduct_name() );
				// compare case-insensitively for convenience
				if ( ObjexxFCL::equali( this_add->first, check_name ) ) {
					not_found = false;
					break;
				}
			}
		}

		if ( not_found ) {
			utility_exit_with_message( "Requested undefined adduct: " + this_add->first + '\n' );
		}
	} // Done with adduct name error-checking
}

/// @brief Apply adducts to residue using a boolean mask
ResidueTypeOP apply_adducts_to_residue( ResidueType const & rsd,
	utility::vector1< bool > & add_mask
)
{
	using numeric::conversions::radians;

	// Use the patching machinery to apply the adducts
	PatchCase temp_patch_case;

	// Starting name
	std::string new_rsd_name( rsd.name() );

	// Throw in all the applicable adducts
	utility::vector1< bool >::iterator mask_iter = add_mask.begin();
	for ( utility::vector1< Adduct >::const_iterator add_iter = rsd.defined_adducts().begin() ,
			end_add_iter = rsd.defined_adducts().end() ; add_iter != end_add_iter ;
			++add_iter, ++mask_iter )  {

		// Skip adducts if dictated by the mask
		if ( !(*mask_iter) ) continue;

		// Add the adduct and it's information
#if defined(WIN32) && !defined(WIN_PYROSETTA)
		PatchOperationOP poop1 = new AddAtomWIN32( add_iter->atom_name(), add_iter->atom_type_name(), add_iter->mm_atom_type_name(), add_iter->atom_charge() );
#else
		PatchOperationOP poop1( new AddAtom( add_iter->atom_name(), add_iter->atom_type_name(), add_iter->mm_atom_type_name(), add_iter->atom_charge() ) );
#endif

		temp_patch_case.add_operation( poop1 );
		PatchOperationOP poop2( new AddBond( add_iter->atom_name(), add_iter->stub_atom1() ) );
		temp_patch_case.add_operation( poop2 );
		PatchOperationOP poop3( new SetICoor( add_iter->atom_name(), radians( add_iter->phi() ),
			radians( add_iter->theta() ), add_iter->d(), add_iter->stub_atom1(), add_iter->stub_atom2(),
			add_iter->stub_atom3() ) );
		TR.Debug << "Making a patch op for residue " << rsd.name() << " adduct " << add_iter->adduct_name() <<
			" phi raw " << add_iter->phi() << " theta raw " << add_iter->theta() << " d raw " << add_iter->d() << std::endl;
		TR.Debug << "Using stub atoms " << add_iter->stub_atom1() << " , " << add_iter->stub_atom2() << " , " <<
			add_iter->stub_atom3() << std::endl;
		temp_patch_case.add_operation( poop3 );

		// Append the adduct name
		new_rsd_name = new_rsd_name + "_adduct:" + add_iter->atom_name();
	}

	// Apply the Patch operations
	ResidueTypeOP new_rsd_type( temp_patch_case.apply( rsd ) );

	// Set the full name
	TR << "Setting new rsd name to " << new_rsd_name << std::endl;
	new_rsd_type->name( new_rsd_name );

	// Set the variant type
	new_rsd_type->add_variant_type( ADDUCT_VARIANT );

	// Set as an adduct-modified type, which is helpful to know
	// in rotamer-building, as adduct variants are generally allowed,
	// while non-adduct variants are in general disallowed
	new_rsd_type->set_adduct_flag( true );

	new_rsd_type->finalize();
	return new_rsd_type;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// adducts
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// @details Generation of new residue types augmented by adduct atoms
/// @note    This is almost perfectly consistent with on_the_fly residue type sets.
///            Why not just make adducts a patch? Then it will all work together nicely.
///            Just need to have the -adduct string (add_map) turned into a PatchSelector
///                 -- rhiju
void
place_adducts( ResidueTypeSet & rsd_type_set )
{
	// First parse the command line for requested adducts
	utility::options::StringVectorOption & add_set
		= basic::options::option[ basic::options::OptionKeys::packing::adducts ];

	// No adducts, skip out
	if ( add_set.size() == 0 ) return;

	// Convert to a map that takes a string descriptor of an adduct and
	// gives the max number of adducts of that class to apply, so
	// a command line option of -adducts <adduct_type> 2 for a type that has
	// 5 entries in the rsd param file will create all combinations with up to
	// 2 total adducts.

	AdductMap add_map = parse_adduct_string( add_set );

	// Error check each requested adduct from the command line, and
	// complain if there are no examples in any residues.  This function
	// will not return if
	error_check_requested_adducts( add_map, rsd_type_set.base_residue_types() );

	// Set up a starting point map where the int value is the number
	// of adducts of a given type placed
	AdductMap blank_map( add_map );
	for ( AdductMap::iterator add_iter = blank_map.begin(),
			end_iter = blank_map.end() ;
			add_iter != end_iter ; ++add_iter ) {
		add_iter->second = 0;
	}

	// Process the residues in turn
	ResidueTypeCOPs residue_types = ResidueTypeFinder( rsd_type_set ).base_property( DNA ).variant_exceptions( utility::tools::make_vector1( LOWER_TERMINUS_VARIANT, UPPER_TERMINUS_VARIANT, ADDUCT_VARIANT ) ).get_all_possible_residue_types();

	for ( ResidueTypeCOPs::const_iterator iter= residue_types.begin(), iter_end = residue_types.end();
			iter != iter_end; ++iter ) {
		ResidueType const & rsd( **iter );
		if ( !rsd.finalized() ) continue;
		AdductMap count_map( blank_map );
		utility::vector1< bool > add_mask( rsd.defined_adducts().size(), false  );
		create_adduct_combinations( rsd_type_set, rsd, add_map, count_map, add_mask, rsd.defined_adducts().begin() );
	}

	// utility_exit_with_message( "Debug stop point \n" );

}

/// @brief Create correct combinations of adducts for a residue type
void
create_adduct_combinations(
	ResidueTypeSet & rsd_type_set,
	ResidueType const & rsd,
	AdductMap ref_map,
	AdductMap count_map,
	utility::vector1< bool > add_mask,
	utility::vector1< Adduct >::const_iterator work_iter
)
{

	if ( work_iter == rsd.defined_adducts().end() ) {
		// Skip the 'no adduct' case - that has already been
		// made when reading in files
		if ( std::find( add_mask.begin(), add_mask.end(), true ) == add_mask.end() ) {
			return;
		}
		// Make this combo and return;
		//  std::cout << "Making an adduct" << std::endl;

		utility::vector1< Adduct >::const_iterator add_iter = rsd.defined_adducts().begin() ;
		BOOST_FOREACH ( bool make, add_mask ) {
			std::cout << "Adduct " << add_iter->adduct_name() << " make is " << make << std::endl;
			++add_iter;
		}

		// Farm this out to a helper function
		rsd_type_set.add_custom_residue_type( apply_adducts_to_residue( rsd, add_mask ) );

		return;
	}

	// Traverse the 'make' branch for this adduct if:
	// 1. The adduct is in the map of requested adducts
	// 2. we haven't exceeded the count limit for this adduct
	AdductMap::iterator test_iter =
		ref_map.find( work_iter->adduct_name() );

	if ( test_iter != ref_map.end() &&
			count_map[ test_iter->first ] < ref_map[ test_iter->first ]   ) {
		AdductMap new_count_map( count_map );
		new_count_map[ work_iter->adduct_name() ]++;
		utility::vector1< bool > new_add_mask( add_mask );
		// This following line may not work if the Adducts are no longer
		// stored in a vector
		new_add_mask[ work_iter - rsd.defined_adducts().begin() + 1 ] = true;
		create_adduct_combinations( rsd_type_set, rsd, ref_map, new_count_map, new_add_mask, work_iter+1 );
	}

	// Always traverse the 'do not make' for this adduct
	// The count is not incremented, and the mask is left at the default (false)
	AdductMap new_count_map( count_map );
	utility::vector1< bool > new_add_mask( add_mask );
	create_adduct_combinations( rsd_type_set, rsd, ref_map, new_count_map, new_add_mask, work_iter+1 );

}




} // chemical
} // core
