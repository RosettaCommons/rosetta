// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/import_pose/options/RNA_MinimizerOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/import_pose/options/RNA_MinimizerOptions.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.import_pose.options.RNA_MinimizerOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace utility::tag;

namespace core {
namespace import_pose {
namespace options {

//Constructor
RNA_MinimizerOptions::RNA_MinimizerOptions():
	minimize_rounds_( 2 ),
	deriv_check_( false ),
	skip_o2prime_trials_( false ),
	vary_bond_geometry_( false ),
	minimizer_use_coordinate_constraints_( false ),
	min_type_( "lbfgs_armijo_nonmonotone" ), //Parin S. Jan 12, 2012
	minimize_bps_( false ),
	minimize_all_protein_( false ),
	minimize_protein_sc_( false ),
	protein_packing_( false ),
	protein_pack_all_( false ),
	protein_packing_distance_( 10.0 ),
	use_nblist_( true ),
	nblist_auto_update_( true )
{}

/// @brief clone the options
RNA_MinimizerOptionsOP
RNA_MinimizerOptions::clone() const
{
	return utility::pointer::make_shared< RNA_MinimizerOptions >( *this );
}

///////////////////////////////////////////////////////////////////
void
RNA_MinimizerOptions::initialize_from_command_line() {
	RNA_BasicOptions::initialize_from_options( basic::options::option );
	initialize_from_options( basic::options::option );
}

void
RNA_MinimizerOptions::initialize_from_options( utility::options::OptionCollection const & opts ) {

	RNA_BasicOptions::initialize_from_options( opts );

	set_minimize_rounds( opts[ OptionKeys::rna::denovo::minimize::minimize_rounds ]() );
	set_vary_bond_geometry( opts[ OptionKeys::rna::vary_geometry ]() || opts[ OptionKeys::relax::cartesian ]() );
	set_deriv_check( opts[ OptionKeys::rna::denovo::minimize::deriv_check ]() );
	set_skip_o2prime_trials( opts[ OptionKeys::rna::denovo::minimize::skip_o2prime_trials]() );

	set_extra_minimize_res( opts[ OptionKeys::rna::denovo::minimize::extra_minimize_res ]() ) ;
	set_extra_minimize_chi_res( opts[ OptionKeys::rna::denovo::minimize::extra_minimize_chi_res ]() ) ;

	set_minimizer_use_coordinate_constraints( opts[ OptionKeys::rna::denovo::minimize::minimizer_use_coordinate_constraints ]() );
	set_min_type( opts[ OptionKeys::rna::denovo::minimize::min_type ]() );
	if ( opts[ OptionKeys::rna::denovo::minimize::skip_coord_constraints]() ) set_minimizer_use_coordinate_constraints( false );
	set_minimize_bps( opts[ OptionKeys::rna::denovo::minimize::minimize_bps ]() ) ;
	set_minimize_all_protein( opts[ OptionKeys::rna::denovo::minimize::minimize_all_protein ]() );
	set_minimize_protein_sc( opts[ OptionKeys::rna::denovo::minimize::minimize_protein_sc ]() );

	set_protein_packing( opts[ OptionKeys::rna::denovo::minimize::protein_packing ]() );
	set_protein_pack_all( opts[ OptionKeys::rna::denovo::minimize::protein_pack_all ]() );
	set_protein_packing_distance( opts[ OptionKeys::rna::denovo::minimize::protein_packing_distance ]() );

	set_use_nblist( opts[ OptionKeys::rna::denovo::minimize::minimize_use_nblist ]() );
	set_nblist_auto_update( opts[ OptionKeys::rna::denovo::minimize::minimize_nblist_auto_update ]() );
	set_max_iter( opts[ OptionKeys::rna::denovo::minimize::minimize_max_iter ]() );

}

void
RNA_MinimizerOptions::initialize_from_tag( utility::tag::TagCOP const & tag ) {

	RNA_BasicOptions::initialize_from_tag( tag );

	if ( tag->hasOption( "minimize_rounds" ) ) {
		set_minimize_rounds( tag->getOption< core::Size >( "minimize_rounds", minimize_rounds_ ) );
	}
	if ( tag->hasOption( "vary_geometry" ) ) {
		set_vary_bond_geometry( tag->getOption< bool >( "vary_geometry", vary_bond_geometry_ ) );
	}
	if ( tag->hasOption( "deriv_check" ) ) {
		set_deriv_check( tag->getOption< bool >( "deriv_check", deriv_check_ ) );
	}
	if ( tag->hasOption( "skip_o2prime_trials" ) ) {
		set_skip_o2prime_trials( tag->getOption< bool >( "skip_o2prime_trials", skip_o2prime_trials_ ) );
	}

	// AMW: oddly this is just an integer vector. that seems bad. shouldn't this be
	// based on resnum-and-chain? AMW TODO
	if ( tag->hasOption( "extra_minimize_res" ) ) {
		std::istringstream ss;
		ss.str( tag->getOption< std::string >( "extra_minimize_res", "" ) );
		utility::vector1< core::Size > vec;
		while ( ss ) {
			core::Size n;
			ss >> n;
			vec.push_back( n );
		}
		set_extra_minimize_res( vec );
	}
	if ( tag->hasOption( "extra_minimize_chi_res" ) ) {
		std::istringstream ss;
		ss.str( tag->getOption< std::string >( "extra_minimize_chi_res", "" ) );
		utility::vector1< core::Size > vec;
		while ( ss ) {
			core::Size n;
			ss >> n;
			vec.push_back( n );
		}
		set_extra_minimize_chi_res( vec );
	}

	if ( tag->hasOption( "minimizer_use_coordinate_constraints" ) ) {
		set_minimizer_use_coordinate_constraints( tag->getOption< bool >( "minimizer_use_coordinate_constraints", minimizer_use_coordinate_constraints_ ) );
	}
	if ( tag->hasOption( "min_type" ) ) {
		set_min_type( tag->getOption< std::string >( "min_type", min_type_ ) );
	}
	if ( tag->hasOption( "skip_coord_constraints" ) ) {
		set_minimizer_use_coordinate_constraints( false );
	}
	if ( tag->hasOption( "minimize_bps" ) ) {
		set_minimize_bps( tag->getOption< bool >( "minimize_bps", minimize_bps_ ) );
	}
	if ( tag->hasOption( "minimize_all_protein" ) ) {
		set_minimize_all_protein( tag->getOption< bool >( "minimize_all_protein", minimize_all_protein_ ) );
	}
	if ( tag->hasOption( "minimize_protein_sc" ) ) {
		set_minimize_protein_sc( tag->getOption< bool >( "minimize_protein_sc", minimize_protein_sc_ ) );
	}

	if ( tag->hasOption( "protein_packing" ) ) {
		set_protein_packing( tag->getOption< bool >( "protein_packing", protein_packing_ ) );
	}
	if ( tag->hasOption( "protein_pack_all" ) ) {
		set_protein_pack_all( tag->getOption< bool >( "protein_pack_all", protein_pack_all_ ) );
	}
	if ( tag->hasOption( "protein_packing_distance" ) ) {
		set_protein_packing_distance( tag->getOption< core::Real >( "protein_packing_distance", protein_packing_distance_ ) );
	}

	if ( tag->hasOption( "minimize_use_nblist" ) ) {
		set_use_nblist( tag->getOption< bool >( "minimize_use_nblist", use_nblist_ ) );
	}
	if ( tag->hasOption( "minimize_nblist_auto_update" ) ) {
		set_nblist_auto_update( tag->getOption< bool >( "minimize_nblist_auto_update", nblist_auto_update_ ) );
	}
	set_max_iter( tag->getOption< core::Size >( "minimize_max_iter", 2000 ) );
}


void
RNA_MinimizerOptions::list_attributes( AttributeList & attlist ) {

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "minimize_rounds", xsct_non_negative_integer, "Number of min cycles", "2" )
		+ XMLSchemaAttribute::attribute_w_default( "vary_geometry", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "deriv_check", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "skip_o2prime_trials", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "extra_minimize_res", xsct_nnegative_int_wsslist, "XRW TODO", "" )
		+ XMLSchemaAttribute::attribute_w_default( "extra_minimize_chi_res", xsct_nnegative_int_wsslist, "XRW TODO", "" )
		+ XMLSchemaAttribute::attribute_w_default( "minimizer_use_coordinate_constraints", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "min_type", xs_string, "XRW TODO", "lbfgs_armijo_nonmonotone" ) // AMW TODO possible min types enum regex?
		+ XMLSchemaAttribute::attribute_w_default( "skip_coord_constraints", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "minimize_bps", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "minimize_all_protein", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "minimize_protein_sc", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "protein_packing", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "protein_pack_all", xsct_rosetta_bool, "XRW TODO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "protein_packing_distance", xsct_real, "XRW TODO", "10.0" )
		+ XMLSchemaAttribute::attribute_w_default( "minimize_use_nblist", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "minimize_nblist_auto_update", xsct_rosetta_bool, "XRW TODO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "minimize_max_iter", xsct_non_negative_integer, "XRW TODO", "2000" );

	RNA_BasicOptions::list_attributes( attlist );
}

void
RNA_MinimizerOptions::list_options_read( utility::options::OptionKeyList & opts ) {
	using namespace basic::options;
	RNA_BasicOptions::list_options_read( opts );
	opts + OptionKeys::rna::denovo::minimize::minimize_rounds
		+ OptionKeys::rna::vary_geometry
		+ OptionKeys::relax::cartesian
		+ OptionKeys::rna::denovo::minimize::deriv_check
		+ OptionKeys::rna::denovo::minimize::skip_o2prime_trials
		+ OptionKeys::rna::denovo::minimize::extra_minimize_res
		+ OptionKeys::rna::denovo::minimize::extra_minimize_chi_res
		+ OptionKeys::rna::denovo::minimize::minimizer_use_coordinate_constraints
		+ OptionKeys::rna::denovo::minimize::min_type
		+ OptionKeys::rna::denovo::minimize::skip_coord_constraints
		+ OptionKeys::rna::denovo::minimize::minimize_bps
		+ OptionKeys::rna::denovo::minimize::minimize_all_protein
		+ OptionKeys::rna::denovo::minimize::minimize_protein_sc
		+ OptionKeys::rna::denovo::minimize::protein_packing
		+ OptionKeys::rna::denovo::minimize::protein_pack_all
		+ OptionKeys::rna::denovo::minimize::protein_packing_distance
		+ OptionKeys::rna::denovo::minimize::minimize_use_nblist
		+ OptionKeys::rna::denovo::minimize::minimize_nblist_auto_update
		+ OptionKeys::rna::denovo::minimize::minimize_max_iter;
}

} //options
} //denovo
} //protocols
