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
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.import_pose.options.RNA_MinimizerOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

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
	protein_packing_distance_( 10.0 )
{}

//Destructor
RNA_MinimizerOptions::~RNA_MinimizerOptions() = default;


/// @brief copy constructor
RNA_MinimizerOptions::RNA_MinimizerOptions( RNA_MinimizerOptions const & src ) :
	ResourceOptions( src ),
	RNA_BasicOptions( src )
{
	*this = src;
}

/// @brief clone the options
RNA_MinimizerOptionsOP
RNA_MinimizerOptions::clone() const
{
	return RNA_MinimizerOptionsOP( new RNA_MinimizerOptions( *this ) );
}

///////////////////////////////////////////////////////////////////
void
RNA_MinimizerOptions::initialize_from_command_line() {
	initialize_from_options( basic::options::option );
}

void
RNA_MinimizerOptions::initialize_from_options( utility::options::OptionCollection const & opts ) {

	RNA_BasicOptions::initialize_from_options( opts );

	set_minimize_rounds( opts[ OptionKeys::rna::denovo::minimize::minimize_rounds ]() );
	set_vary_bond_geometry( opts[ OptionKeys::rna::vary_geometry ]() );
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
}

void
RNA_MinimizerOptions::list_options_read( utility::options::OptionKeyList & opts ) {
	using namespace basic::options;
	RNA_BasicOptions::list_options_read( opts );
	opts + OptionKeys::rna::denovo::minimize::minimize_rounds
		+ OptionKeys::rna::vary_geometry
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
		+ OptionKeys::rna::denovo::minimize::protein_packing_distance;
}

} //options
} //denovo
} //protocols
