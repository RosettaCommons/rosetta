// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/mover/options/StepWiseMoveSelectorOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/mover/options/StepWiseMoveSelectorOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.monte_carlo.mover.options.StepWiseMoveSelectorOptions" );

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {
namespace options {

//Constructor
StepWiseMoveSelectorOptions::StepWiseMoveSelectorOptions():
	allow_internal_hinge_moves_( true ),
	allow_internal_local_moves_( true ),
	add_delete_frequency_( 0.5 ),
	from_scratch_frequency_( 0.0 ),
	docking_frequency_( 0.2 ),
	submotif_frequency_( 0.2 ),
	switch_focus_frequency_( 0.2 ),
	skip_bulge_frequency_( 0.0 ),
	vary_loop_length_frequency_( 0.0 ),
	filter_complex_cycles_( true ),
	allow_submotif_split_( false ),
	force_submotif_without_intervening_bulge_( false )
{
}

//Destructor
StepWiseMoveSelectorOptions::~StepWiseMoveSelectorOptions()
{}

/// @brief copy constructor
StepWiseMoveSelectorOptions::StepWiseMoveSelectorOptions( StepWiseMoveSelectorOptions const & src ) :
	ResourceOptions( src )
{
	*this = src;
}

/// @brief clone the options
StepWiseMoveSelectorOptionsOP
StepWiseMoveSelectorOptions::clone() const
{
	return StepWiseMoveSelectorOptionsOP( new StepWiseMoveSelectorOptions( *this ) );
}

///////////////////////////////////////////////////////////////////
void
StepWiseMoveSelectorOptions::initialize_from_command_line() {
	initialize_from_options_collection( option );
}

void
StepWiseMoveSelectorOptions::initialize_from_options_collection( utility::options::OptionCollection const & the_options ) {
	set_allow_internal_hinge_moves( the_options[ OptionKeys::stepwise::monte_carlo::allow_internal_hinge_moves ]() );
	set_allow_internal_local_moves( the_options[ OptionKeys::stepwise::monte_carlo::allow_internal_local_moves ]() );
	set_add_delete_frequency( the_options[ OptionKeys::stepwise::monte_carlo::add_delete_frequency ]() );
	set_from_scratch_frequency( the_options[ OptionKeys::stepwise::monte_carlo::from_scratch_frequency ]() );

	set_docking_frequency( the_options[ OptionKeys::stepwise::monte_carlo::docking_frequency ]() );
	if ( the_options[ OptionKeys::stepwise::monte_carlo::intermolecular_frequency ].user() ) {
		TR << TR.Red << "Use -docking_frequency instead of -intermolecular_frequency -- will be deprecated soon." << TR.Reset << std::endl;
		set_docking_frequency( the_options[ OptionKeys::stepwise::monte_carlo::intermolecular_frequency ]() );
	}
	if ( the_options[ OptionKeys::stepwise::monte_carlo::allow_skip_bulge ].user() ) { // old option.
		set_skip_bulge_frequency( the_options[ OptionKeys::stepwise::monte_carlo::allow_skip_bulge ]() ? 0.2 : 0.0 );
	}

	set_submotif_frequency( the_options[ OptionKeys::stepwise::monte_carlo::submotif_frequency ]() );
	set_switch_focus_frequency( the_options[ OptionKeys::stepwise::monte_carlo::switch_focus_frequency ]() );
	if ( the_options[ OptionKeys::stepwise::monte_carlo::skip_bulge_frequency ].user() ) set_skip_bulge_frequency( the_options[ OptionKeys::stepwise::monte_carlo::skip_bulge_frequency ]() );
	set_vary_loop_length_frequency( the_options[ OptionKeys::stepwise::monte_carlo::vary_loop_length_frequency ]() );

	// hey what about force_unique_moves?
	filter_complex_cycles_ = !basic::options::option[ basic::options::OptionKeys::score::loop_close::allow_complex_loop_graph ]();
	allow_submotif_split_ = the_options[ OptionKeys::stepwise::monte_carlo::allow_submotif_split ]();
	force_submotif_without_intervening_bulge_ = the_options[ OptionKeys::stepwise::monte_carlo::force_submotif_without_intervening_bulge ]();
}

void
StepWiseMoveSelectorOptions::list_options_read( utility::options::OptionKeyList & opts ) {
	using namespace basic::options;
	opts + OptionKeys::stepwise::monte_carlo::allow_internal_hinge_moves
		+ OptionKeys::stepwise::monte_carlo::allow_internal_local_moves
		+ OptionKeys::stepwise::monte_carlo::add_delete_frequency
		+ OptionKeys::stepwise::monte_carlo::from_scratch_frequency
		+ OptionKeys::stepwise::monte_carlo::docking_frequency
		+ OptionKeys::stepwise::monte_carlo::intermolecular_frequency
		+ OptionKeys::stepwise::monte_carlo::allow_skip_bulge
		+ OptionKeys::stepwise::monte_carlo::submotif_frequency
		+ OptionKeys::stepwise::monte_carlo::switch_focus_frequency
		+ OptionKeys::stepwise::monte_carlo::skip_bulge_frequency
		+ OptionKeys::stepwise::monte_carlo::vary_loop_length_frequency
		+ OptionKeys::score::loop_close::allow_complex_loop_graph
		+ OptionKeys::stepwise::monte_carlo::allow_submotif_split
		+ OptionKeys::stepwise::monte_carlo::force_submotif_without_intervening_bulge;
}

} //options
} //mover
} //monte_carlo
} //stepwise
} //protocols
