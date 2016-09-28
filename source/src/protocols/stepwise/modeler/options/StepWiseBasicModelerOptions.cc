// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/options/StepWiseBasicModelerOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/options/StepWiseBasicModelerOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

//Req'd on WIN32
#include <utility/tag/Tag.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.options.StepWiseBasicModelerOptions" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace options {

//Constructor
StepWiseBasicModelerOptions::StepWiseBasicModelerOptions():
	StepWiseBasicOptions()
{
	initialize_variables();
}

//Destructor
StepWiseBasicModelerOptions::~StepWiseBasicModelerOptions()
{}

/// @brief copy constructor
StepWiseBasicModelerOptions::StepWiseBasicModelerOptions( StepWiseBasicModelerOptions const & src ) :
	ResourceOptions( src ),
	StepWiseBasicOptions( src )
{
	*this = src;
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseBasicModelerOptions::initialize_variables(){
	StepWiseBasicOptions::initialize_variables();
	use_green_packer_ = false; // perhaps deprecate
	verbose_ = false; // perhaps deprecate
	choose_random_ = false;
	dump_ = false;
	skip_minimize_ = false;
	disallow_realign_ = false;
	coordinate_constraints_during_minimize_ = true;
	virtualize_packable_moieties_in_screening_pose_ = false;
}

/// @brief clone the options
StepWiseBasicModelerOptionsOP
StepWiseBasicModelerOptions::clone() const
{
	return StepWiseBasicModelerOptionsOP( new StepWiseBasicModelerOptions( *this ) );
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseBasicModelerOptions::initialize_from_command_line(){

	StepWiseBasicOptions::initialize_from_command_line();

	choose_random_ = option[ basic::options::OptionKeys::stepwise::choose_random ]() ;
	if ( choose_random_ && num_pose_minimize() == 0 ) set_num_pose_minimize( 1 );
	verbose_ = option[ OptionKeys::stepwise::VERBOSE ]();
	use_green_packer_ = option[ basic::options::OptionKeys::stepwise::use_green_packer ]();
	choose_random_ = option[ basic::options::OptionKeys::stepwise::choose_random ]() ;
	dump_ = option[ basic::options::OptionKeys::stepwise::dump ]();
	skip_minimize_ = option[ basic::options::OptionKeys::stepwise::skip_minimize ]();
	virtualize_packable_moieties_in_screening_pose_ = option[ basic::options::OptionKeys::stepwise::virtualize_packable_moieties_in_screening_pose ]();
}

} //options
} //modeler
} //stepwise
} //protocols
