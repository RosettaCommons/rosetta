// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/options/StepWiseBasicModelerOptions.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/options/StepWiseBasicModelerOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

//Req'd on WIN32
#include <utility/tag/Tag.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.options.StepWiseBasicModelerOptions" );

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
	// If you add a variable, initialize it here, and include in operator = definition below!
	void
	StepWiseBasicModelerOptions::initialize_variables(){
		StepWiseBasicOptions::initialize_variables();
		sampler_num_pose_kept_ = 0; // signal to use separate defaults in protein/RNA modeler
		use_green_packer_ = false; // perhaps deprecate
		verbose_ = false; // perhaps deprecate
		choose_random_ = false;
		dump_ = false;
		cluster_rmsd_ = 0.0; // signal for clusterers to use their default value
		skip_minimize_ = false;
		disallow_realign_ = false;
		choose_random_ = false;
	}

	/// @brief clone the options
	StepWiseBasicModelerOptionsOP
	StepWiseBasicModelerOptions::clone() const
	{
		return new StepWiseBasicModelerOptions( *this );
	}

	/////////////////////////////////////////////////////////////////////////////////////
	StepWiseBasicModelerOptions &
	StepWiseBasicModelerOptions::operator = ( StepWiseBasicModelerOptions const & src )
	{
		sampler_num_pose_kept_ = src.sampler_num_pose_kept_;
		use_green_packer_ = src.use_green_packer_;
		verbose_ = src.verbose_;
		choose_random_ = src.choose_random_;
		dump_ = src.dump_;
		cluster_rmsd_ = src.cluster_rmsd_;
		skip_minimize_ = src.skip_minimize_;
		disallow_realign_ = src.disallow_realign_;
		choose_random_ = src.choose_random_;
		return *this;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseBasicModelerOptions::initialize_from_command_line(){

		StepWiseBasicOptions::initialize_from_command_line();

		sampler_num_pose_kept_ = option[ basic::options::OptionKeys::stepwise::rna::sampler_num_pose_kept ]();
		choose_random_ = option[ basic::options::OptionKeys::stepwise::choose_random ]() ;
		if ( choose_random_ && num_pose_minimize() == 0 ) set_num_pose_minimize( 1 );
		verbose_ = option[ OptionKeys::stepwise::VERBOSE ]();
		use_green_packer_ = option[ basic::options::OptionKeys::stepwise::use_green_packer ]();
		choose_random_ = option[ basic::options::OptionKeys::stepwise::choose_random ]() ;
		dump_ = option[ basic::options::OptionKeys::stepwise::dump ]();
		if ( option[ basic::options::OptionKeys::cluster::radius ].user() ) cluster_rmsd_ = option[ basic::options::OptionKeys::cluster::radius ]();
		skip_minimize_ = option[ basic::options::OptionKeys::stepwise::skip_minimize ]();
	}

} //options
} //modeler
} //stepwise
} //protocols
