// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/modeler_options/StepWiseBasicModelerOptions.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/modeler_options/StepWiseBasicModelerOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

//Req'd on WIN32
#include <utility/tag/Tag.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TR( "protocols.stepwise.sampling.modeler_options.StepWiseBasicModelerOptions" );

namespace protocols {
namespace stepwise {
namespace sampling {
namespace modeler_options {

	//Constructor
	StepWiseBasicModelerOptions::StepWiseBasicModelerOptions()
	{
		initialize_variables();
	}

	//Destructor
	StepWiseBasicModelerOptions::~StepWiseBasicModelerOptions()
	{}

	/// @brief copy constructor
	StepWiseBasicModelerOptions::StepWiseBasicModelerOptions( StepWiseBasicModelerOptions const & src ) :
		ResourceOptions ( src )
	{
		*this = src;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	// If you add a variable, initialize it here, and include in operator = definition below!
	void
	StepWiseBasicModelerOptions::initialize_variables(){
		silent_file_ = "";
		sampler_num_pose_kept_ = 0; // signal to use separate defaults in protein/RNA sampling
		num_pose_minimize_ = 0; // signal to minimize all.
		use_green_packer_ = false; // perhaps deprecate
		verbose_ = false; // perhaps deprecate
		choose_random_ = false;
		num_random_samples_ = 20;
		dump_ = false;
		atr_rep_screen_ = true;
		rmsd_screen_ = 0.0;
		cluster_rmsd_ = 0.0; // signal for clusterers to use their default value
		skip_minimize_ = false;
		output_minimized_pose_list_ = false;
		disallow_realign_ = false;
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
		silent_file_ = src.silent_file_;
		sampler_num_pose_kept_ = src.sampler_num_pose_kept_;
		num_pose_minimize_ = src.num_pose_minimize_;
		use_green_packer_ = src.use_green_packer_;
		verbose_ = src.verbose_;
		num_random_samples_ = src.num_random_samples_;
		choose_random_ = src.choose_random_;
		dump_ = src.dump_;
		atr_rep_screen_ = src.atr_rep_screen_;
		rmsd_screen_ = src.rmsd_screen_;
		cluster_rmsd_ = src.cluster_rmsd_;
		skip_minimize_ = src.skip_minimize_;
		output_minimized_pose_list_ = src.output_minimized_pose_list_;
		disallow_realign_ = src.disallow_realign_;
		return *this;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseBasicModelerOptions::initialize_from_command_line(){
		sampler_num_pose_kept_ = option[ basic::options::OptionKeys::stepwise::rna::sampler_num_pose_kept ]();
		if ( option[ basic::options::OptionKeys::stepwise::choose_random ]() ) num_pose_minimize_ = 1;
		if ( option[ basic::options::OptionKeys::stepwise::num_pose_minimize ].user() ) num_pose_minimize_ = option[ basic::options::OptionKeys::stepwise::num_pose_minimize ]();
		verbose_ = option[ OptionKeys::stepwise::VERBOSE ]();
		use_green_packer_ = option[ basic::options::OptionKeys::stepwise::use_green_packer ]();
		choose_random_ = option[ basic::options::OptionKeys::stepwise::choose_random ]() ;
		num_random_samples_ = option[ basic::options::OptionKeys::stepwise::num_random_samples ]();
		dump_ = option[ basic::options::OptionKeys::stepwise::dump ]();
		atr_rep_screen_ = option[ basic::options::OptionKeys::stepwise::atr_rep_screen ]();
		rmsd_screen_ = option[ basic::options::OptionKeys::stepwise::rmsd_screen ]();
		if ( option[ basic::options::OptionKeys::cluster::radius ].user() ) cluster_rmsd_ = option[ basic::options::OptionKeys::cluster::radius ]();
		skip_minimize_ = option[ basic::options::OptionKeys::stepwise::skip_minimize ]();
	}

} //modeler_options
} //sampling
} //stepwise
} //protocols
