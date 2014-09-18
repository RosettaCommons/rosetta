// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/options/StepWiseBasicOptions.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/options/StepWiseBasicOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer TR( "protocols.stepwise.options.StepWiseBasicOptions" );

namespace protocols {
namespace stepwise {
namespace options {

	//Constructor
	StepWiseBasicOptions::StepWiseBasicOptions()
	{
		initialize_variables();
	}

	//Destructor
	StepWiseBasicOptions::~StepWiseBasicOptions()
	{}

	/// @brief copy constructor
	StepWiseBasicOptions::StepWiseBasicOptions( StepWiseBasicOptions const & src ) :
		ResourceOptions ( src )
	{
		*this = src;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	// If you add a variable, initialize it here, and include in operator = definition below!
	void
	StepWiseBasicOptions::initialize_variables(){
		silent_file_ = "";
		sampler_silent_file_ = "";
		num_pose_minimize_ = 0; // signal to minimize all.
		num_random_samples_ = 20;
		atr_rep_screen_ = true;
		rmsd_screen_ = 0.0;
		output_minimized_pose_list_ = false;
	}

	/// @brief clone the options
	StepWiseBasicOptionsOP
	StepWiseBasicOptions::clone() const
	{
		return new StepWiseBasicOptions( *this );
	}

	/////////////////////////////////////////////////////////////////////////////////////
	StepWiseBasicOptions &
	StepWiseBasicOptions::operator = ( StepWiseBasicOptions const & src )
	{
		silent_file_ = src.silent_file_;
		sampler_silent_file_ = src.sampler_silent_file_;
		num_pose_minimize_ = src.num_pose_minimize_;
		num_random_samples_ = src.num_random_samples_;
		atr_rep_screen_ = src.atr_rep_screen_;
		rmsd_screen_ = src.rmsd_screen_;
		output_minimized_pose_list_ = src.output_minimized_pose_list_;
		return *this;
	}


	/////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseBasicOptions::initialize_from_command_line(){
		sampler_silent_file_ = option[ basic::options::OptionKeys::stepwise::sampler_silent_file ]();
		num_random_samples_ = option[ basic::options::OptionKeys::stepwise::num_random_samples ]();
		atr_rep_screen_ = option[ basic::options::OptionKeys::stepwise::atr_rep_screen ]();
		rmsd_screen_ = option[ basic::options::OptionKeys::stepwise::rmsd_screen ]();
		if ( option[ basic::options::OptionKeys::stepwise::num_pose_minimize ].user() ) num_pose_minimize_ = option[ basic::options::OptionKeys::stepwise::num_pose_minimize ]();
	}


} //options
} //stepwise
} //protocols
