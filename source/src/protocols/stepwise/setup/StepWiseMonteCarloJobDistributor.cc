// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/setup/StepWiseMonteCarloJobDistributor.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/setup/StepWiseMonteCarloJobDistributor.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <core/io/silent/util.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.setup.StepWiseMonteCarloJobDistributor" );

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Bare-bones job distributor for StepWiseMonteCarlo. Also should work for single moves.
//
// Just looks to see if tag S_0, S_1, etc. is already in silent_file up to
//  nstruct.
//
// Written as a separated class to allow alternative job setup, e.g., conformational space annealing,
//   genetic algorithms, or other 'population-based' approaches to optimization.
//
// Could perhaps be combined with Rosetta's standard JobDistributor -- would just take some time to unify.
//
//                       -- rhiju, 2014.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace setup {

	//Constructor
 	StepWiseMonteCarloJobDistributor::StepWiseMonteCarloJobDistributor( stepwise::monte_carlo::StepWiseMonteCarloOP stepwise_monte_carlo,
																																			std::string const silent_file,
																																			core::Size const nstruct ):
		StepWiseJobDistributor( stepwise_monte_carlo, silent_file, nstruct ),
		count_( 1 ),
		out_tag_( "" ),
		init_tag_is_done_( false )
	{
	}

	//Destructor
	StepWiseMonteCarloJobDistributor::~StepWiseMonteCarloJobDistributor()
	{}

	void
	StepWiseMonteCarloJobDistributor::initialize( core::pose::Pose const & pose ) {
		start_pose_ = pose.clone();
		count_ = 1;
	}

	void
	StepWiseMonteCarloJobDistributor::move_forward_to_next_model() {
		while ( count_ <= nstruct_ && !get_out_tag() ){
			count_++;
			if ( !get_out_tag() ) TR << "Already done: " << out_tag_ << std::endl;
		}
	}

	bool
	StepWiseMonteCarloJobDistributor::has_another_job() {
		move_forward_to_next_model();
		return ( count_ <= nstruct_ );
	}

	void
	StepWiseMonteCarloJobDistributor::apply( core::pose::Pose & pose ) {

		runtime_assert( has_another_job() );
		TR << std::endl << TR.Green << "Embarking on structure " << count_ << " of " << nstruct_ << TR.Reset << std::endl;

		runtime_assert( start_pose_ != 0 ); // initialized.
		pose = *start_pose_;
		stepwise_monte_carlo_->set_model_tag( out_tag_ );

 		stepwise_monte_carlo_->apply( pose );

		if (! stepwise_monte_carlo_->options()->output_minimized_pose_list()){
			monte_carlo::output_to_silent_file( out_tag_, silent_file_,
														 pose, get_native_pose(), superimpose_over_all_,
														 true /*rms_fill*/ );
		}
		tag_is_done_[ out_tag_ ] = true;
	}

	///////////////////////////////////////////////////////////////////////////////
	// used to be in a util.cc
	bool
	StepWiseMonteCarloJobDistributor::get_out_tag(){
		using namespace core::io::silent;
		if ( !init_tag_is_done_ ){
			tag_is_done_ = initialize_tag_is_done( silent_file_ );
			init_tag_is_done_ = true;
		}

		out_tag_ = "S_" + ObjexxFCL::lead_zero_string_of( count_, 6 );
		if ( tag_is_done_[ out_tag_ ] ) {
			return false;
		}
		return true; //ok, not done, so do it.
	}


} //setup
} //stepwise
} //protocols
