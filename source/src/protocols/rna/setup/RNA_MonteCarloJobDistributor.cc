// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/setup/RNA_MonteCarloJobDistributor.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rna/setup/RNA_MonteCarloJobDistributor.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/rna/denovo/RNA_FragmentMonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/rna/denovo/RNA_DeNovoProtocol.hh>
#include <core/io/silent/util.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <utility/file/FileName.hh>

static basic::Tracer TR( "protocols.rna.setup.RNA_MonteCarloJobDistributor" );

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Bare-bones job distributor for RNA_MonteCarlo. Also should work for single moves.
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
namespace rna {
namespace setup {

//Constructor
RNA_MonteCarloJobDistributor::RNA_MonteCarloJobDistributor(
	stepwise::monte_carlo::StepWiseMonteCarloOP stepwise_monte_carlo,
	std::string const & silent_file,
	core::Size const nstruct ):
	RNA_JobDistributor( stepwise_monte_carlo, silent_file, nstruct ),
	count_( 1 ),
	out_tag_( "" ),
	init_tag_is_done_( false )
{
}

RNA_MonteCarloJobDistributor::RNA_MonteCarloJobDistributor(
	protocols::rna::denovo::RNA_FragmentMonteCarloOP rna_fragment_monte_carlo,
	std::string const & silent_file,
	core::Size const nstruct ):
	RNA_JobDistributor( rna_fragment_monte_carlo, silent_file, nstruct ),
	count_( 1 ),
	out_tag_( "" ),
	init_tag_is_done_( false )
{
}

//Destructor
RNA_MonteCarloJobDistributor::~RNA_MonteCarloJobDistributor() = default;

void
RNA_MonteCarloJobDistributor::initialize( core::pose::Pose const & pose ) {
	start_pose_ = pose.clone();
	count_ = 1;
}

void
RNA_MonteCarloJobDistributor::move_forward_to_next_model() {
	while ( count_ <= nstruct_ && !get_out_tag() ) {
		count_++;
		if ( !get_out_tag() ) TR << "Already done: " << out_tag_ << std::endl;
	}
}

bool
RNA_MonteCarloJobDistributor::has_another_job() {
	move_forward_to_next_model();
	return ( count_ <= nstruct_ );
}

void
RNA_MonteCarloJobDistributor::apply( core::pose::Pose & pose ) {

	runtime_assert( stepwise_monte_carlo_ || rna_fragment_monte_carlo_ );

	runtime_assert( has_another_job() );
	TR << std::endl << TR.Green << "Embarking on structure " << count_ << " of " << nstruct_ << TR.Reset << std::endl;

	runtime_assert( start_pose_ != nullptr ); // initialized.
	pose = *start_pose_;

	// AMW TODO unify -- right now we just switch on pointer nullity to do one of two behaviors!
	if ( stepwise_monte_carlo_ ) {
		stepwise_monte_carlo_->set_model_tag( out_tag_ );
		stepwise_monte_carlo_->set_out_file_prefix( utility::file::FileName( silent_file_ ).base() );

		// Maybe here is where we check for a possibly checkpointed file for structure out_tag_

		stepwise_monte_carlo_->apply( pose );

		if ( ! stepwise_monte_carlo_->options()->output_minimized_pose_list() ) {
			stepwise::monte_carlo::output_to_silent_file( out_tag_, silent_file_,
				pose, get_native_pose(), superimpose_over_all_,
				true /*rms_fill*/ );
		}
	} else {
		rna_fragment_monte_carlo_->set_out_file_tag( out_tag_ );
		rna_fragment_monte_carlo_->apply( pose );




	}
	tag_is_done_[ out_tag_ ] = true;
}

///////////////////////////////////////////////////////////////////////////////
// used to be in a util.cc
bool
RNA_MonteCarloJobDistributor::get_out_tag(){
	using namespace core::io::silent;
	if ( !init_tag_is_done_ ) {
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
} //rna
} //protocols
