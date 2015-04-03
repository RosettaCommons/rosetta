// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/RNA_ChainClosureScreener.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/RNA_ChainClosureScreener.hh>
#include <protocols/stepwise/screener/util.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosureChecker.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerBase.hh>
#include <protocols/moves/CompositionMover.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.screener.RNA_ChainClosureScreener" );

using namespace core;


namespace protocols {
namespace stepwise {
namespace screener {

  RNA_ChainClosureScreener::RNA_ChainClosureScreener( modeler::rna::checker::RNA_ChainClosureCheckerOP chain_closure_checker ):
		SampleApplier( chain_closure_checker->pose() ),
		chain_closure_checker_( chain_closure_checker ),
		screening_pose_( chain_closure_checker->pose() ),
		just_do_closure_check_( false )
	{}

	//Constructor -- this one is useful if we don't want to close the chain but check if it is already closed.
  RNA_ChainClosureScreener::RNA_ChainClosureScreener( modeler::rna::checker::RNA_ChainClosureCheckerOP chain_closure_checker,
																							pose::Pose & screening_pose,
																							bool const just_do_closure_check /*= false */ ):
		SampleApplier( chain_closure_checker->pose() ),
		chain_closure_checker_( chain_closure_checker ),
		screening_pose_( screening_pose ),
		just_do_closure_check_( just_do_closure_check )
	{}

	//Destructor
	RNA_ChainClosureScreener::~RNA_ChainClosureScreener()
	{}

	/////////////////////////////////////////
	bool
	RNA_ChainClosureScreener::check_screen() {
		if ( just_do_closure_check_ ) return chain_closure_checker_->check_loop_closed( screening_pose_ );
		return chain_closure_checker_->check_screen( screening_pose_ );
	}

	/////////////////////////////////////////
	void
	RNA_ChainClosureScreener::add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover ){
		update_mover->add_mover( chain_closure_checker_ );
		restore_mover->add_mover( 0 );
	}

	////////////////////////////////////////////////////////////////////////////
	// this is also used in StepWiseResiduePairScreener and in principle this class
	//  could derive from that parent, but we also want SampleApplier functionality...
	void
	RNA_ChainClosureScreener::fast_forward( sampler::StepWiseSamplerBaseOP sampler ){
		fast_forward_to_next_residue_pair( sampler,
																			 chain_closure_checker_->five_prime_res(),
																			 chain_closure_checker_->five_prime_res() + 1); // in screener util.
	}

} //screener
} //stepwise
} //protocols
