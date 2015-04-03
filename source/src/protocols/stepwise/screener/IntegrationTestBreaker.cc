// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/IntegrationTestBreaker.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/stepwise/screener/IntegrationTestBreaker.hh>
#include <protocols/stepwise/screener/AlignRMSD_Screener.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueList.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.screener.IntegrationTestBreaker" );

using namespace protocols::stepwise::sampler;
using namespace protocols::stepwise::sampler::rigid_body;

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	IntegrationTestBreaker::IntegrationTestBreaker( StepWiseScreenerOP screener_whose_counts_to_check,
																									StepWiseScreenerOP final_screener /*total_count -- for turning on align screen*/,
																									AlignRMSD_ScreenerOP align_rmsd_screener ):
		screener_whose_counts_to_check_( screener_whose_counts_to_check ),
		final_screener_( final_screener ),
		align_rmsd_screener_( align_rmsd_screener )
	{}

	//Destructor
	IntegrationTestBreaker::~IntegrationTestBreaker()
	{}

	//////////////////////////////////////////////////////////////////////////
	// If you see an issue with changes to integration tests, take
	// a close look at options -- note that it should 'spawn' get_sampler_options()
	// which is fixed for integration tests.
	bool
	IntegrationTestBreaker::check_screen() {
		if ( screener_whose_counts_to_check_ && screener_whose_counts_to_check_->count() >= 100 ) 	return false;

		if ( align_rmsd_screener_ != 0 ){
			if ( final_screener_->count() >= 10 ) align_rmsd_screener_->set_do_screen( true );
			if ( align_rmsd_screener_->pass_count() >= 10 ) return false;
		}
		return true;
	}

	////////////////////////////////////////////////////////////////////////////
	void
	IntegrationTestBreaker::fast_forward( sampler::StepWiseSamplerBaseOP sampler ){
		if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_LIST ){
			RigidBodyStepWiseSamplerWithResidueList & rigid_body_rotamer_with_copy_dofs = *( static_cast< RigidBodyStepWiseSamplerWithResidueList * >( sampler.get() ) );
			rigid_body_rotamer_with_copy_dofs.fast_forward();
		} else if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ){
			RigidBodyStepWiseSamplerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyStepWiseSamplerWithResidueAlternatives * >( sampler.get() ) );
			rigid_body_rotamer_with_residue_alternatives.fast_forward();
		}
	}

} //screener
} //stepwise
} //protocols
