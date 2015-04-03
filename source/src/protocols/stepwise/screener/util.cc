// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/StepWiseScreenerUtil.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/util.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerBase.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueList.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.screener.util" );

namespace protocols {
namespace stepwise {
namespace screener {

	///////////////////////////////////////////////////////////////////////////////////////
	// used by StepWiseResiduePairScreener, and also ChainClosureScreener.
	void
	fast_forward_to_next_residue_pair( sampler::StepWiseSamplerBaseOP sampler,
																		 Size const res1,
																		 Size const res2 ){

		using namespace sampler;
		using namespace sampler::rigid_body;

		if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_LIST ){
			RigidBodyStepWiseSamplerWithResidueList & rigid_body_rotamer_with_residue_list = *( static_cast< RigidBodyStepWiseSamplerWithResidueList * >( sampler.get() ) );
			rigid_body_rotamer_with_residue_list.fast_forward_to_next_rigid_body();
			return;
		} else if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ){
			RigidBodyStepWiseSamplerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyStepWiseSamplerWithResidueAlternatives * >( sampler.get() ) );
			rigid_body_rotamer_with_residue_alternatives.fast_forward_to_next_residue_pair( res1, res2 );
			return;
		}

	}

} //screener
} //stepwise
} //protocols
