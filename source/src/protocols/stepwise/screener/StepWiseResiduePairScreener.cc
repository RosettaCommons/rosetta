// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/StepWiseResiduePairScreener.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/StepWiseResiduePairScreener.hh>
#include <protocols/stepwise/screener/util.hh>
#include <protocols/stepwise/sampler/StepWiseSampler.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.StepWiseResiduePairScreener" );

namespace protocols {
namespace stepwise {
namespace screener {

//Constructor
StepWiseResiduePairScreener::StepWiseResiduePairScreener( Size const res1, Size const res2 ):
	res1_( res1 ),
	res2_( res2 )
{}

//Destructor
StepWiseResiduePairScreener::~StepWiseResiduePairScreener()
{}

////////////////////////////////////////////////////////////////////////////
void
StepWiseResiduePairScreener::fast_forward( sampler::StepWiseSamplerOP sampler ){
	fast_forward_to_next_residue_pair( sampler, res1_, res2_ ); // in screener util.
}

} //screener
} //stepwise
} //protocols
