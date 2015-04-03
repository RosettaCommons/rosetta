// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/copy_dofs/ResidueListStepWiseSampler.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/copy_dofs/ResidueListStepWiseSampler.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>
#include <core/conformation/Residue.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.sampler.rigid_body.ResidueListStepWiseSampler" );

namespace protocols {
namespace stepwise {
namespace sampler {
namespace copy_dofs {

	//Constructor
	ResidueListStepWiseSampler::ResidueListStepWiseSampler( utility::vector1< core::conformation::ResidueOP > copy_dofs ):
		copy_dofs_( copy_dofs )
	{
	}

	//Destructor
	ResidueListStepWiseSampler::~ResidueListStepWiseSampler()
	{}

	core::conformation::ResidueOP
	ResidueListStepWiseSampler::get_residue_at_origin(){
		return copy_dofs_[ id() ];
	}

} //copy_dofs
} //sampler
} //stepwise
} //protocols
