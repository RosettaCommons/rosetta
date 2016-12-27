// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/protein/ProteinBetaAntiParallelStepWiseSampler.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_protein_ProteinBetaAntiParallelStepWiseSampler_HH
#define INCLUDED_protocols_sampler_protein_ProteinBetaAntiParallelStepWiseSampler_HH

#include <protocols/stepwise/sampler/jump/JumpStepWiseSampler.hh>
#include <protocols/stepwise/sampler/protein/ProteinBetaAntiParallelStepWiseSampler.fwd.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace protein {

class ProteinBetaAntiParallelStepWiseSampler: public jump::JumpStepWiseSampler {

public:

	//constructor
	ProteinBetaAntiParallelStepWiseSampler( core::pose::Pose const & pose,
		Size const moving_residue );

	//destructor
	~ProteinBetaAntiParallelStepWiseSampler();

public:

	/// @brief Name of the class
	virtual std::string get_name() const { return "ProteinBetaAntiParallelStepWiseSampler"; }

	/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
	virtual toolbox::SamplerPlusPlusType type() const { return toolbox::PROTEIN_BETA_ANTIPARALLEL; }

private:

	Size
	get_antiparallel_beta_jumps( core::pose::Pose const & pose, int const sample_res );

};

} //protein
} //sampler
} //stepwise
} //protocols

#endif
