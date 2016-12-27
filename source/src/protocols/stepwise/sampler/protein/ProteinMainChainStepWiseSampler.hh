// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/protein/ProteinMainChainStepWiseSampler.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_protein_ProteinMainChainStepWiseSampler_HH
#define INCLUDED_protocols_sampler_protein_ProteinMainChainStepWiseSampler_HH

#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>
#include <protocols/stepwise/sampler/protein/ProteinMainChainStepWiseSampler.fwd.hh>
#include <core/types.hh>
#include <core/id/TorsionID.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace protein {

class ProteinMainChainStepWiseSampler: public StepWiseSamplerSized {

public:

	//constructor
	ProteinMainChainStepWiseSampler( utility::vector1< core::id::TorsionID > const & which_torsions,
		utility::vector1< utility::vector1< core::Real > > const & main_chain_torsion_set_lists,
		bool const choose_random = false );

	ProteinMainChainStepWiseSampler();

	//destructor
	~ProteinMainChainStepWiseSampler();

	/// @brief Get the total number of rotamers in sampler
	virtual core::Size size() const { return main_chain_torsion_set_lists_.size(); }

	/// @brief Apply the i-th rotamer to pose
	virtual void apply( core::pose::Pose&, core::Size const );

	/// @brief Name of the class
	virtual std::string get_name() const { return "ProteinMainChainStepWiseSampler"; }

	/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
	virtual toolbox::SamplerPlusPlusType type() const { return toolbox::PROTEIN_MAIN_CHAIN; }

	/// @brief return OP to the subsampler that controls exactly this torsion_id (assume only one).
	virtual
	toolbox::SamplerPlusPlusOP
	find( core::id::TorsionID const & torsion_id );

private:

	utility::vector1< core::id::TorsionID > const which_torsions_;
	utility::vector1< utility::vector1< core::Real > > const main_chain_torsion_set_lists_;

};

} //protein
} //sampler
} //stepwise
} //protocols

#endif
