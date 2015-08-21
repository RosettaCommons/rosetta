// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/copy_dofs/ResidueListStepWiseSampler.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_rigid_body_ResidueListStepWiseSampler_HH
#define INCLUDED_protocols_sampler_rigid_body_ResidueListStepWiseSampler_HH

#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueListStepWiseSampler.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace copy_dofs {

class ResidueListStepWiseSampler: public StepWiseSamplerSized {

public:

	//constructor
	ResidueListStepWiseSampler( utility::vector1< core::conformation::ResidueOP > copy_dofs );

	//destructor
	~ResidueListStepWiseSampler();

public:

	core::conformation::ResidueOP get_residue_at_origin();

	/// @brief Get the total number of rotamers in sampler
	virtual core::Size size() const{ return copy_dofs_.size(); }

	/// @brief Apply the i-th rotamer to pose
	// do nothing. job is to return a residue.
	virtual void apply( core::pose::Pose&, core::Size const ){}

	/// @brief Name of the class
	virtual std::string get_name() const { return "ResidueListStepWiseSampler"; }

	/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
	virtual StepWiseSamplerType type() const { return RESIDUE_LIST; }

private:

	utility::vector1< core::conformation::ResidueOP > copy_dofs_;

};

} //copy_dofs
} //sampler
} //stepwise
} //protocols

#endif
