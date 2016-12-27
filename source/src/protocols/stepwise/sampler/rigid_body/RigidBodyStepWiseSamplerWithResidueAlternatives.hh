// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_rigid_body_RigidBodyStepWiseSamplerWithResidueAlternatives_HH
#define INCLUDED_protocols_sampler_rigid_body_RigidBodyStepWiseSamplerWithResidueAlternatives_HH

#include <protocols/stepwise/sampler/StepWiseSamplerComb.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.fwd.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSampler.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSamplerComb.fwd.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rigid_body {

class RigidBodyStepWiseSamplerWithResidueAlternatives: public StepWiseSamplerComb {

public:

	//constructor
	RigidBodyStepWiseSamplerWithResidueAlternatives(
		protocols::stepwise::sampler::copy_dofs::ResidueAlternativeStepWiseSamplerCombOP residue_alternative_rotamer,
		RigidBodyStepWiseSamplerOP rigid_body_rotamer );

	//destructor
	~RigidBodyStepWiseSamplerWithResidueAlternatives();

public:

	void
	fast_forward();

	void
	fast_forward_to_next_rigid_body();

	void
	fast_forward_to_next_translation();

	void
	fast_forward_to_next_euler_gamma();

	ValueList const & get_rigid_body_values();

	// from rigid body rotamer
	core::kinematics::Stub get_stub();

	// from residue list rotamer
	core::conformation::Residue const & get_residue_at_origin();

	// from residue list rotamer
	core::conformation::Residue const & get_residue_at_origin( core::Size const seqpos );

	protocols::stepwise::sampler::copy_dofs::ResidueAlternativeStepWiseSamplerCombOP residue_alternatives_rotamer();
	RigidBodyStepWiseSamplerOP rigid_body_rotamer();

	/// @brief Name of the class
	virtual std::string get_name() const { return "RigidBodyStepWiseSamplerWithResidueAlternatives"; }

	/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
	virtual toolbox::SamplerPlusPlusType type() const { return toolbox::RIGID_BODY_WITH_RESIDUE_ALTERNATIVES; }

	void
	apply_rigid_body_only( core::pose::Pose & pose );

	void
	fast_forward_to_next_residue_pair( core::Size const i, core::Size const j);

	void
	fast_forward_to_next_residue( core::Size const i );

	core::Vector
	get_xyz( core::Size const seqpos, std::string const & atom_name  );

	core::conformation::ResidueCOP
	get_residue( core::Size const seqpos );

private:

	protocols::stepwise::sampler::copy_dofs::ResidueAlternativeStepWiseSamplerCombOP residue_alternatives_rotamer_;
	RigidBodyStepWiseSamplerOP rigid_body_rotamer_;

	std::map< core::Size, core::conformation::ResidueOP > transformed_residues;


};

} //rigid_body
} //sampler
} //stepwise
} //protocols

#endif
