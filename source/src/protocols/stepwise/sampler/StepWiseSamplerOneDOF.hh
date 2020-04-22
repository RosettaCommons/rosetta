// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerOneDOF.hh
/// @brief Generate rotamer for one DOF angle.
/// @author Andy Watkins


#ifndef INCLUDED_stepwise_sampler_StepWiseSamplerOneDOF_HH
#define INCLUDED_stepwise_sampler_StepWiseSamplerOneDOF_HH

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerOneDOF.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerOneValue.hh>
#include <core/id/AtomID.hh>

// Project headers
#include <core/id/DOF_ID.hh>

namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSamplerOneDOF : public StepWiseSamplerOneValue {
public:
	using sampler::StepWiseSampler::TorsionList;

	StepWiseSamplerOneDOF();

	StepWiseSamplerOneDOF(
		core::id::DOF_ID const & tor_id,
		TorsionList const & DOFs
	);

	~StepWiseSamplerOneDOF() override;

	/// @brief Apply the current rotamer to pose //
	void apply( core::pose::Pose & pose ) override { apply( pose, id_ ); }

	/// @brief Apply the i-th rotamer to pose
	void apply( core::pose::Pose & pose, core::Size const i ) override;

	/// @brief Set the allowed DOFs in sampler
	virtual void set_DOFs( TorsionList const & setting ) {
		set_values( setting );
	}

	/// @brief Set the atom id being sampled
	virtual void set_atom_id( core::id::AtomID const & setting ) {
		DOF_id_ = core::id::DOF_ID( setting, DOF_id_.type() );
	}

	/// @brief Set the DOF_ID of the sampler
	virtual void set_DOF_id( core::id::DOF_ID const & setting ) {
		DOF_id_ = setting;
	}

	/// @brief Name of the class
	std::string get_name() const override;

	/// @brief Type of class (see enum in SamplerPlusPlusTypes.hh)
	toolbox::SamplerPlusPlusType type() const override { return toolbox::ONE_DOF; }

private:

	core::id::DOF_ID DOF_id_;

};

} //sampler
} //stepwise
} //protocols

#endif
