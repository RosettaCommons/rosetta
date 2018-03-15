// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerRingConformer.hh
/// @brief Generate rotamer for one torsion angle.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_stepwise_sampler_StepWiseSamplerRingConformer_HH
#define INCLUDED_stepwise_sampler_StepWiseSamplerRingConformer_HH

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerRingConformer.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerOneValue.hh>

// Project headers
#include <core/id/TorsionID.hh>

#include <core/chemical/rings/RingConformerSet.fwd.hh>

namespace core { namespace chemical { namespace rings { struct RingConformer; } } }

namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSamplerRingConformer : public StepWiseSamplerSized {
public:
	using sampler::StepWiseSampler::TorsionList;

	StepWiseSamplerRingConformer();

	StepWiseSamplerRingConformer(
		core::Size const ring_num,
		core::Size const moving_rsd,
		core::chemical::rings::RingConformerSetCOP ring_conformers
	);

	virtual ~StepWiseSamplerRingConformer();

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & pose ) { apply( pose, id_ ); }

	/// @brief Apply the i-th rotamer to pose
	virtual void apply( core::pose::Pose & pose, Size const i );

	/// @brief Set the allowed torsions in sampler
	// virtual void set_torsions( TorsionList const & setting ) {
	//  set_values( setting );
	// }

	/// @brief Set the residue id being sampled
	// virtual void set_rsd_id( core::Size const setting ) {
	//  torsion_id_ =
	//   core::id::TorsionID( setting, torsion_id_.type(),
	//   torsion_id_.torsion() );
	// }

	/// @brief Set the TorsionID of the sampler
	// virtual void set_torsion_id( core::id::TorsionID const & setting ) {
	//  torsion_id_ = setting;
	// }

	core::Size size() const;

	/// @brief Name of the class
	virtual std::string get_name() const;

	/// @brief Type of class (see enum in SamplerPlusPlusTypes.hh)
	virtual toolbox::SamplerPlusPlusType type() const { return toolbox::RING_CONFORMERS; }

private:

	core::Size ring_num_ = 0;
	core::Size moving_rsd_ = 0;
	utility::vector1< core::chemical::rings::RingConformer > ring_conformers_;

};

} //sampler
} //stepwise
} //protocols

#endif
