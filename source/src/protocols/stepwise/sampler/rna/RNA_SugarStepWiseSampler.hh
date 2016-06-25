// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/rna/RNA_SugarStepWiseSampler.hh
/// @brief Generate sugar pucker rotamers for RNA.
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_sampler_rna_RNA_SugarStepWiseSampler_HH
#define INCLUDED_protocols_sampler_rna_RNA_SugarStepWiseSampler_HH

// Unit headers
#include <protocols/stepwise/sampler/rna/RNA_SugarStepWiseSampler.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>
#include <core/chemical/rna/util.hh>
#include <core/id/DOF_ID_Map.hh>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core::chemical::rna;

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

class RNA_SugarStepWiseSampler : public StepWiseSamplerSized {
public:
	RNA_SugarStepWiseSampler(
		core::Size const rsd_id,
		PuckerState const pucker_state
	);

	/// @brief Initialization
	void init();

	/// @brief Apply the i-th rotamer to pose
	void apply( core::pose::Pose & pose, core::Size const i );

	/// @brief Get the total number of rotamers in sampler
	core::Size size() const {
		runtime_assert( is_init() );
		return pucker_states_.size();
	}

	/// @brief Set the residue id being sampled
	void set_rsd_id( core::Size const setting ) { rsd_id_ = setting; }

	/// @brief Get the current pucker state.
	PuckerState pucker() const {
		runtime_assert( is_init() );
		return pucker_states_[id()];
	}

	/// @brief Set the pucker_state (WHATEVER / NORTH / SOUTH)
	void set_pucker_state( PuckerState const setting ) {
		set_and_reinit( pucker_state_, setting );
	}

	/// @brief Set if the sampler will skip pucker applying when input pose has
	//  same pucker assginment as sampler.
	void set_skip_same_pucker( bool const setting ) {
		skip_same_pucker_ = setting;
	}

	/// @brief Set if using RNA_IdealCoord to sample puckers
	void set_idealize_coord( bool const setting ) {
		idealize_coord_ = setting;
	}

	/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
	virtual StepWiseSamplerType type() const { return RNA_SUGAR; }

	/// @brief Name of the class
	virtual std::string get_name() const;

private:
	utility::vector1<PuckerState> pucker_states_;

	core::Size rsd_id_;
	PuckerState pucker_state_;

	bool skip_same_pucker_, idealize_coord_;

	std::map < core::id::DOF_ID , core::Real > north_pucker_dof_key_values_;

	std::map < core::id::DOF_ID , core::Real > south_pucker_dof_key_values_;

	bool north_pucker_dofs_have_not_been_initialized_;

	bool south_pucker_dofs_have_not_been_initialized_;
};

} //rna
} //sampler
} //stepwise
} //protocols

#endif
