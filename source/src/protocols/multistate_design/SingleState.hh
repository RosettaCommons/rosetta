// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SingleState.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_protocols_multistate_design_SingleState_hh
#define INCLUDED_protocols_multistate_design_SingleState_hh

#include <protocols/multistate_design/SingleState.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/multistate_design/SingleStateFitnessFunction.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>

#ifdef WIN32
#include <core/pose/Pose.hh>
#endif


namespace protocols {
namespace multistate_design {

class SingleState : public utility::pointer::ReferenceCount {

public:
	SingleState();
	virtual ~SingleState();

	SingleState( core::pose::Pose const & pose, bool is_positive );
	SingleState( SingleState const & other );

	virtual void set_best_score( core::Real score ) { best_score_ = score; }
	virtual core::Real best_score() const { return best_score_; }
	virtual bool is_positive_state() const { return is_positive_state_; }

	virtual core::pose::Pose const & pose() const;
	virtual core::pose::Pose & nonconst_pose();

	void fitness_function(SingleStateFitnessFunctionCOP fitness_function) { fitness_function_ = fitness_function; };
	SingleStateFitnessFunctionCOP fitness_function() const { return fitness_function_; }

private:
	core::pose::PoseOP pose_p_;
	bool is_positive_state_;
	core::Real best_score_;
	SingleStateFitnessFunctionCOP fitness_function_;
};

} // namespace multistate_design
} // namespace protocols

#endif
