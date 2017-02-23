// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/evaluator/SixDTransRotPotentialEvaluator.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_loop_graph_SixDTransRotPotentialEvaluator_HH
#define INCLUDED_core_scoring_loop_graph_SixDTransRotPotentialEvaluator_HH

#include <core/scoring/loop_graph/evaluator/LoopClosePotentialEvaluator.hh>
#include <core/scoring/loop_graph/evaluator/SixDTransRotPotentialEvaluator.fwd.hh>
#include <core/scoring/loop_graph/evaluator/SixDTransRotPotential.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>

namespace core {
namespace scoring {
namespace loop_graph {
namespace evaluator {

class SixDTransRotPotentialEvaluator: public LoopClosePotentialEvaluator {

public:

	//constructor
	SixDTransRotPotentialEvaluator( Size const & takeoff_pos,
		Size const & landing_pos,
		pose::Pose const & pose,
		core::Real const & loop_fixed_cost,
		SixDTransRotPotential const & potential );

	//destructor
	~SixDTransRotPotentialEvaluator();

public:

private:

	/// @brief evaluate 6D potential for pose, figuring out atom ids for stubs, etc.
	core::Real
	evaluate( core::pose::Pose const & pose );

	void
	figure_out_if_loop_involves_current_pose( core::pose::Pose const & pose );

	virtual
	void
	get_f1_f2( Vector & f1, Vector & f2, bool const takeoff ) const;

private:

	core::Size const takeoff_pos_, landing_pos_;
	core::id::AtomID takeoff_atom_id_, landing_atom_id_;
	SixDTransRotPotential const & potential_;
	core::kinematics::Stub stub1_;
	core::kinematics::Jump j_;

};

} //evaluator
} //loop_graph
} //scoring
} //core

#endif
