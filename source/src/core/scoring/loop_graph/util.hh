// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/util.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_loop_graph_util_HH
#define INCLUDED_core_scoring_loop_graph_util_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/scoring/loop_graph/util.fwd.hh>
#include <core/scoring/loop_graph/LoopCycle.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/loop_graph/evaluator/LoopClosePotentialEvaluator.fwd.hh>

namespace core {
namespace scoring {
namespace loop_graph {

	void
	get_loop_atom( core::Size const & res,
								 core::pose::Pose const & pose,
								 id::AtomID & atom_id,
								 Vector & xyz,
								 bool const takeoff /* as opposed to landing */ );

	void
	get_loop_atom( core::Size const & res,
								 core::pose::Pose const & pose,
								 std::string const & atom_name,
								 id::AtomID & atom_id,
								 Vector & xyz );

	evaluator::LoopClosePotentialEvaluatorCOP
	get_loop_close_potential( pose::Pose const & pose,
														LoopCycle const & loop_cycle,
														core::Real const & loop_fixed_cost,
														bool const use_6D_potential /* = false */ );

	evaluator::LoopClosePotentialEvaluatorCOP
	get_6D_trans_rot_potential_evaluator( LoopCycle const & loop_cycle,
																				core::Real const & loop_fixed_cost,
																				pose::Pose const & pose );

} //loop_graph
} //scoring
} //core

#endif
