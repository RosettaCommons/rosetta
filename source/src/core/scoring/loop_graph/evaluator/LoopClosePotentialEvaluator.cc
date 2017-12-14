// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/evaluator/LoopClosePotentialEvaluator.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/loop_graph/evaluator/LoopClosePotentialEvaluator.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.scoring.loop_graph.evaluator.LoopClosePotentialEvaluator" );

namespace core {
namespace scoring {
namespace loop_graph {
namespace evaluator {

//Constructor
LoopClosePotentialEvaluator::LoopClosePotentialEvaluator( core::Real const & loop_fixed_cost ):
	loop_closure_energy_( 0.0 ),
	involves_current_pose_( false ),
	loop_fixed_cost_( loop_fixed_cost )
{}

//Destructor
LoopClosePotentialEvaluator::~LoopClosePotentialEvaluator() = default;

} //evaluator
} //loop_graph
} //scoring
} //core
