// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/loop_graph/LoopCloseEnergy.cc
/// @brief  Cost of bringing two chains together.
/// @author Rhiju Das


// Unit headers
#include <core/scoring/loop_graph/LoopCloseEnergy.hh>
#include <core/scoring/loop_graph/LoopCloseEnergyCreator.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/loop_graph/LoopGraph.hh>
#include <core/scoring/loop_graph/evaluator/LoopClosePotentialEvaluator.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.io.hh>

// C++
static THREAD_LOCAL basic::Tracer TR( "core.scoring.loop_graph.LoopCloseEnergy" );

namespace core {
namespace scoring {
namespace loop_graph {


/// @details This must return a fresh instance of the LoopCloseEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
LoopCloseEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new LoopCloseEnergy );
}

ScoreTypes
LoopCloseEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( loop_close );
	return sts;
}


/// c-tor
LoopCloseEnergy::LoopCloseEnergy() :
	parent( methods::EnergyMethodCreatorOP( new LoopCloseEnergyCreator ) )
{}

/// clone
methods::EnergyMethodOP
LoopCloseEnergy::clone() const
{
	return methods::EnergyMethodOP( new LoopCloseEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
void
LoopCloseEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	update_loop_atoms_and_lengths( pose );
}

/////////////////////////////////////////////////////////////////////////////
void
LoopCloseEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	update_loop_atoms_and_lengths( pose );
}

/////////////////////////////////////////////////////////////////////////////
void
LoopCloseEnergy::update_loop_atoms_and_lengths( pose::Pose & pose ) const {
	if ( !loop_graph_ ) loop_graph_ = core::scoring::loop_graph::LoopGraphOP( new scoring::loop_graph::LoopGraph );
	loop_graph_->update( pose );
}

///////////////////////////////////////////////////////////////////////////////
void
LoopCloseEnergy::finalize_total_energy(
	pose::Pose &,
	ScoreFunction const &,
	EnergyMap & totals
) const {

	using namespace core::scoring::loop_graph;
	using namespace core::scoring::constraints;

	totals[ loop_close          ] = loop_graph_->total_energy();

} // finalize_total_energy

///////////////////////////////////////////////////////////////////////////////
void
LoopCloseEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const &,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{

	using namespace conformation;
	using namespace core::scoring::loop_graph;

	Vector f1, f2;

	for ( Size n = 1; n <= loop_graph_->num_current_pose_loops(); n++ ) {

		evaluator::LoopClosePotentialEvaluatorCOP loop_score_evaluator = loop_graph_->loop_score_evaluator( n );

		if ( atom_id == loop_score_evaluator->current_pose_takeoff_atom() ) {
			loop_score_evaluator->get_f1_f2( f1, f2, true /* takeoff */ );
		  F1 += weights[ loop_close ] * f1;
			F2 += weights[ loop_close ] * f2;
		}

		if ( atom_id == loop_score_evaluator->current_pose_landing_atom() ) {
			loop_score_evaluator->get_f1_f2( f1, f2, false /* takeoff */ );
			F1 += weights[ loop_close ] * f1;
			F2 += weights[ loop_close ] * f2;
		}

	}
} // eval atom derivative

core::Size
LoopCloseEnergy::version() const
{
	return 2;
}


} // loop_graph
} // scoring
} // core
