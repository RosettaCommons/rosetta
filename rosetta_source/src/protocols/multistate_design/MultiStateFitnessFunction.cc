// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MultiStateFitnessFunction.tmpl.hh
/// @brief
/// @author ashworth

// Unit headers
#include <protocols/multistate_design/MultiStateFitnessFunction.hh>

#include <protocols/multistate_design/SingleState.hh>
#include <protocols/multistate_design/MultiStateEntity.hh>

#include <protocols/viewer/viewers.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <basic/MetricValue.fwd.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <protocols/genetic_algorithm/Entity.hh>
#include <protocols/genetic_algorithm/FitnessFunction.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/multistate_design/MultiStateAggregateFunction.hh>
#include <protocols/multistate_design/MultiStateFitnessFunction.fwd.hh>
#include <protocols/multistate_design/SingleState.fwd.hh>
#include <protocols/multistate_design/SingleStateEntityData.fwd.hh>
#include <protocols/multistate_design/SingleStateEntityData.hh>
#include <protocols/multistate_design/SingleStateFitnessFunction.fwd.hh>
#include <protocols/multistate_design/SingleStateFitnessFunction.hh>
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <ObjexxFCL/format.hh>
#include <algorithm>
#include <cassert>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

namespace protocols {
namespace multistate_design {

MultiStateFitnessFunction::MultiStateFitnessFunction()
	: genetic_algorithm::FitnessFunction(),
		aggregate_function_(NULL),
		scorefxn_(0),
		best_fitness_(0.)
{}

MultiStateFitnessFunction::~MultiStateFitnessFunction(){}

void
MultiStateFitnessFunction::add_state( core::pose::Pose const & pose, bool is_positive )
{
	add_state( new SingleState( pose, is_positive ) );
}

void
MultiStateFitnessFunction::add_state( SingleStateOP state )
{
	states_.push_back( state );
	// for real-time pose tracking of best trait set vs. positive state pose(s) (graphics)
	if ( state->is_positive_state() ) {
		core::pose::PoseOP pose = new core::pose::Pose;
		*pose = state->pose();
// ja this is annoying during iterative protocols because there is no way(?) to close old viewers
//		protocols::viewer::add_conformation_viewer( pose->conformation(), "Best fitness" );
		best_entity_positive_states_.push_back( pose );
	}
}

core::Real
MultiStateFitnessFunction::evaluate( protocols::genetic_algorithm::Entity & entity )
{
	runtime_assert(aggregate_function_);

	protocols::multistate_design::MultiStateEntity* multi_state_entity =
		dynamic_cast< protocols::multistate_design::MultiStateEntity* >( &entity );
	if (multi_state_entity) multi_state_entity->single_state_entity_data().resize(states().size());

	utility::vector1<core::Real> single_state_fitnesses(states_.size());
	for (core::Size i = 1; i <= states_.size(); ++i) {
		single_state_fitnesses[i] += evaluate(entity, i);
	}

	core::Real fitness = aggregate_function_->evaluate(single_state_fitnesses, *this );
	entity.set_fitness(fitness);

	if ( fitness < best_fitness_ ) {
		best_fitness_ = fitness;
		// real-time pose tracking of best traits vs. positive state pose(s) (graphics)
		utility::vector1< core::pose::PoseOP >::iterator pose( best_entity_positive_states_.begin() );
		for ( SingleStateOPs::iterator s( states_.begin() ), end( states_.end() ); s != end; ++s ) {
			if ( (*s)->is_positive_state() ) {
				**pose = (*s)->pose();
				++pose;
			}
		}
	}

	return fitness;
}

core::Real
MultiStateFitnessFunction::evaluate_positive_states( protocols::genetic_algorithm::Entity & entity )
{
	core::Real fitness(0.);
	for ( SingleStateOPs::iterator s( states_.begin() ), end( states_.end() ); s != end; ++s ) {
		if ( !(*s)->is_positive_state() ) continue;
		fitness += evaluate( entity, *s );
	}
	return fitness;
}

void
MultiStateFitnessFunction::set_aggregate_function( MultiStateAggregateFunction::COP aggregate_function )
{
	aggregate_function_ = aggregate_function;
}

MultiStateAggregateFunction::COP
MultiStateFitnessFunction::aggregate_function() const
{
	return aggregate_function_;
}

void
MultiStateFitnessFunction::set_scorefxn( core::scoring::ScoreFunctionCOP sf ) { scorefxn_ = sf; }

core::scoring::ScoreFunctionCOP
MultiStateFitnessFunction::scorefxn() const { return scorefxn_; }

///@brief true const (read only) access to states
SingleStateCOPs
MultiStateFitnessFunction::const_states( bool positive_only /* = false */ ) const
{
	SingleStateCOPs const_states;
	for ( SingleStateOPs::const_iterator s( states_.begin() ), end( states_.end() );
				s != end; ++s ) {
		if ( positive_only && !(*s)->is_positive_state() ) continue;
		const_states.push_back( *s );
	}
	return const_states;
}

SingleStateCOPs
MultiStateFitnessFunction::positive_states() const { return const_states( true ); }

core::Size
MultiStateFitnessFunction::num_states( bool pos_neg ) const
{
	core::Size n(0);
	for ( SingleStateOPs::const_iterator s( states_.begin() ), end( states_.end() );
				s != end; ++s ) {
		if ( (*s)->is_positive_state() != pos_neg ) continue;
		++n;
	}
	return n;
}

void
MultiStateFitnessFunction::add_metric_value_getter(
	std::string const & name,
	protocols::toolbox::pose_metric_calculators::MetricValueGetter const & metric_value_getter
)
{
	metric_value_getters_[name] = metric_value_getter;
}

SingleStateOPs &
MultiStateFitnessFunction::states() { return states_; }

std::map< std::string, protocols::toolbox::pose_metric_calculators::MetricValueGetter > const &
MultiStateFitnessFunction::metric_value_getters() const { return metric_value_getters_; }

} // namespace multistate_design
} // namespace protocols

