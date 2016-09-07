// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MultiStateFitnessFunction.tmpl.hh
/// @brief
/// @author ashworth

// Unit headers
#include <protocols/multistate_design/MultiStateFitnessFunction.hh>

#include <protocols/multistate_design/SingleState.hh>
#include <protocols/multistate_design/MultiStateEntity.hh>


#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.hh>

#include <protocols/multistate_design/SingleStateEntityData.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace multistate_design {

MultiStateFitnessFunction::MultiStateFitnessFunction()
: genetic_algorithm::FitnessFunction(),
	aggregate_function_(/* NULL */),
	scorefxn_(/* 0 */),
	best_fitness_(0.)
{}

MultiStateFitnessFunction::~MultiStateFitnessFunction()= default;

void
MultiStateFitnessFunction::add_state( core::pose::Pose const & pose, bool is_positive )
{
	add_state( SingleStateOP( new SingleState( pose, is_positive ) ) );
}

void
MultiStateFitnessFunction::add_state( SingleStateOP state )
{
	states_.push_back( state );
	// for real-time pose tracking of best trait set vs. positive state pose(s) (graphics)
	if ( state->is_positive_state() ) {
		core::pose::PoseOP pose( new core::pose::Pose );
		*pose = state->pose();
		// ja this is annoying during iterative protocols because there is no way(?) to close old viewers
		//  protocols::viewer::add_conformation_viewer( pose->conformation(), "Best fitness" );
		best_entity_positive_states_.push_back( pose );
	}
}

core::Real
MultiStateFitnessFunction::evaluate( protocols::genetic_algorithm::Entity & entity )
{
	runtime_assert(aggregate_function_ != nullptr);

	if ( dynamic_cast< protocols::multistate_design::MultiStateEntity * >( &entity ) ) {
		protocols::multistate_design::MultiStateEntity & multi_state_entity =
			static_cast< protocols::multistate_design::MultiStateEntity & >( entity );
		multi_state_entity.single_state_entity_data().resize(states().size());
	}

	utility::vector1<core::Real> single_state_fitnesses(states_.size());
	for ( core::Size i = 1; i <= states_.size(); ++i ) {
		single_state_fitnesses[i] += evaluate(entity, i);
	}

	core::Real fitness = aggregate_function_->evaluate(single_state_fitnesses, *this );
	entity.set_fitness(fitness);

	if ( fitness < best_fitness_ ) {
		best_fitness_ = fitness;
		// real-time pose tracking of best traits vs. positive state pose(s) (graphics)
		auto pose( best_entity_positive_states_.begin() );
		for (auto & state : states_) {
			if ( state->is_positive_state() ) {
				**pose = state->pose();
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
	for (auto & state : states_) {
		if ( !state->is_positive_state() ) continue;
		fitness += evaluate( entity, state ? 1 : 0 ); // FIXME: is this correct? OP -> core::Size
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

/// @brief true const (read only) access to states
SingleStateCOPs
MultiStateFitnessFunction::const_states( bool positive_only /* = false */ ) const
{
	SingleStateCOPs const_states;
	for (const auto & state : states_) {
		if ( positive_only && !state->is_positive_state() ) continue;
		const_states.push_back( state );
	}
	return const_states;
}

SingleStateCOPs
MultiStateFitnessFunction::positive_states() const { return const_states( true ); }

core::Size
MultiStateFitnessFunction::num_states( bool pos_neg ) const
{
	core::Size n(0);
	for (const auto & state : states_) {
		if ( state->is_positive_state() != pos_neg ) continue;
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

