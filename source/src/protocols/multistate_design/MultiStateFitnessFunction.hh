// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MultiStateFitnessFunction.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_protocols_multistate_design_MultiStateFitnessFunction_hh
#define INCLUDED_protocols_multistate_design_MultiStateFitnessFunction_hh

#include <protocols/genetic_algorithm/Entity.fwd.hh>
#include <protocols/genetic_algorithm/FitnessFunction.hh>

#include <protocols/multistate_design/SingleState.fwd.hh>
// AUTO-REMOVED #include <protocols/multistate_design/MultiStateEntity.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/multistate_design/MultiStateAggregateFunction.hh>
// AUTO-REMOVED #include <protocols/toolbox/PoseMetricCalculators/MetricValueGetter.hh>

#include <core/types.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>

#include <protocols/multistate_design/MultiStateFitnessFunction.fwd.hh>
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.fwd.hh>
#include <utility/vector1.hh>
#include <map>

#ifdef WIN32
	#include <protocols/multistate_design/SingleState.hh>
	#include <core/pose/Pose.hh>
#endif


namespace protocols {
namespace multistate_design {

class MultiStateFitnessFunction : public genetic_algorithm::FitnessFunction {

public:
	typedef utility::pointer::shared_ptr< MultiStateFitnessFunction > OP;

	MultiStateFitnessFunction();
	virtual ~MultiStateFitnessFunction();

	virtual void add_state( SingleStateOP state );
	virtual void add_state( core::pose::Pose const & pose, bool is_positive );

	virtual core::Real evaluate( protocols::genetic_algorithm::Entity & entity );
	virtual core::Real evaluate( protocols::genetic_algorithm::Entity & entity, core::Size single_state_num ) = 0;
	virtual core::Real evaluate_positive_states( protocols::genetic_algorithm::Entity & entity );

	virtual void set_aggregate_function( MultiStateAggregateFunction::COP aggregate_function );
	virtual MultiStateAggregateFunction::COP aggregate_function() const;

	virtual void set_scorefxn( core::scoring::ScoreFunctionCOP sf );
	virtual core::scoring::ScoreFunctionCOP scorefxn() const;

	///@brief true const (read only) access to states
	virtual SingleStateCOPs const_states( bool positive_only = false ) const;
	virtual SingleStateCOPs positive_states() const;

	virtual core::Size num_states() const { return states_.size(); }
	virtual core::Size num_states( bool pos_neg ) const;
	virtual core::Size num_positive_states() const { return num_states( true ); }
	virtual core::Size num_negative_states() const { return num_states( false ); }

	virtual void add_metric_value_getter(
		std::string const & name,
		protocols::toolbox::pose_metric_calculators::MetricValueGetter const & metric_value_getter
	);

protected:
	typedef std::map<std::string, protocols::toolbox::pose_metric_calculators::MetricValueGetter> MetricValueGetterMap;

	virtual SingleStateOPs & states();
	MetricValueGetterMap const & metric_value_getters() const;

private:
	SingleStateOPs states_;
	MultiStateAggregateFunction::COP aggregate_function_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	// pose(s) to track most fit trait set in real time (such as for graphics)
	utility::vector1< core::pose::PoseOP > best_entity_positive_states_;
	core::Real best_fitness_;
	MetricValueGetterMap metric_value_getters_;
};

} // namespace multistate_design
} // namespace protocols

#endif
