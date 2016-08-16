// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MetricCalculatorFitnessFunction.cc
/// @brief
/// @author Colin A. Smith

#include <protocols/multistate_design/MetricCalculatorFitnessFunction.hh>

#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <basic/MetricValue.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace multistate_design {

MetricCalculatorFitnessFunction::MetricCalculatorFitnessFunction(
	std::string const & calculator_name,
	std::string const & key
) :
	SingleStateFitnessFunction(),
	calculator_name_(calculator_name),
	key_(key)
{
	if ( ! core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( calculator_name_ ) ) {
		basic::Error() << "Tried to tie MetricCalculatorFitnessFunction to PoseMetricCalculator " <<
			calculator_name_ << " but this calculator does not exist." << std::endl;
		utility_exit();
	}
}

core::Real
MetricCalculatorFitnessFunction::evaluate(
	core::pose::Pose const & pose
) const
{
	basic::MetricValue<core::Real> fitness_value;
	pose.metric(calculator_name_, key_, fitness_value);

	return fitness_value.value();
}

} // namespace multistate_design
} // namespace protocols
