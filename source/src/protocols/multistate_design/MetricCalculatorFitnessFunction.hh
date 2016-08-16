// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MetricCalculatorFitnessFunction.hh
/// @brief
/// @author Colin A. Smith

#ifndef INCLUDED_protocols_multistate_design_MetricCalculatorFitnessFunction_hh
#define INCLUDED_protocols_multistate_design_MetricCalculatorFitnessFunction_hh

#include <protocols/multistate_design/MetricCalculatorFitnessFunction.fwd.hh>
#include <protocols/multistate_design/SingleStateFitnessFunction.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace multistate_design {

class MetricCalculatorFitnessFunction : public SingleStateFitnessFunction {

public:
	MetricCalculatorFitnessFunction() : SingleStateFitnessFunction() {}
	MetricCalculatorFitnessFunction(
		std::string const & calculator_name,
		std::string const & key
	);

	virtual ~MetricCalculatorFitnessFunction() {}

	virtual core::Real evaluate( core::pose::Pose const & pose ) const;

private:
	std::string calculator_name_;
	std::string key_;
};

} // namespace multistate_design
} // namespace protocols

#endif
