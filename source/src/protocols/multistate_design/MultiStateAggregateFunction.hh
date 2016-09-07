// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MultiStateAggregateFunction.hh
/// @brief
/// @author Colin A. Smith

#ifndef INCLUDED_protocols_multistate_design_MultiStateAggregateFunction_hh
#define INCLUDED_protocols_multistate_design_MultiStateAggregateFunction_hh

#include <protocols/multistate_design/MultiStateFitnessFunction.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.fwd.hh>


namespace protocols {
namespace multistate_design {

class MultiStateAggregateFunction : public utility::pointer::ReferenceCount {

public:
	typedef utility::pointer::shared_ptr< MultiStateAggregateFunction > OP;
	typedef utility::pointer::shared_ptr< MultiStateAggregateFunction const > COP;

	MultiStateAggregateFunction() {}
	~MultiStateAggregateFunction() override = default;

	virtual
	core::Real
	evaluate(
		utility::vector1< core::Real > const & single_state_fitnesses,
		MultiStateFitnessFunction & fitness_function
	) const;

private:
};

} // namespace multistate_design
} // namespace protocols

#endif
