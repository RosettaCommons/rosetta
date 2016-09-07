// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/NumericalDerivCheckResult.hh
/// @brief  Declaration for nuerical derivative check results classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#include <core/optimization/NumericalDerivCheckResult.hh>

#include <utility/vector1.hh>


namespace core {
namespace optimization {

/// @details Auto-generated virtual destructor
SimpleDerivCheckResult::~SimpleDerivCheckResult() = default;

NumericalDerivCheckResult::NumericalDerivCheckResult() : send_to_stdout_( true ) {}
NumericalDerivCheckResult::~NumericalDerivCheckResult() = default;

void
NumericalDerivCheckResult::add_deriv_data(
	NumDerivCheckDataOP deriv_check_data
)
{
	deriv_check_results_.push_back( deriv_check_data );
}

Size
NumericalDerivCheckResult::n_deriv_check_results() const
{
	return deriv_check_results_.size();
}

NumDerivCheckData const &
NumericalDerivCheckResult::deriv_check_result( Size ind ) const
{
	return *deriv_check_results_[ ind ];
}


} // namespace optimization
} // namespace core

