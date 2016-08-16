// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/ScalarWeightedFunc.hh
/// @brief Weighted constraint function that encapsulates other constraints
/// @author James Thompson, Greg Taylor
/// @details The ScalarWeightedFunc is a class used to scale other
/// constraints by multiplying it by a constant value.
/// This is useful in cases where constraint templates are similar
/// to each other and should be down-weighted to avoid double counting
/// for the same structural information.


#ifndef INCLUDED_core_scoring_func_ScalarWeightedFunc_fwd_hh
#define INCLUDED_core_scoring_func_ScalarWeightedFunc_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace func {

class ScalarWeightedFunc;
typedef utility::pointer::shared_ptr< ScalarWeightedFunc > ScalarWeightedFuncOP;
typedef utility::pointer::shared_ptr< ScalarWeightedFunc const > ScalarWeightedFuncCOP;


} // constraints
} // scoring
} // core

#endif
