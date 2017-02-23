// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/evaluator/GaussianChainFuncPotentialEvaluator.fwd.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_loop_graph_GaussianChainFuncPotentialEvaluator_FWD_HH
#define INCLUDED_core_scoring_loop_graph_GaussianChainFuncPotentialEvaluator_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace loop_graph {
namespace evaluator {

class GaussianChainFuncPotentialEvaluator;
typedef utility::pointer::shared_ptr< GaussianChainFuncPotentialEvaluator > GaussianChainFuncPotentialEvaluatorOP;
typedef utility::pointer::shared_ptr< GaussianChainFuncPotentialEvaluator const > GaussianChainFuncPotentialEvaluatorCOP;

} //evaluator
} //loop_graph
} //scoring
} //core

#endif
