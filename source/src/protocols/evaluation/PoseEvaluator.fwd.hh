// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PoseEvaluator.fwd.hh
/// @brief
/// @details
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_evaluation_PoseEvaluator_fwd_hh
#define INCLUDED_protocols_evaluation_PoseEvaluator_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace evaluation {

class PoseEvaluator;
class MetaPoseEvaluator;
typedef MetaPoseEvaluator PoseEvaluators;

template <class T>
class SingleValuePoseEvaluator;

typedef utility::pointer::shared_ptr< PoseEvaluator > PoseEvaluatorOP;
typedef utility::pointer::shared_ptr< PoseEvaluator const > PoseEvaluatorCOP;

typedef utility::pointer::shared_ptr< MetaPoseEvaluator > MetaPoseEvaluatorOP;
typedef MetaPoseEvaluatorOP PoseEvaluatorsOP;
typedef utility::pointer::shared_ptr< MetaPoseEvaluator const > MetaPoseEvaluatorCOP;
typedef MetaPoseEvaluatorCOP PoseEvaluatorsCOP;

}
}


#endif
