// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/methods/chainbreak_eval.hh
/// @brief  methods for chainbreak evaluation
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_methods_chainbreak_eval_hh
#define INCLUDED_protocols_forge_methods_chainbreak_eval_hh

// type headers
#include <core/types.hh>

// project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace methods {


/// @brief evaluate linear chainbreak at a position
/// @remarks Copies the Pose, if necessary swaps the Pose with a cut fold tree,
///  then evaluates the chainbreak.
core::Real
linear_chainbreak(
	core::pose::Pose const & pose,
	core::Size const pos
);


/// @brief evaluate linear chainbreak at a position
/// @remarks If necessary, will evaluate using a copy of the Pose with a cut
///  fold tree.
core::Real
linear_chainbreak(
	core::pose::Pose & pose,
	core::Size const pos
);


/// @brief evaluate overlap chainbreak at a position
/// @remarks Copies the Pose, if necessary swaps the Pose with a cut fold tree,
///  then evaluates the chainbreak.
core::Real
overlap_chainbreak(
	core::pose::Pose const & pose,
	core::Size const pos
);


/// @brief evaluate overlap chainbreak at a position
/// @remarks If necessary, will evaluate using a copy of the Pose with a cut
///  fold tree.
core::Real
overlap_chainbreak(
	core::pose::Pose & pose,
	core::Size const pos
);


/// @brief evaluate quadratic chainbreak at a position
/// @remarks Copies the Pose, if necessary swaps the Pose with a cut fold tree,
///  then evaluates the chainbreak.
core::Real
quadratic_chainbreak(
	core::pose::Pose const & pose,
	core::Size const pos
);


/// @brief evaluate quadratic chainbreak at a position
/// @remarks If necessary, will evaluate using a copy of the Pose with a cut
///  fold tree.
core::Real
quadratic_chainbreak(
	core::pose::Pose & pose,
	core::Size const pos
);


} // methods
} // forge
} // protocols


#endif /* INCLUDED_protocols_forge_methods_chainbreak_eval_HH */
