// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/Region.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_CORE_UTIL_CHAINBREAKUTIL_HH
#define INCLUDED_CORE_UTIL_CHAINBREAKUTIL_HH

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

//Auto Headers
#include <utility/vector1.hh>

namespace core {
namespace util {

class ChainbreakUtil {
public:
	/// @brief Returns true if `pose` has a chainbreak, false otherwise
	bool has_chainbreak(const core::pose::Pose& pose) const;

private:
	mutable core::scoring::ScoreFunctionOP score_;
};

}  // namespace util
}  // namespace core

#endif  // CORE_UTIL_CHAINBREAK_UTIL_HH_
