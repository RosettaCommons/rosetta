// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/util/ChainbreakUtil.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <core/util/ChainbreakUtil.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/util.hh>

//Auto Headers
#include <utility/vector1.hh>

namespace core {
namespace util {

// Static member initialization
core::scoring::ScoreFunctionOP ChainbreakUtil::score_ = 0;

bool ChainbreakUtil::has_chainbreak(const core::pose::Pose& pose) {
  using core::pose::Pose;

  if (!score_) {
    score_ = new core::scoring::ScoreFunction();
		score_->set_weight( core::scoring::linear_chainbreak, 1.0 );
  }

  Pose copy(pose);
  core::pose::correctly_add_cutpoint_variants(copy);

  score_->score(copy);
  return (*score_)[core::scoring::linear_chainbreak] > 0;
}

}  // namesapce util
}  // namespace core
