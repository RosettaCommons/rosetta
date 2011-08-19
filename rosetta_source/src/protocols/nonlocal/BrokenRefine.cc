// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/BrokenRefine.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/BrokenRefine.hh>

// C/C++ headers
#include <string>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <protocols/filters/Filter.hh>

// Package headers
#include <protocols/nonlocal/BrokenBase.hh>
#include <protocols/nonlocal/Policy.hh>
#include <protocols/nonlocal/PolicyFactory.hh>
#include <protocols/nonlocal/RationalMonteCarlo.hh>
#include <protocols/nonlocal/SingleFragmentMover.hh>
#include <protocols/nonlocal/util.hh>

namespace protocols {
namespace nonlocal {

BrokenRefine::BrokenRefine(core::fragment::FragSetOP fragments,
                           core::kinematics::MoveMapOP movable) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::Size;

  // Consider only the top k best fragments in each library
  Size num_fragments = option[OptionKeys::abinitio::number_3mer_frags];
  Size nres = fragments->max_pos();

  // fragment selection policies
  PolicyOP uniform = PolicyFactory::get_policy("uniform", fragments, num_fragments);
  PolicyOP smooth = PolicyFactory::get_policy("smooth", fragments, num_fragments);

  for (Size i = 1; i <= 8; ++i) {
    PolicyOP policy = (i <= 5) ? uniform : smooth;
    add_mover(make_minimizer(score_function(i, nres)));
    add_mover(new RationalMonteCarlo(
        new SingleFragmentMover(fragments, movable, policy),
        score_function(i, nres),
        num_cycles(i),
        option[OptionKeys::abinitio::temperature](),
        true));
  }
}

core::Size BrokenRefine::num_cycles(int stage) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  assert(stage >= 0);

  double m = option[OptionKeys::abinitio::increase_cycles]();
  double cycles = (stage <= 5) ? 2000 * m : 4000 * m;
  return static_cast<core::Size>(cycles);
}

std::string BrokenRefine::get_name() const {
  return "BrokenRefine";
}

}  // namespace nonlocal
}  // namespace protocols
