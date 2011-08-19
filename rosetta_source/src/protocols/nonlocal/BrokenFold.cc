// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/BrokenFold.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/BrokenFold.hh>

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
#include <protocols/filters/Filter.hh>

// Package headers
#include <protocols/nonlocal/Policy.hh>
#include <protocols/nonlocal/PolicyFactory.hh>
#include <protocols/nonlocal/RationalMonteCarlo.hh>
#include <protocols/nonlocal/SingleFragmentMover.hh>
#include <protocols/nonlocal/util.hh>

namespace protocols {
namespace nonlocal {

BrokenFold::BrokenFold(core::fragment::FragSetOP fragments_lg,
                       core::fragment::FragSetOP fragments_sm,
                       core::kinematics::MoveMapOP movable) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::Size;
  using protocols::moves::MoverOP;

  // Consider only the top k best fragments in each library
  Size num_small = option[OptionKeys::abinitio::number_3mer_frags];
  Size num_large = option[OptionKeys::abinitio::number_9mer_frags];

  // fragment selection policies
  PolicyOP uniform_lg = PolicyFactory::get_policy("uniform", fragments_lg, num_large);
  PolicyOP uniform_sm = PolicyFactory::get_policy("uniform", fragments_sm, num_small);
  PolicyOP smooth_sm  = PolicyFactory::get_policy("smooth",  fragments_sm, num_small);

  assert(fragments_lg->max_pos() == fragments_sm->max_pos());
  Size nres = fragments_sm->max_pos();

  // large fragments, uniformly chosen
  for (Size i = 1; i <= 3; ++i) {
    MoverOP mover = new RationalMonteCarlo(
        new SingleFragmentMover(fragments_lg, movable, uniform_lg),
        score_function(i, nres),
        num_cycles(i),
        temperature(i),
        true);

    add_mover(mover);
  }

  // small fragments, uniformly chosen
  for (Size i = 4; i <= 6; ++i) {
    MoverOP mover = new RationalMonteCarlo(
        new SingleFragmentMover(fragments_sm, movable, uniform_sm),
        score_function(i, nres),
        num_cycles(i),
        temperature(i),
        true);

    add_mover(mover);
  }

  // small fragments, smoothly chosen
  for (Size i = 7; i <= 9; ++i) {
    MoverOP mover = new RationalMonteCarlo(
        new SingleFragmentMover(fragments_sm, movable, smooth_sm),
        score_function(i, nres),
        num_cycles(i),
        temperature(i),
        true);

    add_mover(mover);
  }
}

std::string BrokenFold::get_name() const {
  return "BrokenFold";
}

core::Size BrokenFold::num_cycles(int i) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  core::Size amount;
  switch (i) {
    case 1:
    case 2:
    case 3:
      amount = 5000;
      break;
    case 4:
    case 5:
    case 6:
      amount = 7500;
      break;
    case 7:
    case 8:
    case 9:
      amount = 10000;
      break;
  }
  return static_cast<core::Size>(
      amount * option[OptionKeys::abinitio::increase_cycles]());
}

core::Real BrokenFold::temperature(int i) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  core::Real factor;
  switch (i) {
    case 1:
    case 2:
    case 3:
      factor = 1.0;
    case 4:
    case 5:
    case 6:
      factor = 0.8;
    case 7:
    case 8:
    case 9:
      factor = 0.6;
  }
  return factor * option[OptionKeys::abinitio::temperature]();
}

}  // namespace nonlocal
}  // namespace protocols
