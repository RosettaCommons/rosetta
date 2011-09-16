// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/RigidBodyMotionMover.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/moves/RigidBodyMotionMover.hh>

// C/C++ headers
#include <iostream>
#include <iterator>
#include <string>

// External headers
#include <boost/unordered/unordered_map.hpp>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <numeric/random/random.hh>

// Project headers
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <protocols/loops/Loops.hh>

// Package headers
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace moves {

typedef protocols::loops::Loops Regions;
typedef boost::unordered_map<int, core::kinematics::Jump> Jumps;

static basic::Tracer TR("protocols.medal.RigidBodyMotionMover");

RigidBodyMotionMover::RigidBodyMotionMover(const Regions& regions, const Jumps& jumps)
    : regions_(regions), jumps_(jumps) {}

void RigidBodyMotionMover::apply(core::pose::Pose& pose) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::kinematics::Jump;
  using std::endl;

  // randomly select a jump to manipulate
  Jumps::iterator i = jumps_.begin();
  std::advance(i, numeric::random::random_range(0, jumps_.size() - 1));
  int jump_num = i->first;
  Jump jump = i->second;

  // perturb the rigid body transformation
  jump.gaussian_move(1, option[OptionKeys::rigid::translation](), option[OptionKeys::rigid::rotation]());
  pose.set_jump(jump_num, jump);
}

std::string RigidBodyMotionMover::get_name() const {
  return "RigidBodyMotionMover";
}

MoverOP RigidBodyMotionMover::clone() const {
  return new RigidBodyMotionMover(*this);
}

}  // namespace moves
}  // namespace protocols
