// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/medal/MedalExchangeMover.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_MEDAL_MEDAL_EXCHANGE_MOVER_HH_
#define PROTOCOLS_MEDAL_MEDAL_EXCHANGE_MOVER_HH_

// Unit header
#include <protocols/medal/MedalExchangeMover.fwd.hh>

// C/C++ headers
#include <string>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <protocols/medal/MedalMover.hh>

namespace protocols {
namespace medal {

/// @class 3-stage protocol:
/// 1. Fragment insertion moves with no chainbreak and short-range constraints
/// 2. Fragment insertion moves with no chainbreak and medium-range constraints
/// 3. Alternating fragment insertion and rigid body moves with chainbreak and all constraints
class MedalExchangeMover : public MedalMover {
 public:
  MedalExchangeMover();
  void apply(core::pose::Pose& pose);

  // -- jd2 -- //
  std::string get_name() const;
  protocols::moves::MoverOP clone() const;
  protocols::moves::MoverOP fresh_instance() const;
};

}  // namespace medal
}  // namespace protocols

#endif  // PROTOCOLS_MEDAL_MEDAL_EXCHANGE_MOVER_HH_
