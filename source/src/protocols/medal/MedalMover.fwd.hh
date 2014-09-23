// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/medal/MedalMover.fwd.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_MEDAL_MEDALMOVER_FWD_HH
#define INCLUDED_PROTOCOLS_MEDAL_MEDALMOVER_FWD_HH

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace medal {

class MedalMover;
typedef utility::pointer::shared_ptr<MedalMover> MedalMoverOP;
typedef utility::pointer::shared_ptr<MedalMover const> MedalMoverCOP;

}  // namespace medal
}  // namespace protocols

#endif  // PROTOCOLS_MEDAL_MEDAL_MOVER_FWD_HH_
