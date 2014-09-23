// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/medal/MedalFreeMover.fwd.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_MEDAL_MEDALEXCHANGEMOVER_FWD_HH
#define INCLUDED_PROTOCOLS_MEDAL_MEDALEXCHANGEMOVER_FWD_HH

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace medal {

class MedalFreeMover;
typedef utility::pointer::shared_ptr<MedalFreeMover> MedalFreeMoverOP;
typedef utility::pointer::shared_ptr<MedalFreeMover const> MedalFreeMoverCOP;

}  // namespace medal
}  // namespace protocols

#endif  // PROTOCOLS_MEDAL_MEDAL_FREE_MOVER_FWD_HH_
