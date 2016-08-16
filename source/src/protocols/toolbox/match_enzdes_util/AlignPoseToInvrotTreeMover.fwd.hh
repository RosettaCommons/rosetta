// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/match_enzdes_util/AlignPoseToInvrotTreeMover.fwd.hh
/// @brief  Forward declaration for inverse rotamer tree node base
/// @author Florian Richter, flosopher@gmail.com, mar 2012

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_AlignPoseToInvrotTreeMover_fwd_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_AlignPoseToInvrotTreeMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

class AlignPoseToInvrotTreeMover;

typedef utility::pointer::shared_ptr< AlignPoseToInvrotTreeMover > AlignPoseToInvrotTreeMoverOP;
typedef utility::pointer::shared_ptr< AlignPoseToInvrotTreeMover const > AlignPoseToInvrotTreeMoverCOP;


}
}
}

#endif
