// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LoopCreationMover.fwd.hh
///
/// @brief
/// @author Tim Jacobs


#ifndef INCLUDED_devel_sewing_LoopCreationMover_FWD_HH
#define INCLUDED_devel_sewing_LoopCreationMover_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace loop_creation {

class LoopCreationMover;
typedef utility::pointer::shared_ptr< LoopCreationMover > LoopCreationMoverOP;
typedef utility::pointer::shared_ptr< LoopCreationMover const > LoopCreationMoverCOP;

} //loop creation
} //devel

#endif


