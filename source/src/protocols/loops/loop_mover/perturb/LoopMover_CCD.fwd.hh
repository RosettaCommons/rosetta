// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loop_mover/perturb/LoopMover_CCD.fwd.hh
/// @brief  LoopMover_Perutrb_CCD forward declaration
/// @author Brian Weitzner


#ifndef INCLUDED_protocols_loops_loop_mover_perturb_LoopMover_CCD_fwd_hh
#define INCLUDED_protocols_loops_loop_mover_perturb_LoopMover_CCD_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loops {
namespace loop_mover {
namespace perturb {

// Forward
class LoopMover_Perturb_CCD;

typedef utility::pointer::shared_ptr< LoopMover_Perturb_CCD > LoopMover_Perturb_CCDOP;
typedef utility::pointer::shared_ptr< LoopMover_Perturb_CCD const > LoopMover_Perturb_CCDCOP;

} //namespace perturb
} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_perturb_LoopMover_CCD_fwd_hh

