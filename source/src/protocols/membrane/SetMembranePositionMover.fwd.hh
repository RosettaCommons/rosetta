// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/membrane/SetMembranePositionMover.fwd.hh
///
/// @brief  Sets the membrane position normal and center
/// @details Sets the membrane position normal and center
///    CAUTION: ONLY FOR FLEXIBLE MEMBRANE AND FIXED PROTEIN!!!
///    Last Modified: 6/28/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_SetMembranePositionMover_fwd_hh
#define INCLUDED_protocols_membrane_SetMembranePositionMover_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {

/// @brief Membrane Position Rotation/Translation Mover (by stub & jump)
class SetMembranePositionMover;
typedef utility::pointer::shared_ptr< SetMembranePositionMover > SetMembranePositionMoverOP;
typedef utility::pointer::shared_ptr< SetMembranePositionMover const > SetMembranePositionMoverCOP;

/// @brief Membrane Position Rotation Mover (by jump)
class SetMembraneNormalMover;
typedef utility::pointer::shared_ptr< SetMembraneNormalMover > SetMembraneNormalMoverOP;
typedef utility::pointer::shared_ptr< SetMembraneNormalMover const > SetMembraneNormalMoverCOP;

/// @brief Memrbane Position Translaiton Mover (by stub)
class SetMembraneCenterMover;
typedef utility::pointer::shared_ptr< SetMembraneCenterMover > SetMembraneCenterMoverOP;
typedef utility::pointer::shared_ptr< SetMembraneCenterMover const > SetMembraneCenterMoverCOP;

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_SetMembranePositionMover_fwd_hh
