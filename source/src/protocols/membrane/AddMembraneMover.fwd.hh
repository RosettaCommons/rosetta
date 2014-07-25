// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/AddMembraneMover.fwd.hh
///
/// @brief      Add Membrane Representation to the Pose
/// @details	Given a pose, setup membrane topology, lips info,
///				and a membrane virtual residue in the pose. All of this information
///				is coordinated via the MembraneInfo object maintained in
///				the Pose's conformation. After applying AddMembraneMover
///				to the pose, pose.conformation().is_membrane() should always
///				return true.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (6/21/14)

#ifndef INCLUDED_protocols_membrane_AddMembraneMover_fwd_hh
#define INCLUDED_protocols_membrane_AddMembraneMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {

class AddMembraneMover;
typedef utility::pointer::owning_ptr< AddMembraneMover > AddMembraneMoverOP;
typedef utility::pointer::owning_ptr< AddMembraneMover const > AddMembraneMoverCOP;

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMembraneMover_fwd_hh
