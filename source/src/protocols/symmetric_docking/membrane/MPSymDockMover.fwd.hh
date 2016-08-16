// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/symmetric_docking/membrane/MPSymDockMover.fwd.hh
///
/// @brief      Membrane symmetric docking protocol
/// @details    Given an asymmetric starting pose, setup the pose for symmetry,
///             add a membrane representation for the full symmetric complex,
///             position the pose at the center of mass of all transmembrane spans,
///             and proceed with the standard symmetric docking protocol.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 2/9/15

#ifndef INCLUDED_protocols_symmetric_docking_MPSymDockMover_fwd_hh
#define INCLUDED_protocols_symmetric_docking_MPSymDockMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace symmetric_docking {
namespace membrane {

class MPSymDockMover;
typedef utility::pointer::shared_ptr< MPSymDockMover > MPSymDockMoverOP;
typedef utility::pointer::shared_ptr< MPSymDockMover const > MPSymDockMoverCOP;

} // membrane
} // symmetric_docking
} // protocols

#endif // INCLUDED_protocols_symmetric_docking_MPSymDockMover_fwd_hh
