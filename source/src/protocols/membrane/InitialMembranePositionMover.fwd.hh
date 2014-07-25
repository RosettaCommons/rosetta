// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	    protocols/membrane/InitialMembranePositionMover.fwd.hh
///
/// @brief      Compute the Initial Position of the membrane
/// @details	Compute the initial position of the membrane from
///				sequence or structure (can be specified by the user at construction
///				or as a seutp cmd flag)
///
///				NOTE: Requires a membrane pose!
///				NOTE: sequence not yet implemented
///				Last Modified: 6/21/14
///
/// @author		Rebecca Alford (rflaford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_InitialMembranePositionMover_fwd_hh
#define INCLUDED_protocols_membrane_InitialMembranePositionMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {
		
class InitialMembranePositionMover;
typedef utility::pointer::owning_ptr< InitialMembranePositionMover > InitialMembranePositionMoverOP;
typedef utility::pointer::owning_ptr< InitialMembranePositionMover const > InitialMembranePositionMoverCOP;
		
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_InitialMembranePositionMover_fwd_hh
